import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import statsmodels.api as sm
from scipy import stats
from scipy.stats import chi2
import glob
import warnings


class PDIAnalyzer:
    """Analyzes accuracy vs PDI using multiple model fits."""

    def __init__(self, data_path="../Data/*.csv"):
        self.data = self._load_data(data_path)
        self.df_cell = self._prepare_data()
        self.models = {}
        self._predictor_map = {
            'Linear': ['PDI100'],
            'Quadratic': ['PDI100', 'PDI1002'],
            'Cubic': ['PDI100', 'PDI1002', 'PDI1003'],
            'Log': ['lPDIc']
        }

    def _load_data(self, path):
        """Load and combine CSV files."""
        files = glob.glob(path)
        if len(files) == 1:
            return pd.read_csv(files[0])

        data_list = []
        for i, file in enumerate(files):
            df = pd.read_csv(file)
            df['source'] = i
            data_list.append(df)
        return pd.concat(data_list, ignore_index=True)

    def _prepare_data(self):
        """Clean and aggregate data."""
        # Filter and transform
        df = (self.data
              .query("Trial_type == 'experiment'")
              [['PART_ID', 'CSI', 'Corr']]
              .rename(columns={'CSI': 'PDI', 'Corr': 'CORR'})
              .assign(
            PDI=lambda x: pd.to_numeric(x['PDI']) * 16.67,
            PART_ID=lambda x: x['PART_ID'].astype('category'),
            CORR=lambda x: x['CORR'].astype(int)
        )
              .query("PDI >= 0"))

        # Aggregate by participant and PDI
        df_cell = (df.groupby(['PART_ID', 'PDI'])
                   .agg(correct_n=('CORR', 'sum'), total_n=('CORR', 'count'))
                   .reset_index())

        # Sanity check
        if (df_cell['total_n'] != 11).any():
            warnings.warn("Some cells don't have 11 trials.")

        # Create predictors
        pdi_mean = df_cell['PDI'].mean()
        log_pdi_mean = np.log1p(df_cell['PDI']).mean()

        return df_cell.assign(
            PDI100=(df_cell['PDI'] - pdi_mean) / 100,
            PDI1002=lambda x: x['PDI100'] ** 2,
            PDI1003=lambda x: x['PDI100'] ** 3,
            lPDIc=np.log1p(df_cell['PDI']) - log_pdi_mean,
            prop_correct=lambda x: x['correct_n'] / x['total_n']
        )

    def fit_models(self):
        """Fit all four models."""
        model_specs = {
            'Linear': ['PDI100'],
            'Quadratic': ['PDI100', 'PDI1002'],
            'Cubic': ['PDI100', 'PDI1002', 'PDI1003'],
            'Log': ['lPDIc']
        }

        for name, predictors in model_specs.items():
            self.models[name] = self._fit_binomial_model(predictors)

        return self.models

    def _fit_binomial_model(self, predictors):
        """Fit binomial GLM with GEE for clustering."""
        try:
            X = sm.add_constant(self.df_cell[predictors])
            model = sm.GEE(
                self.df_cell['prop_correct'], X,
                groups=self.df_cell['PART_ID'],
                family=sm.families.Binomial(),
                cov_struct=sm.cov_struct.Exchangeable()
            )
            return model.fit()
        except:
            # Fallback to regular GLM
            model = sm.GLM(self.df_cell['prop_correct'], X,
                           family=sm.families.Binomial())
            return model.fit()

    def analyze_trends(self):
        """Perform statistical tests and comparisons."""
        if 'Linear' not in self.models:
            raise ValueError("Models must be fitted first.")

        # Linear trend analysis
        lin_model = self.models['Linear']
        b_lin = lin_model.params['PDI100']
        se_lin = lin_model.bse['PDI100']
        z_lin = b_lin / se_lin
        p_lin = 2 * (1 - stats.norm.cdf(abs(z_lin)))

        print("=== Linear Trend Test ===")
        print(f"Slope per +100ms: b={b_lin:.4f}, SE={se_lin:.4f}, z={z_lin:.2f}, p={p_lin:.4g}")
        print(
            f"Odds ratio: {np.exp(b_lin):.3f} [{np.exp(b_lin - 1.96 * se_lin):.3f}, {np.exp(b_lin + 1.96 * se_lin):.3f}]")

        # Polynomial term tests
        self._test_polynomial_terms()

        # AIC comparison
        self._compare_models()

    def _test_polynomial_terms(self):
        """Test significance of higher-order terms."""
        print("\n=== Polynomial Terms ===")

        for model_name, param in [('Quadratic', 'PDI1002'), ('Cubic', 'PDI1003')]:
            if model_name in self.models and param in self.models[model_name].params:
                model = self.models[model_name]
                b = model.params[param]
                se = model.bse[param]
                z = b / se
                p = 2 * (1 - stats.norm.cdf(abs(z)))
                print(f"{param}: b={b:.4f}, SE={se:.4f}, z={z:.2f}, p={p:.4g}")

    def _weighted_r2(self, y_true, y_pred, weights=None):
        """Compute weighted R^2 (coefficient of determination) on probabilities.

        Parameters
        - y_true: array-like of true proportions
        - y_pred: array-like of predicted proportions
        - weights: array-like of weights (e.g., total_n per cell); if None, equal weights are used
        """
        y_true = np.asarray(y_true, dtype=float)
        y_pred = np.asarray(y_pred, dtype=float)
        if weights is None:
            weights = np.ones_like(y_true, dtype=float)
        else:
            weights = np.asarray(weights, dtype=float)
        m = np.isfinite(y_true) & np.isfinite(y_pred) & np.isfinite(weights)
        if not np.any(m):
            return np.nan
        y_true = y_true[m]
        y_pred = y_pred[m]
        weights = weights[m]
        w_sum = weights.sum()
        if w_sum <= 0:
            return np.nan
        y_bar = np.sum(weights * y_true) / w_sum
        ss_tot = np.sum(weights * (y_true - y_bar) ** 2)
        ss_res = np.sum(weights * (y_true - y_pred) ** 2)
        if ss_tot <= 0:
            return np.nan
        return 1.0 - ss_res / ss_tot

    def _compare_models(self):
        """Compare models using AIC and weighted R^2 (on fitted probabilities)."""
        print("\n=== Model Comparison (AIC, R^2) ===")

        comp = self._rank_models()

        # Print table
        for _, row in comp.iterrows():
            aic_str = f"{row['AIC']:.1f}" if np.isfinite(row['AIC']) else "NA"
            daic_str = f", ΔAIC={row['ΔAIC']:.1f}" if 'ΔAIC' in comp.columns and np.isfinite(row.get('ΔAIC', np.nan)) else ""
            r2_str = f"{row['R2']:.3f}" if np.isfinite(row['R2']) else "NA"
            print(f"{row['Model']:<10} AIC={aic_str}{daic_str}, R^2={r2_str}")

        # Best model message
        if 'ΔAIC' in comp.columns:
            print(f"\nBest model (AIC): {comp.iloc[0]['Model']}")
        else:
            print(f"\nBest model (R^2): {comp.iloc[0]['Model']}")

    def _rank_models(self):
        """Create a DataFrame ranking models by AIC (if available) else R^2.
        Returns a DataFrame sorted by criterion (best first)."""
        # Prepare true values and weights
        y_true = self.df_cell['prop_correct'].values
        weights = self.df_cell['total_n'].values if 'total_n' in self.df_cell.columns else np.ones_like(y_true)

        rows = []
        for name, model in self.models.items():
            # Predictions and R^2
            r2 = np.nan
            try:
                cols = self._predictor_map.get(name, [])
                X = sm.add_constant(self.df_cell[cols])
                y_pred = np.asarray(model.predict(X), dtype=float)
                r2 = self._weighted_r2(y_true, y_pred, weights)
            except Exception:
                r2 = np.nan

            # AIC if available
            aic = getattr(model, 'aic', np.nan)

            rows.append({'Model': name, 'AIC': aic, 'R2': r2})

        comp = pd.DataFrame(rows)

        # Sort by AIC if any available, otherwise by R^2 descending
        if np.isfinite(comp['AIC']).any():
            comp = comp.sort_values('AIC')
            comp['ΔAIC'] = comp['AIC'] - np.nanmin(comp['AIC'])
        else:
            comp = comp.sort_values('R2', ascending=False)
        comp = comp.reset_index(drop=True)
        return comp

    def plot_fits(self, figsize=(10, 6)):
        """Plot all model fits with confidence intervals."""
        if not self.models:
            raise ValueError("Models must be fitted first.")

        # Create prediction grid
        pdi_range = np.linspace(self.df_cell['PDI'].min(), self.df_cell['PDI'].max(), 200)
        grid = self._create_prediction_grid(pdi_range)

        # Plot setup
        fig, ax = plt.subplots(figsize=figsize)
        colors = {'Linear': '#2c7fb8', 'Quadratic': '#f03b20',
                  'Cubic': '#6a3d9a', 'Log': '#33a02c'}

        # Plot each model
        for name, model in self.models.items():
            predictions = self._predict_with_ci(model, grid, name)

            ax.fill_between(pdi_range, predictions['lo'], predictions['hi'],
                            alpha=0.15, color=colors[name])
            ax.plot(pdi_range, predictions['fit'],
                    color=colors[name], linewidth=2, label=name)

        ax.set_xlabel('PDI (ms)')
        ax.set_ylabel('Accuracy')
        ax.set_title('Model Fits: Accuracy vs PDI (95% CI)')
        ax.legend()
        ax.grid(True, alpha=0.3)
        plt.tight_layout()
        plt.show()

    def _create_prediction_grid(self, pdi_range):
        """Create prediction grid with all predictors."""
        pdi_mean = self.df_cell['PDI'].mean()
        log_pdi_mean = np.log1p(self.df_cell['PDI']).mean()

        pdi100 = (pdi_range - pdi_mean) / 100
        return pd.DataFrame({
            'PDI': pdi_range,
            'PDI100': pdi100,
            'PDI1002': pdi100 ** 2,
            'PDI1003': pdi100 ** 3,
            'lPDIc': np.log1p(pdi_range) - log_pdi_mean
        })

    def _predict_with_ci(self, model, grid, model_name):
        """Make predictions with confidence intervals."""
        predictor_map = {
            'Linear': ['PDI100'],
            'Quadratic': ['PDI100', 'PDI1002'],
            'Cubic': ['PDI100', 'PDI1002', 'PDI1003'],
            'Log': ['lPDIc']
        }

        try:
            X_pred = sm.add_constant(grid[predictor_map[model_name]])
            pred = model.predict(X_pred)
            se_pred = np.sqrt(np.diag(X_pred @ model.cov_params() @ X_pred.T))

            # Convert to probability scale
            fit = 1 / (1 + np.exp(-pred))
            lo = 1 / (1 + np.exp(-(pred - 1.96 * se_pred)))
            hi = 1 / (1 + np.exp(-(pred + 1.96 * se_pred)))

            return {'fit': fit, 'lo': lo, 'hi': hi}
        except:
            # Fallback
            n = len(grid)
            return {'fit': np.full(n, 0.5), 'lo': np.full(n, 0.4), 'hi': np.full(n, 0.6)}

    def _hosmer_lemeshow(self, mu, n, y_count, n_groups=10):
        """Hosmer–Lemeshow goodness-of-fit test for grouped/binomial data.

        Parameters
        - mu: predicted probabilities (array-like)
        - n: number of trials per observation (array-like)
        - y_count: observed successes (array-like)
        - n_groups: number of groups (deciles by default)
        Returns: dict with chi2, df, p, table (per-group summary)
        """
        mu = np.asarray(mu, dtype=float)
        n = np.asarray(n, dtype=float)
        y_count = np.asarray(y_count, dtype=float)
        order = np.argsort(mu)
        mu, n, y_count = mu[order], n[order], y_count[order]
        k = len(mu)
        if k < n_groups:
            n_groups = max(2, k)
        # Define group boundaries (approximately equal size)
        edges = np.linspace(0, k, n_groups + 1).astype(int)
        stats_rows = []
        hl = 0.0
        for g in range(n_groups):
            s, e = edges[g], edges[g + 1]
            if e <= s:
                continue
            ng = n[s:e]
            mug = mu[s:e]
            yg = y_count[s:e]
            Og = np.sum(yg)
            Eg = np.sum(ng * mug)
            Og0 = np.sum(ng - yg)
            Eg0 = np.sum(ng * (1.0 - mug))
            # Guard against division by zero
            if Eg > 0:
                term1 = (Og - Eg) ** 2 / Eg
            else:
                term1 = 0.0
            if Eg0 > 0:
                term2 = (Og0 - Eg0) ** 2 / Eg0
            else:
                term2 = 0.0
            hl += term1 + term2
            stats_rows.append({
                'group': g + 1,
                'n_obs': int(np.sum(ng)),
                'obs_succ': float(Og),
                'exp_succ': float(Eg),
                'obs_fail': float(Og0),
                'exp_fail': float(Eg0),
                'mean_mu': float(np.average(mug, weights=ng))
            })
        df = max(1, len(stats_rows) - 2)
        p = 1.0 - chi2.cdf(hl, df)
        return {'chi2': float(hl), 'df': int(df), 'p': float(p), 'table': pd.DataFrame(stats_rows)}

    def residual_analysis(self, show_plots=True, save_prefix=None, n_groups=10, figsize=(12, 10)):
        """Residual diagnostics for the best model (by AIC if available, else R^2).

        - Computes raw, Pearson, and deviance residuals.
        - Performs Hosmer–Lemeshow GOF and overdispersion test.
        - Generates diagnostic plots.
        Returns a dictionary with key statistics and residual arrays.
        """
        if not self.models:
            raise ValueError("Models must be fitted first.")

        comp = self._rank_models()
        best_name = comp.iloc[0]['Model']
        model = self.models[best_name]
        cols = self._predictor_map.get(best_name, [])
        X = sm.add_constant(self.df_cell[cols])
        y_prop = np.asarray(self.df_cell['prop_correct'], dtype=float)
        n = np.asarray(self.df_cell['total_n'], dtype=float) if 'total_n' in self.df_cell.columns else np.ones_like(y_prop)

        # Predicted means (probabilities)
        mu = np.asarray(model.predict(X), dtype=float)
        # guard against 0/1
        eps = 1e-9
        mu = np.clip(mu, eps, 1 - eps)
        y_prop = np.clip(y_prop, eps, 1 - eps)

        # Residuals
        raw_resid = y_prop - mu
        var_prop = np.clip(mu * (1.0 - mu) / np.maximum(n, 1.0), eps, None)
        pearson_resid = raw_resid / np.sqrt(var_prop)

        # Deviance residuals for grouped binomial (using proportions)
        # D_i = 2*n_i [ y*log(y/mu) + (1-y)*log((1-y)/(1-mu)) ]; r_i = sign(y-mu)*sqrt(D_i)
        D_i = 2.0 * n * (y_prop * np.log(y_prop / mu) + (1.0 - y_prop) * np.log((1.0 - y_prop) / (1.0 - mu)))
        D_i = np.where(np.isfinite(D_i) & (D_i >= 0), D_i, 0.0)
        deviance_resid = np.sign(raw_resid) * np.sqrt(D_i)

        # Overdispersion (Pearson X^2 / df)
        X2 = float(np.nansum(pearson_resid ** 2))
        p_params = int(len(getattr(model, 'params', [])))
        df = max(1, len(y_prop) - p_params)
        phi = X2 / df
        overdisp_p = 1.0 - chi2.cdf(X2, df)

        # Hosmer–Lemeshow test
        hl = self._hosmer_lemeshow(mu=mu, n=n, y_count=y_prop * n, n_groups=n_groups)

        # Correlation of residuals with PDI (trend check)
        r_trend, p_trend = stats.pearsonr(pearson_resid, np.asarray(self.df_cell['PDI'], dtype=float))

        print("\n=== Residual analysis (best model: %s) ===" % best_name)
        print(f"Overdispersion: X^2={X2:.2f}, df={df}, phi={phi:.3f}, p={overdisp_p:.4g}")
        print(f"Hosmer–Lemeshow: chi2={hl['chi2']:.2f}, df={hl['df']}, p={hl['p']:.4g}")
        print(f"Residuals vs PDI correlation: r={r_trend:.3f}, p={p_trend:.4g}")

        # Plots
        if show_plots or save_prefix is not None:
            fig, axes = plt.subplots(2, 2, figsize=figsize)

            # 1) Pearson residuals vs fitted
            axes[0, 0].scatter(mu, pearson_resid, alpha=0.7, edgecolor='k', linewidth=0.3)
            axes[0, 0].axhline(0, color='gray', linestyle='--', linewidth=1)
            axes[0, 0].set_xlabel('Fitted probability')
            axes[0, 0].set_ylabel('Pearson residual')
            axes[0, 0].set_title('Residuals vs Fitted')
            axes[0, 0].grid(True, alpha=0.3)

            # 2) Pearson residuals vs PDI
            axes[0, 1].scatter(self.df_cell['PDI'], pearson_resid, alpha=0.7, edgecolor='k', linewidth=0.3)
            axes[0, 1].axhline(0, color='gray', linestyle='--', linewidth=1)
            axes[0, 1].set_xlabel('PDI (ms)')
            axes[0, 1].set_ylabel('Pearson residual')
            axes[0, 1].set_title('Residuals vs PDI')
            axes[0, 1].grid(True, alpha=0.3)

            # 3) QQ plot of Pearson residuals
            osm, osr = stats.probplot(pearson_resid, dist='norm')
            axes[1, 0].scatter(osm[0], osm[1], alpha=0.7, edgecolor='k', linewidth=0.3)
            # Add reference line
            slope, intercept = np.polyfit(osm[0], osm[1], 1)
            xx = np.array([np.min(osm[0]), np.max(osm[0])])
            axes[1, 0].plot(xx, slope * xx + intercept, color='red', linestyle='--', linewidth=1)
            axes[1, 0].set_title('QQ-plot (Pearson residuals)')
            axes[1, 0].set_xlabel('Theoretical quantiles')
            axes[1, 0].set_ylabel('Sample quantiles')
            axes[1, 0].grid(True, alpha=0.3)

            # 4) Histogram of Pearson residuals
            axes[1, 1].hist(pearson_resid, bins=20, color='#2c7fb8', alpha=0.8, edgecolor='white')
            axes[1, 1].set_title('Histogram (Pearson residuals)')
            axes[1, 1].set_xlabel('Residual')
            axes[1, 1].set_ylabel('Count')
            axes[1, 1].grid(True, alpha=0.3)

            plt.suptitle(f'Residual diagnostics: {best_name}', y=1.02)
            plt.tight_layout()

            if save_prefix is not None:
                fig.savefig(f"{save_prefix}_residuals.png", dpi=200, bbox_inches='tight')
                # Also export HL table
                try:
                    hl['table'].to_csv(f"{save_prefix}_hosmer_lemeshow_table.csv", index=False)
                except Exception:
                    pass
            if show_plots:
                plt.show()
            else:
                plt.close(fig)

        return {
            'best_model': best_name,
            'raw_residuals': raw_resid,
            'pearson_residuals': pearson_resid,
            'deviance_residuals': deviance_resid,
            'fitted': mu,
            'overdispersion': {'X2': X2, 'df': df, 'phi': phi, 'p': overdisp_p},
            'hosmer_lemeshow': hl,
            'trend_test': {'r': r_trend, 'p': p_trend}
        }


# Usage
def main():
    """Run the complete analysis."""
    analyzer = PDIAnalyzer()
    analyzer.fit_models()
    analyzer.analyze_trends()
    analyzer.plot_fits()


if __name__ == "__main__":
    main()