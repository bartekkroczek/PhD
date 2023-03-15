using CSV
using GLM
using Chain
using Statistics
using VegaLite
using DataFrames
using DataFramesMeta


corr_by_CSI |>
@vlplot(
    :line,
    x = :CSI,
    y = {:mean_corr, scale = {zero = false}},
    width = 400,
    height = 400
)

part_ids = []
fvals = []
betas = []

for f in filter(x -> endswith(x, ".csv"), readdir("Data/Data_A/"))
    beh = DataFrame(CSV.File("Data/Data_A/$f"))
    corr_by_CSI = @chain beh begin
        @subset(:Trial_type .== "experiment")
        @select(:CSI, :Corr)
        @by(:CSI, :mean_corr = mean(:Corr))
    end
    nullmodel = lm(@formula(mean_corr ~ 1), corr_by_CSI)
    model = lm(@formula(mean_corr ~ CSI), corr_by_CSI)
    ft = ftest(nullmodel.model, model.model)
    fval = round(ft.pval[2], digits=5)
    beta = round.(coef(model), digits=5)[2]
    push!(part_ids, f)
    append!(betas, beta)
    append!(fvals, fval)
end

df = DataFrame(PART_ID = part_ids, beta = betas, fval = fvals)


files = filter(x -> endswith(x, ".csv"), readdir("Data/Data_A/"));
df = vcat([DataFrame(CSV.File("Data/Data_A/$f")) for f in files]...);

@chain df begin
    @subset(:Trial_type .== "experiment")
    # @select(:PART_ID, :CSI, :Corr)
    @by(:PART_ID, :CSI, :mean_corr = mean(:Corr))
    # @by(:PART_ID, beta = coef(lm(@formula(mean_corr ~ 1), _)))
end