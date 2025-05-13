n = 10;  % Filter order
f = [0 0.4 0.5 1];  % Normalized frequencies
a = [1 1 0 0];  % Desired response
b = firls(n, f, a);  % FIR filter coefficients

% Save coefficients for R
save('fir_coeff.txt', 'b', '-ascii');

