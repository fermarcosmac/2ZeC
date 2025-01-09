function spec_tol = get_optimal_spec_tol(SNR)
% Calculates 2ZeC's optimal spectral tolerance hyperparameter by means of
% the empirical equation found with synthetic dataset
%
% @params:
%   SNR: impulse response signal-to-noise ratio [dB] (double)
% @returns
%   spec_tol: optimal spectral tolerance

    spec_tol = 10^((-4.1123-0.521*SNR+0.0035*SNR^2)/10);

end

