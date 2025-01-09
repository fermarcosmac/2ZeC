function state = init_state(h,f_lims,fs,p,spec_tol)
% Initializes 2ZeC state vector (specified below)
%
% @params:
%   h: original impulse response (double array)
%   f_lims: frequency limits, specified as a vector of the format [f0, f1] (double)
%   fs: impulse response sample rate (double)
%   p: hyperparameter defining p-norm (double)
%   spec_tol: hyperparameter defining spectral tolerance (double)
%
% @returns:
%   state: initial 2ZeC state vector
%
%    state vector is defined as:
%      1 - sl: left stride
%      2 - sr: right stride
%      3 - is: start index
%      4 - ie: end index
%      5 - f0: start frequency
%      6 - f1: end frequency
%      7 - fs: sampling frequency
%      8 - p : p-norm for error calculations
%      9 - spec_tol: spectral tolerance

    N = length(h);
    sl = floor(N/4);
    sr = sl;
    is = 1;
    ie = N;
    f0 = f_lims(1);
    f1 = f_lims(2);

    state = [sl,sr,is,ie,f0,f1,fs,p,spec_tol];

end

