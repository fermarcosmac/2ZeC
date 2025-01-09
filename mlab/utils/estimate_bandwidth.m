function f_lims = estimate_bandwidth(h,fs)
% Estimates the bandwidth of original impulse response h. Given the
% recommended use of the infinity norm with 2ZeC, bandwidth estimation has
% proven not to be that relevant. Uncomment code for a simple PSD-based
% estimate
%
% @params:
%   h: original impulse response (double array)
%   fs: original impulse sample rate (int)
%
% @returns:
%   f_lims: estimated frequency limits, as a vector of format [f0, f1] (double)
    
    f_lims = [0, round(fs/2)];

    %[phh,f] = pwelch(h,8,[],[],fs);                     % PSD estimate
    %[peak,~] = max(phh);                                % Max power
    %
    %passband_freqs = f(phh > (peak/2));                 % -3dB band
    %f_lims = [passband_freqs(1) passband_freqs(end)];   % Band limits

end