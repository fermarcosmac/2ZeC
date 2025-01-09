function noisy_signal = add_gaussian_noise(signal,snr)
% Adds gaussian noise to input signal with specified signal-to-noise ratio
%
% @params:
%   signal: original signal (double array)
%   snr: desired signal-to-noise power ratio [dB] (double)
%
% @returns:
%   noisy_signal: signal contaminated with additive noise

    % Noise synthesis
    signal_power = mean(detrend(signal).^2);    % only AC power
    noise_power = signal_power/(10^(snr/ 10));
    noise = sqrt(noise_power)*randn(size(signal));

    % Return
    noisy_signal = signal + noise;

end