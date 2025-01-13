function [h_cropped,t_lims] = MIRACLE(h,~,varargin)
% Implements the Impulse Response truncation used for the construction of
% the MIRACLE dataset:
%
% A. Kujawski, A. J. R. Pelling, and E. Sarradj, “MIRACLE-microphone array
% impulse response dataset for acoustic learning.”
%
% Which keeps the first samples such that 99.9% of the total energy is
% preserved. Here, a slight modification is introduced, keeping the exact
% sample that fulfills the energy condition instead of the next power of 2,
% since otherwise short responses are barely cropped at all.
%
% @params:
%   h: original impulse response (double array)
%   varargin: can hold original impulse's SNR (double)
%
% @returns:
%   h_cropped: truncated impulse response (double array)
%   t_lims: vector containing the truncation limits, in samples (integers)

    % Get crop limits (s0,s1) that preserve 99.0% of total energy
    s0 = 1;

    cum_energy = cumsum(abs(h).^2);
    cum_energy = cum_energy./cum_energy(end);

    valid_indices = find(cum_energy>0.999);
    s11 = valid_indices(1);
    %s1 = min(2^nextpow2(s11),length(h_ref));   % Uncomment for literal implementation of the algorithm
    s1 = s11;

    % Return
    h_cropped = h(s0:s1);
    t_lims = [s0,s1];

end

