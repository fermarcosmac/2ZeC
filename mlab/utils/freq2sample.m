function sample = freq2sample(f,fs,NDFT)
% Finds DFT sample corresponding to a given frequency
%
% @params:
%   f: frequency to transform [Hz] (double)
%   fs: sample rate
%   NDTF: number of DFT points used (int)
%
% @returns:
%   sample: closest DTF sample corresponding to frequency f (int)

    sample = round((NDFT/fs)*f);
    
end