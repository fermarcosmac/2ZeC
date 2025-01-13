function [h_cropped,t_lims,f_lims] = twoZeC(h,fs,p,spec_tol,varargin)
% Implements 2ZeC algorithm for generic Impulse Response truncation
%
% @params:
%   h: original impulse response (array)
%   fs: impulse response sample rate (double)
%   p: hyperparameter defining p-norm (double)
%   spec_tol: hyperparameter defining spectral tolerance (double)
%   varargin: can hold frequency limits as a vector of the format [f0, f1] (double)
%
% @returns:
%   h_cropped: truncated impulse response (double array)
%   t_lims: vector containing the truncation limits, in samples (integers)
%   f_lims: vector containing the used frequency limits (passed as parameter or estimated)

    % Normalize response
    norm_factor = 1/max(h);
    h = h*norm_factor;

    % Get bandwidth
    NDFT = length(h);
    H = abs(fft(h,NDFT)/sqrt(NDFT));
    if isempty(varargin) || isempty(varargin{1})
        f_lims = estimate_bandwidth(h,fs);
    else
        f_lims = varargin{1};
        assert_valid_f_lims(f_lims);
    end

    % Set initial stride-lengths and crop indices
    state = init_state(h,f_lims,fs,p,spec_tol);

    % Main loop
    timeout = 1; % timeout when 3 iters with static state
    while ( (state(1) >= 1) || (state(2) >= 1) ) && mod(timeout,4)
        prev_state = state;
        state = step(state,h,H,"left");
        state = step(state,h,H,"right");
        if state == prev_state
            timeout = timeout + 1;
        end
    end
     
    % Return
    is = state(3);
    ie = state(4);
    
    t_lims = [is ie];
    h_cropped = h(is:ie)/norm_factor;
end



