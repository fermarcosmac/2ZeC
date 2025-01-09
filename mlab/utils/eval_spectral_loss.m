function [loss,state] = eval_spectral_loss(state,h,H,side)
% Evaluates the spectral loss function to determine whether current
% truncation state is tolerable, and updates state vector accordingly
%
% @params:
%   state: 2ZeC state vector -> see init_state() function
%   h: original impulse response (double array)
%   H: DTF modulus of h (double array)
%   side: string -> "left" or "right"
%
% @returns:
%   loss: spectral loss value (double)
%   state: updated 2ZeC state vector

    % Get current crop limits
    sl = state(1);
    sr = state(2);
    is = state(3);
    ie = state(4);

    % Get interest band samples
    f0 = state(5);
    f1 = state(6);
    fs = state(7);
    NDFT = length(H);
    f0_s = freq2sample(f0,fs,NDFT) + 1;
    f1_s = freq2sample(f1,fs,NDFT) + 1;
    band_samples = f0_s:f1_s;

    % Compute loss
    p = state(8);
    w = h(is:ie);
    W = abs(fft(w,NDFT)/sqrt(NDFT));
    error_vector = W(band_samples) - H(band_samples);

    % Consider 0-norm and p-norm
    if p==0
        loss = sum(error_vector~=0);
    else
        loss = norm(error_vector,p);
    end

    % If error is too high, step back and shorten step
    spec_tol = state(9);
    if loss > spec_tol*norm(H(band_samples),p)
        switch side
            case "left"
                is = is - sl;
                sl = floor(sl/2);
            case "right"
                ie = ie + sr;
                sr = floor(sr/2);
        end
    end

    % Update state
    state(1) = sl;
    state(2) = sr;
    state(3) = is;
    state(4) = ie;

end