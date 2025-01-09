function state = step(state,h,H,side)
% Step function is a part of the 2ZeC IR truncation algorithm, which
% defines the actions taken at each iteration.
%
% @params:
%   state: 2ZeC state vector -> see init_state() function
%   h: original impulse response (double array)
%   H: DTF modulus of h (double array)
%   side: string -> "left" or "right"
%
% @returns:
%   state: updated 2ZeC state vector

    sl = state(1);
    sr = state(2);
    is = state(3);
    ie = state(4);

    switch side
        case "left"
            crop_condition = (sl > 1) || ((sl <= 1) && (sr <= 1));
            is_new = min([is + sl , ie - 1]);
            ie_new = ie;
        case "right"
            crop_condition = (sr > 1) || ((sr <= 1) && (sl <= 1));
            is_new = is;
            ie_new = max([ie - sr , is + 1]);
        otherwise
            error("[2ZeC] Unknown side: must be \'left\' or \'right\'");
    end

    if crop_condition
        is = is_new;
        ie = ie_new;

        state(3) = is;
        state(4) = ie;

        [~,state] = eval_spectral_loss(state,h,H,side);
    end



end

