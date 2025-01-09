function assert_valid_f_lims(f_lims)
% Asserts validity of user-defined f_lims (frequency limits)
%
% @params:
%   f_lims: vector containing frequency limits, in format [f0, f1]

    assert(isnumeric(f_lims) && isvector(f_lims) && numel(f_lims) == 2, ...
        'f_lims must be a 1x2 numeric vector.');
    assert(f_lims(1) <= f_lims(2), ...
        'Elements of f_lims must be in increasing order.');
    
end