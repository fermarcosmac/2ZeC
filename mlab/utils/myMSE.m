function mse = myMSE(x_ref,x_deg)
% Computes Mean Squared Error (MSE) between two vectors
%
% @params:
%   x_ref: reference vector (double array)
%   x_deg: degraded vector (double array)

    % Check lengths
    assert(length(x_ref)==length(x_deg),"Compared vectors must be of equal lengths.")
    N = length(x_ref);

    % Compute MSE
    mse = (1/N) * sum( (x_ref-x_deg).^2 );
    
end

