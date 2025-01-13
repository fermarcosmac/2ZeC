function [h_cropped,t_lims] = SWED(h,~,varargin)
% Implements the Sliding-Window Energy Detection (SWED) approach to
% impulsive signal time-limit detection proposed in:
% 
% J. L. Blanco-Murillo and V. Yagüe-Jiménez, “A method for in-
% formed selection of memory-length and nonlinearity-order parameters
% in volterra–wiener systems from exponential sweep excitations,” Multi-
% dimensional Systems and Signal Processing, vol. 29, pp. 1861–1893, 10
% 2018
%
% @params:
%   h: original impulse response (double array)
%   varargin: can hold original impulse's SNR (double)
%
% @returns:
%   h_cropped: truncated impulse response (double array)
%   t_lims: vector containing the truncation limits, in samples (integers)


    % Get SNR if provided
    if ~isempty(varargin)
        SNR = varargin{1};
    else
        SNR = 40; % Default value
    end

    % Set energy detection threshold
    threshold = (1/(2*length(h)))/(1+SNR);

    % Rolling-window parameters
    win_len = 8;                        % Hyperparameter: arbitrary
    in_packet = false;                  % Flag indicating if we are in a packet or not
    packet_lims = [];                   % This will be a matrix with the limits of all found packets
    stride = round(win_len/4);          % Inter-window stride. Means there is window overlapping

    % Power detection loop
    for i = 1:stride:length(h)
        % Extract window
        end_idx = min(i+win_len,length(h));     % Make sure that windows do not exceed signal bounds
        win = h(i:end_idx);                     % Current window
        
        % Compute window power
        win_pow = sum(win.^2)/length(win);

        % Make a decision
        if (win_pow < threshold) && ~in_packet
            % Surfing a dead-zone
        elseif (win_pow < threshold) && in_packet
            % Stepping out of packet
            packet_lims(size(packet_lims,1),2) = i+win_len;
            in_packet = false;
        elseif (win_pow > threshold) && in_packet
            % Surfing a packet
        else
            % Entering a packet
            packet_lims = [packet_lims ; [i+stride , 0]];
            in_packet = true;
        end
    end

    % Keep largest energy packet
    [~,idx] = max(packet_lims(:,1)-packet_lims(:,2));
    s0 = max(packet_lims(idx,1),1);
    s1 = min(packet_lims(idx,2),length(h));

    % Return
    t_lims = [s0, s1];
    h_cropped = h(s0:s1);
    

end

