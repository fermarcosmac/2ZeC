%% Example of use: 2ZeC
% This example script reads an impulse response from the .\data directory,
% adds gaussian noise with the specified SNRs and crops it using either
% 2ZeC, SWED or MIRACLE methods. To shift between truncation algorithms,
% comment/uncomment lines 46-48. Time-and-frequency domain results are
% plotted, as well as a comparison of a rendered wideband signal
% (.\data\test_signal.wav).

clear, close all
addpath("utils\");

% User-defined parameters first
ir_file = "example_h_bp.wav";
test_signal_file = "test_signal.wav";
data_dir = ".\data\";
algorithm_fcn_name = 'twoZeC';
metric_names = ["MSE" "SDR"];
hyperparam_names = ["p" "spec_tol" "SNR" "f_lims"];

% 2ZeC hyperparameters (spec_tol)
p = Inf;
SNRs = [40 20 10];
f_lims = [0 20e3];

% Retrieve impulse response from IR's directory
ir_path = strcat(data_dir,ir_file);
[h_ref,fs] = audioread(ir_path);

% Retrieve test signal (MLS)
test_signal_path = strcat(data_dir,test_signal_file);
[test_signal, fs2] = audioread(test_signal_path);

% Initialize figure
nplots = length(SNRs);
figure(1), clf
figure(2), clf

for i = 1:nplots
    % Get SNR
    SNR = SNRs(i);

    % Set optimal spectral tolerance
    spec_tol = get_optimal_spec_tol(SNR);

    % Add noise (if any)
    h_ref = add_gaussian_noise(h_ref,SNR);
    
    % Call IR truncation algorithm (2ZeC, SWED or MIRACLE) and return cropped response + limits in original IR
    [h_crop,t_lims,f_lims] = twoZeC(h_ref,fs,p,spec_tol,f_lims);
    %[h_cropped,t_lims] = SWED(h,0,SNR);
    %[h_cropped,t_lims] = MIRACLE(h,0,SNR)
    
    % Plot the original IR (time & freq), the limits and the computed loss
    figure(1)
    subplot(nplots,2,2*(i-1)+1)
    plot(h_ref), hold on, xline(t_lims(1),'r'), xline(t_lims(2),'r');
    grid on
    xlim tight
    ylabel("Amplitude","Interpreter","latex")
    if i == nplots
       xlabel("Samples","Interpreter","latex") 
    end
    if i == 1
        title("Impulse response","Interpreter","latex")
    end
    subplot(nplots,2,2*(i-1)+2)
    NDFT = length(h_ref);
    H_ref = 20*log10(abs(fft(h_ref,NDFT)));
    H_crop = 20*log10(abs(fft(h_crop,NDFT)));
    f = linspace(0,fs,NDFT)*1e-3;
    h1 = plot(f,H_ref,'DisplayName','$H(e^{j\omega})$'); hold on
    h2 = plot(f,H_crop,'DisplayName','$\hat{H}(e^{j\omega})$');
    grid on, xlim([0,fs/2*1e-3]);
    ylabel("dB","Interpreter","latex")
    if i == nplots
        xlabel("f(kHz)","Interpreter","latex")
    end
    if i == 1
        title("Frequency response","Interpreter","latex")
    end

    % Render test signal through both IRs (original and truncated)
    h_pad = [zeros(t_lims(1)-1,1) ; h_crop];
    nfft = max(length(h_ref),length(h_crop)) + length(test_signal) + 1;
    H_ref = fft(h_ref,nfft);
    H_pad = fft(h_pad,nfft);
    X_test = fft(test_signal,nfft);
    y_ref = ifft(H_ref.*X_test,nfft);
    y_crop = ifft(H_pad.*X_test,nfft);

    % Plot rendered test signals comparison  
    figure(2)
    subplot(nplots,1,i)
    plot(y_ref,'DisplayName','$y[n]$'), hold on, grid on
    plot(y_crop,'DisplayName','$\hat{y}[n]$')
    xlabel('Samples','Interpreter','latex')
    ylabel('Amplitude','Interpreter','latex')
    title('Output signals comparison','Interpreter','latex')
    legend('Interpreter','latex')

    % Get error metric
    mse_output = myMSE(y_ref,y_crop);
    fprintf('Output MSE for SNR = %2.2f dB: %4.4f\n',SNR,mse_output)
end

