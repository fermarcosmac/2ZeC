%% Example of use: 2ZeC
close all

% User-defined parameters first
ir_file = "h_bp_0006.wav";
data_dir = ".\impulses\";
algorithm_fcn_name = 'twoZeC';
metric_names = ["MSE" "SDR"];
hyperparam_names = ["p" "spec_tol" "SNR" "f_lims"];

% 2ZeC hyperparameters
p = 2;
spec_tol = 1e-1;
SNRs = [40 20 10];
f_lims = [0 20e3];
hyperparams = {p, spec_tol, SNR, f_lims};

% Add path to functions
addpath("utils\");
addpath("classes\");

% Retrieve impulse response from IR's directory
data_path = strcat(data_dir,ir_file);
[h_ref,fs] = audioread(data_path);

% Initialize figure
nplots = length(SNRs);
figure(1), clf
figure(2), clf

for i = 1:nplots
    % Get SNR
    SNR = SNRs(i);

    % Add noise (if any)
    h_ref = add_gaussian_noise(h_ref,SNR);
    
    % Call 2ZeC and return cropped response + limits in original IR
    [h_crop,t_lims,f_lims] = twoZeC(h_ref,fs,p,spec_tol,f_lims);
    
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
    %xline(f_lims(1),'r'), xline(f_lims(2),'r');
    grid on, xlim([0,fs/2*1e-3]);
    ylabel("dB","Interpreter","latex")
    if i == nplots
        xlabel("f(kHz)","Interpreter","latex")
    end
    if i == 1
        title("Frequency response","Interpreter","latex")
    end
    %legend([h1,h2],'Interpreter','latex')
    
    % Plot test signals comparison 
    TZL = TwoZeCLogger(data_dir,algorithm_fcn_name,metric_names,hyperparam_names,hyperparams);
    TZL.set_current_fs(fs);
    [y_ref, y_crop] = TZL.assessCrop(ir_file,h_crop,h_ref,t_lims);
    
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
    fprintf('Output MSE for SNR = %2.2f dB: %4.4f',SNR,mse_output)
end

