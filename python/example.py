""" 
Example of use: 2ZeC

This example script reads an impulse response from the .\data directory,
adds gaussian noise with the specified SNRs and crops it using either
2ZeC, SWED or MIRACLE methods. To shift between truncation algorithms,
comment/uncomment lines 55-57. Time-and-frequency domain results are
plotted, as well as a comparison of a rendered wideband signal
(.\data\test_signal.wav).

"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.io.wavfile import read as wav_read
from scipy.fft import fft, ifft
from utils import twoZeC, add_gaussian_noise, get_optimal_spec_tol, myMSE, SWED, MIRACLE

# User-defined parameters first
ir_file = "example_h_bp.wav"
test_signal_file = "test_signal.wav"
data_dir = "./data/"
algorithm_fcn_name = 'twoZeC'
metric_names = ["MSE", "SDR"]
hyperparam_names = ["p", "spec_tol", "SNR", "f_lims"]

# 2ZeC hyperparameters
p = np.inf
SNRs = [40, 20, 10]
f_lims = [0, 20e3]

# Retrieve impulse response from IR's directory
ir_path = data_dir + ir_file
fs, h_ref = wav_read(ir_path)
h_ref = h_ref.astype(np.float32)

# Retrieve test signal (MLS)
test_signal_path = data_dir + test_signal_file
fs2, test_signal = wav_read(test_signal_path)
test_signal = test_signal.astype(np.float32)

# Initialize figures
nplots = len(SNRs)
fig1, axs1 = plt.subplots(nplots, 2, figsize=(10, 5 * nplots))
fig2, axs2 = plt.subplots(nplots, 1, figsize=(10, 5 * nplots))

for i, SNR in enumerate(SNRs):
    # Set optimal spectral tolerance
    spec_tol = get_optimal_spec_tol(SNR)

    # Add noise
    h_noisy = add_gaussian_noise(h_ref, SNR)

    # Call IR truncation algorithm (2ZeC, SWED or MIRACLE) and return cropped response + limits in original IR
    h_crop, t_lims, f_lims = twoZeC(h_noisy, fs, p, spec_tol, f_lims)
    #h_crop, t_lims = SWED(h_noisy,0,SNR)
    #h_crop, t_lims = MIRACLE(h_noisy,0,SNR)

    # Plot the original IR (time & freq), the limits, and the computed loss
    axs1[i, 0].plot(h_noisy, label="Noisy IR")
    axs1[i, 0].axvline(x=t_lims[0], color='r', linestyle='--')
    axs1[i, 0].axvline(x=t_lims[1], color='r', linestyle='--')
    axs1[i, 0].grid(True)
    axs1[i, 0].set_xlim([0, len(h_noisy)])
    axs1[i, 0].set_ylabel("Amplitude")
    if i == nplots - 1:
        axs1[i, 0].set_xlabel("Samples")
    if i == 0:
        axs1[i, 0].set_title("Impulse Response")

    NDFT = len(h_noisy)
    H_noisy = 20 * np.log10(np.abs(fft(h_noisy, NDFT)))
    H_crop = 20 * np.log10(np.abs(fft(h_crop, NDFT)))
    f = np.linspace(0, fs / 2, NDFT // 2) * 1e-3
    axs1[i, 1].plot(f, H_noisy[:NDFT // 2], label="Original IR")
    axs1[i, 1].plot(f, H_crop[:NDFT // 2], label="Cropped IR")
    axs1[i, 1].grid(True)
    axs1[i, 1].set_xlim([0, fs / 2 * 1e-3])
    axs1[i, 1].set_ylabel("dB")
    if i == nplots - 1:
        axs1[i, 1].set_xlabel("f (kHz)")
    if i == 0:
        axs1[i, 1].set_title("Frequency Response")

    # Render test signal through both IRs (original and truncated)
    h_pad = np.pad(h_crop, (t_lims[0] - 1, 0), mode='constant')
    nfft = max(len(h_ref), len(h_crop)) + len(test_signal) + 1
    H_ref = fft(h_noisy, nfft)
    H_pad = fft(h_pad, nfft)
    X_test = fft(test_signal, nfft)
    y_ref = ifft(H_ref * X_test, nfft).real
    y_crop = ifft(H_pad * X_test, nfft).real

    # Plot rendered test signals comparison
    axs2[i].plot(y_ref, label="Original Output")
    axs2[i].plot(y_crop, label="Cropped Output")
    axs2[i].grid(True)
    axs2[i].set_xlabel("Samples")
    axs2[i].set_ylabel("Amplitude")
    axs2[i].set_title(f"Output Signals Comparison (SNR={SNR} dB)")
    axs2[i].legend()

    # Get error metric
    mse_output = myMSE(y_ref, y_crop)
    print(f"Output MSE for SNR = {SNR:.2f} dB: {mse_output:.4f}")

plt.tight_layout()
plt.show()