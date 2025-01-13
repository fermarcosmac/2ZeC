import numpy as np
from scipy.fft import fft
from scipy.signal import welch


def twoZeC(h, fs, p, spec_tol, f_lims=None):
    """
    Implements 2ZeC algorithm for generic Impulse Response truncation.

    Parameters:
        h (np.ndarray): Original impulse response.
        fs (float): Impulse response sample rate.
        p (float): Hyperparameter defining p-norm.
        spec_tol (float): Hyperparameter defining spectral tolerance.
        f_lims (list, optional): Frequency limits as [f0, f1].

    Returns:
        h_cropped (np.ndarray): Truncated impulse response.
        t_lims (list): Vector containing truncation limits, in samples.
        f_lims (list): Vector containing the used frequency limits.
    """
    # Normalize response
    norm_factor = 1 / np.max(np.abs(h))
    h = h * norm_factor

    # Get bandwidth
    NDFT = len(h)
    H = np.abs(fft(h, NDFT)) / np.sqrt(NDFT)
    if f_lims is None or len(f_lims) == 0:
        f_lims = estimate_bandwidth(h, fs)
    else:
        assert_valid_f_lims(f_lims)

    # Set initial stride lengths and crop indices
    state = init_state(h, f_lims, fs, p, spec_tol)

    # Main loop
    timeout = 1  # Timeout when 3 iterations with static state
    while ((state[0] >= 1) or (state[1] >= 1)) and timeout % 4:
        prev_state = state.copy()
        state = step(state, h, H, "left")
        state = step(state, h, H, "right")
        if np.array_equal(state, prev_state):
            timeout += 1

    # Return
    is_idx, ie_idx = state[2], state[3]
    t_lims = [is_idx, ie_idx]
    h_cropped = h[is_idx:ie_idx] / norm_factor

    return h_cropped, t_lims, f_lims


def estimate_bandwidth(h, fs):
    """
    Estimates the bandwidth of the original impulse response h.

    Parameters:
        h (np.ndarray): Original impulse response.
        fs (int): Original impulse sample rate.

    Returns:
        f_lims (list): Estimated frequency limits as [f0, f1].
    """
    # Default to Nyquist frequency as bandwidth limits
    f_lims = [0, round(fs / 2)]

    # Uncomment the code below for a PSD-based estimate
    # f, phh = welch(h, nperseg=8, fs=fs)  # PSD estimate
    # peak = np.max(phh)                   # Max power
    # passband_freqs = f[phh > (peak / 2)] # -3dB band
    # f_lims = [passband_freqs[0], passband_freqs[-1]]  # Band limits

    return f_lims


def assert_valid_f_lims(f_lims):
    """
    Asserts the validity of user-defined frequency limits (f_lims).

    Parameters:
        f_lims (list or np.ndarray): Vector containing frequency limits in the format [f0, f1].

    Raises:
        ValueError: If f_lims is not a valid 1x2 numeric vector or if the values are not in increasing order.
    """
    if not (isinstance(f_lims, (list, np.ndarray)) and len(f_lims) == 2 and all(isinstance(f, (int, float)) for f in f_lims)):
        raise ValueError("f_lims must be a 1x2 numeric vector.")
    
    if not (f_lims[0] <= f_lims[1]):
        raise ValueError("Elements of f_lims must be in increasing order.")


def init_state(h, f_lims, fs, p, spec_tol):
    """
    Initializes the 2ZeC state vector.

    Parameters:
        h (np.ndarray): Original impulse response.
        f_lims (list): Frequency limits as [f0, f1].
        fs (float): Impulse response sample rate.
        p (float): Hyperparameter defining p-norm.
        spec_tol (float): Hyperparameter defining spectral tolerance.

    Returns:
        state (list): Initial 2ZeC state vector.

    The state vector is defined as:
        1 - sl: Left stride
        2 - sr: Right stride
        3 - is: Start index
        4 - ie: End index
        5 - f0: Start frequency
        6 - f1: End frequency
        7 - fs: Sampling frequency
        8 - p: p-norm for error calculations
        9 - spec_tol: Spectral tolerance
    """
    N = len(h)
    sl = N // 4  # Left stride
    sr = sl      # Right stride
    is_idx = 0   # Start index (Python uses 0-based indexing)
    ie_idx = N   # End index
    f0 = f_lims[0]
    f1 = f_lims[1]

    # Create state vector
    state = [sl, sr, is_idx, ie_idx, f0, f1, fs, p, spec_tol]
    return state


def step(state, h, H, side):
    """
    Step function for the 2ZeC IR truncation algorithm, defining actions at each iteration.

    Parameters:
        state (list): 2ZeC state vector (see init_state() function for structure).
        h (np.ndarray): Original impulse response.
        H (np.ndarray): DTF modulus of h.
        side (str): Direction of the step, either "left" or "right".

    Returns:
        state (list): Updated 2ZeC state vector.
    """
    sl, sr, is_idx, ie_idx = state[:4]

    if side == "left":
        crop_condition = (sl > 1) or ((sl <= 1) and (sr <= 1))
        is_new = min(is_idx + sl, ie_idx - 1)
        ie_new = ie_idx
    elif side == "right":
        crop_condition = (sr > 1) or ((sr <= 1) and (sl <= 1))
        is_new = is_idx
        ie_new = max(ie_idx - sr, is_idx + 1)
    else:
        raise ValueError("[2ZeC] Unknown side: must be 'left' or 'right'.")

    if crop_condition:
        is_idx = is_new
        ie_idx = ie_new

        # Update state with new indices
        state[2] = is_idx
        state[3] = ie_idx

        # Evaluate spectral loss and update the state
        _, state = eval_spectral_loss(state, h, H, side)

    return state


def eval_spectral_loss(state, h, H, side):
    """
    Evaluates the spectral loss function to determine whether the current truncation state is tolerable
    and updates the state vector accordingly.

    Parameters:
        state (list): 2ZeC state vector (see init_state() function for structure).
        h (np.ndarray): Original impulse response.
        H (np.ndarray): DTF modulus of h.
        side (str): Direction of the step, either "left" or "right".

    Returns:
        loss (float): Spectral loss value.
        state (list): Updated 2ZeC state vector.
    """
    # Get current crop limits
    sl, sr, is_idx, ie_idx = state[:4]
    f0, f1, fs = state[4:7]
    p = state[7]
    spec_tol = state[8]
    NDFT = len(H)

    # Get interest band samples
    f0_s = freq2sample(f0, fs, NDFT)
    f1_s = freq2sample(f1, fs, NDFT)
    band_samples = slice(f0_s, f1_s + 1)  # Include f1_s

    # Compute the loss
    w = h[is_idx:ie_idx]
    W = np.abs(np.fft.fft(w, NDFT) / np.sqrt(NDFT))
    error_vector = W[band_samples] - H[band_samples]

    # Consider 0-norm and p-norm
    if p == 0:
        loss = np.sum(error_vector != 0)
    else:
        loss = np.linalg.norm(error_vector, p)

    # If error is too high, step back and shorten step
    if loss > spec_tol * np.linalg.norm(H[band_samples], p):
        if side == "left":
            is_idx -= sl
            sl = max(1, sl // 2)  # Ensure stride is at least 1
        elif side == "right":
            ie_idx += sr
            sr = max(1, sr // 2)  # Ensure stride is at least 1
        else:
            raise ValueError("[2ZeC] Unknown side: must be 'left' or 'right'.")

    # Update state
    state[0] = sl
    state[1] = sr
    state[2] = is_idx
    state[3] = ie_idx

    return loss, state


def freq2sample(f, fs, NDFT):
    """
    Finds the DFT sample corresponding to a given frequency.

    Parameters:
        f (float): Frequency to transform [Hz].
        fs (float): Sample rate [Hz].
        NDFT (int): Number of DFT points used.

    Returns:
        int: Closest DFT sample corresponding to frequency f.
    """

    return int(round((NDFT/fs)*f))



def get_optimal_spec_tol(SNR):
    """
    Calculates 2ZeC's optimal spectral tolerance hyperparameter using 
    an empirical equation derived from a synthetic dataset.

    Parameters:
        SNR (float): Impulse response signal-to-noise ratio [dB].

    Returns:
        float: Optimal spectral tolerance.
    """
    spec_tol = 10 ** ((-4.1123 - 0.521 * SNR + 0.0035 * SNR**2) / 10)
    return spec_tol


def add_gaussian_noise(signal, snr):
    """
    Adds Gaussian noise to the input signal with the specified signal-to-noise ratio.

    Parameters:
        signal (numpy.ndarray): Original signal.
        snr (float): Desired signal-to-noise power ratio [dB].

    Returns:
        numpy.ndarray: Signal contaminated with additive Gaussian noise.
    """
    # Noise synthesis
    signal_power = np.mean((signal - np.mean(signal))**2)  # Only AC power
    noise_power = signal_power / (10**(snr / 10))
    noise = np.sqrt(noise_power) * np.random.randn(*signal.shape)

    # Return noisy signal
    noisy_signal = signal + noise
    return noisy_signal


def myMSE(x_ref, x_deg):
    """
    Computes the Mean Squared Error (MSE) between two vectors.

    Parameters:
        x_ref (numpy.ndarray): Reference vector.
        x_deg (numpy.ndarray): Degraded vector.

    Returns:
        float: Mean Squared Error between the two vectors.
    """
    # Check if lengths match
    if len(x_ref) != len(x_deg):
        raise ValueError("Compared vectors must be of equal lengths.")
    
    # Compute MSE
    N = len(x_ref)
    mse = np.mean((x_ref - x_deg)**2)
    
    return mse



def SWED(h, *varargin):
    """
    Implements the Sliding-Window Energy Detection (SWED) approach to
    impulsive signal time-limit detection.

    Parameters:
        h (numpy.ndarray): Original impulse response.
        varargin (float, optional): SNR value (default is 40).

    Returns:
        numpy.ndarray: Truncated impulse response.
        tuple: Truncation limits in samples (integers).
    """
    # Get SNR if provided, else set default value
    SNR = varargin[0] if varargin else 40  # Default value is 40 if not provided

    # Set energy detection threshold
    threshold = (1 / (2 * len(h))) / (1 + SNR)

    # Rolling-window parameters
    win_len = 8  # Hyperparameter: arbitrary
    in_packet = False  # Flag indicating if we are in a packet or not
    packet_lims = []  # This will be a list with the limits of all found packets
    stride = round(win_len / 4)  # Inter-window stride. Means there is window overlapping

    # Power detection loop
    for i in range(0, len(h), stride):
        # Extract window
        end_idx = min(i + win_len, len(h))  # Make sure that windows do not exceed signal bounds
        win = h[i:end_idx]  # Current window
        
        # Compute window power
        win_pow = np.sum(win**2) / len(win)

        # Make a decision based on window power
        if win_pow < threshold and not in_packet:
            # Surfing a dead-zone
            pass
        elif win_pow < threshold and in_packet:
            # Stepping out of packet
            packet_lims[-1][1] = i + win_len
            in_packet = False
        elif win_pow > threshold and in_packet:
            # Surfing a packet
            pass
        else:
            # Entering a packet
            packet_lims.append([i + stride, 0])
            in_packet = True

    # Keep largest energy packet
    idx = np.argmax(np.array(packet_lims)[:, 0] - np.array(packet_lims)[:, 1])
    s0 = max(packet_lims[idx][0], 0)
    s1 = min(packet_lims[idx][1], len(h))

    # Return
    t_lims = (s0, s1)
    h_cropped = h[s0:s1]
    
    return h_cropped, t_lims


def MIRACLE(h, *varargin):
    """
    Implements the Impulse Response truncation used for the construction of
    the MIRACLE dataset. The function keeps the first samples such that 99.9% 
    of the total energy is preserved.

    Parameters:
        h (numpy.ndarray): Original impulse response.
        varargin (float, optional): SNR value (not used in this function).

    Returns:
        numpy.ndarray: Truncated impulse response.
        tuple: Truncation limits in samples (integers).
    """
    # Get crop limits (s0, s1) that preserve 99.9% of total energy
    s0 = 0  # Python indexing starts from 0

    cum_energy = np.cumsum(np.abs(h)**2)
    cum_energy /= cum_energy[-1]

    # Find the index where 99.9% of the energy is preserved
    valid_indices = np.where(cum_energy > 0.999)[0]
    s1 = valid_indices[0]

    # Return
    h_cropped = h[s0:s1+1]
    t_lims = (s0, s1)

    return h_cropped, t_lims