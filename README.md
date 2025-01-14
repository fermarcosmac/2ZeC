# 2ZeC: Two-sided Zeno's Crop

An algorithmic proposal for cropping Impulse Respose (IR) data in an automatic fashion. Said impulse responses may typically be the result of (non)linear system identification procedures. It is common to prune the tails of the obtained response by hand, a process which requires time and human knowledge. 2ZeC proposes a simple iterative solution based on spectral error within the bandwidth of interest.

Both MATLAB and Python implementations of the algorithm are provided. For each language, there is an example script which reads the selected Impulse Response from th /data/ directory and applies 2ZeC to prune it. Results are plotted both in time and frequency domain. Additionally, a wideband test signal (Maximum Length Sequence) is rendered through both the original and the pruned IR, and a visualization is also provided.

Other algorithms (discussed in paper) are also provided. To test them, uncomment corresponding lines of the example scripts (detailed in the header).
