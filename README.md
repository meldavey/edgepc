Edge Pulse Caller

This example is the code included in the science paper:
Real-time dynamic single-molecule protein sequencing on an integrated semiconductor device
Science  Oct 13, 2022   https://www.science.org/doi/10.1126/science.abo7651

Here, a simple time series dataset (trace.out) is input into the detector.  The detector establishes the mean and variance from a small initial portion of the data, and then processes the trace as an online time series streaming dataset, meaning that there is no requirement to store the entire trace in RAM.  The algorithm can correctly detect the rise and fall in signal and trace the "in pulse" regions as well as updating the baseline, all in realtime.  It is designed to be embarisingly parallell with no contstraints on neighbor data.  This is how the Quantum-Si pulse caller scales up to 2M and 8M apertures in realtime (30fps on 4 cores at 2.2Ghz).

The system has a bootstrap mode where the first few seconds are used to establish a baseline mean and variance.  After that, only 6 frames are stored in RAM, along with some state information.  This is used to detect pulses as the data streams in.

This version of the code is called the "edge" pulse caller because it is using a simple edge detection scheme to determine statistical significant changes in baseline.  It is also using the maximum change in baseline for both entering and exiting the pulse in order to improve accuracy.

The code plots the raw data along with the detected pulses.   Included are a couple of png files showing all pulses and a zommed in area of the first 5000 frames of the dataset.

