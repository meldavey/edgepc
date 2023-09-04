import numpy as np
import sys
import os
# don't try and display stuff if we are running on headless systems
canDisplay = True
if os.environ.get('DISPLAY', '')  == '':
    print('Running in no display mode')
    import matplotlib
    matplotlib.use('Agg')
    canDisplay = False
import matplotlib.pyplot as plt

from scipy.optimize import curve_fit

def gauss(x,mu,sigma,A):
    return A*np.exp(-(x-mu)**2/2/sigma**2)

def bimodal(x,mu1,sigma1,A1,mu2,sigma2,A2):
    return gauss(x,mu1,sigma1,A1)+gauss(x,mu2,sigma2,A2)

def fit_bimodal(hist):
    num1 = np.max(hist[:50])
    num2 = np.max(hist[50:])
    expected=(25,10,num1,75,5,num2)
    x = np.arange(0,101)
    params,cov=curve_fit(bimodal,x,hist,expected)
    return params

# solves the intersection of two Gaussians, return is an array of one or more intersections
def solve(m1,m2,std1,std2):
  a = 1/(2*std1**2) - 1/(2*std2**2)
  b = m2/(std2**2) - m1/(std1**2)
  c = m1**2 /(2*std1**2) - m2**2 / (2*std2**2) - np.log(std2/std1)
  return np.roots([a,b,c])

# estimates the baseline and noise for the data points considered not in a pulse
# to do this, we generate a histogram of intensity values, where we expect the pulse
# distribution to be one Gaussian, and the background "noise" values to be the other
# we then fit a bimodal model to the data, and find the intersetion of the two Gaussians
# we can then use all the points "not in a pulse" to calculate the true baseline and noise
# these bootstrap values are required for the edge detector
def estimator(data):
    bootstrap_data = data[:500]
    minval = np.min(bootstrap_data)
    maxval = np.max(bootstrap_data)
    hist = np.zeros((101))
    for i in range(len(bootstrap_data)):
        val = int((bootstrap_data[i] - minval) / (maxval - minval) * 100.0)
        hist[val] += 1
    # hist should contain both the background distribution as well as the pulse intensity distribution
    # find the min point that separates the two distributions by fitting two Gaussians
    params = fit_bimodal(hist)
    intersects = solve(params[0], params[3], params[1], params[4])
    intersect = -1
    print('params: %f %f' % (params[0], params[3]))
    for i in intersects:
        if i > params[0] and i < params[3]:
            intersect = i
    if intersect == -1:
        print('solve issues?')
        intersect = 50
    print('intersetion: %s' % intersect)
    cutoff = minval + (intersect/100.0)*(maxval-minval)
    print('cutoff val: %.1f' % cutoff)
    cutoff_vals = [bootstrap_data[i] for i in range(len(bootstrap_data)) if bootstrap_data[i] < cutoff]
    noise = np.std(cutoff_vals)
    baseline = np.mean(cutoff_vals)
    return baseline, noise

# detects pulses in the data using an online 6 frame window
# the algorithm is designed to run in realtime and requires only the most recent 6 frames of data
# along with some state information for each aperture, and the initial background and noise.
# Two sliding 3-frame window averages are calculated, which can be thought of as the average of the
# data on either side of a potential edge. Optimal edges are those edges which maximize the distance
# between the two windows.
# Edges are evaluated to measure distance from the baseline, and also from each other.  A simple
# state machine is used to detect entry into a pulse, optimal pulse entry point, and exit from a
# pulse.
# While not in a pulse, the background estimation is updated to account for slow variations in
# the background signal.  Noise could similarly be updatd in this way with a running variance approach.
def pulses(data, baseline, noise):
    num_points = len(data)
    thresh = noise * 2.0
    p_cur = 5
    in_pulse = False
    pulse_start = 0
    pulse_end = 0
    pulse_intensity = 0
    prev_edge_step = 0
    pulse_list = []

    while p_cur < num_points:
        p = p_cur - 5
        edge1 = np.sum(data[p:p+3]) / 3
        edge2 = np.sum(data[p+3:p+6]) / 3

        baseline_step = edge2 - baseline

        if baseline_step >= thresh:
            if not in_pulse:
                # if this is the first significant step up away from the baseline, we are in a pulse
                in_pulse = True
                prev_edge_step = edge2 - edge1
                pulse_start = p + 3
                pulse_intensity = 0 # data[pulse_start]
                pulse_start_baseline = baseline
                pulse_start_noise = noise
            elif (edge2 - edge1) > prev_edge_step and ((p+3) - pulse_start) < 3:
                # if, within the next 2 frames, we are still in a pulse but the edge step is larger, this is the true start of the frame
                # part of the logic here is that the edges are averaging across 3 frames, so for large steps up, we will start entering
                # into a pulse too early, this logic helps find the largest step up with is the most likely start of a pulse
                prev_edge_step = edge2 - edge1
                pulse_start = p + 3
                pulse_intensity = 0 # data[pulse_start]

        elif in_pulse and baseline_step < thresh*0.5:
            # here, we are in a pulse and observe that the step down has now returned to that baseline, so we can exit the pulse
            in_pulse = False
            prev_edge_step = 0
            pulse_end = p + 3 - 1 # don't include this frame because the pulse already ended
            num_frames = pulse_end - pulse_start + 1
            if num_frames < 10000: # limit to reasonable pulse durations
                pulse_list.append([pulse_start, num_frames, pulse_intensity, pulse_start_baseline, pulse_start_noise])

        # outside of the state transitions, if we are still in a pulse, record the intensity of the frame representing the first frame of the step-up
        # otherwise, we are not in a pulse, and can update the baseline
        if in_pulse:
            pulse_intensity += data[p+3]
        else:
            baseline = 0.999 * baseline + 0.001 * data[p+3]

        p_cur += 1

    return pulse_list

def edge(data):
    # bootstrab up and estimate a few params
    baseline, noise = estimator(data)
    print('baseline: %.1f  noise: %.1f' % (baseline, noise))

    # generate a list of pulses
    # each pulse contains: start (frame), duration (frames), intensity(summed over duration), baseline_at_start, noise_at_start
    pulse_list = pulses(data, baseline, noise)
    return pulse_list

# load up our sample data
with open('trace.out') as f:
    lines = f.readlines()
    data = [int(x) for x in lines]
data = np.array(data)

# perform edge pulse calling on the raw trace sample
pulses = edge(data)
print('found %d pulses' % len(pulses))

# generate pulse trace so we can plot it
x = []
y = []
x.append(0)
y.append(pulses[0][3])
for pulse in pulses:
    x.append(pulse[0]) # start
    pulse_intensity = pulse[2]/pulse[1]
    y.append(pulse[3]) # baseline

    x.append(pulse[0]) # start
    y.append(pulse_intensity)

    x.append(pulse[0] + pulse[1]) # duration
    y.append(pulse_intensity)

    x.append(pulse[0] + pulse[1]) # duration
    y.append(pulse[3]) # baseline

# plot the raw trace and the pulses
fig = plt.figure('trace')
plt.plot(data)
plt.plot(x,y)
plt.show()

