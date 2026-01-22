function StepLocations = FindSteps(input_data, Fs, plotdata)
%FINDSTEPS Detect step events (stair gait) from an IMU time series.
%
% This function returns the *sample indices* where steps occur (not time).
% To convert indices to seconds, divide by Fs.
%
% Intended use:
%   - input_data: typically the foot gyroscope signal (e.g., Gyr_Y in your pipeline).
%   - Fs: sampling frequency in Hz.
%   - plotdata: optional secondary signal for plotting (e.g., Acc_Z) to visually
%     verify detected steps (used only when plotting is enabled below).
%
% Output:
%   StepLocations: vector of step event sample indices.

%% ------------------------------------------------------------------------
% Optional/legacy approach (commented out):
% Peak detection on accelerometer Z (not used in current pipeline).
% ------------------------------------------------------------------------
% low  = 0.75 / (Fs/2);
% band = [0.75 22] / (Fs/2);
% [b, a] = butter(2, band, 'bandpass');
% Filt_data = filtfilt(b, a, input_data);
%
% [all_peaks, locs_peaks] = findpeaks(Filt_data, "MinPeakHeight",3.5, "MinPeakDistance",0);
%
% max_peak = max(all_peaks);
% threshold = max_peak / 1.8;
% first_large_idx = find(all_peaks > threshold, 1, 'first');
% pks_after  = all_peaks(first_large_idx+1:end);
% locs_after = locs_peaks(first_large_idx+1:end);
%
% % Keep only peaks below threshold (smaller peaks)
% valid_idx = pks_after < threshold;
% StepLocations  = locs_after(valid_idx);
% peaklocations  = locs_after / Fs;

%% ------------------------------------------------------------------------
% Current approach: detect "pits" (negative peaks) in the low-pass filtered
% gyroscope signal.
% ------------------------------------------------------------------------

% Low-pass filter the gyroscope signal to reduce noise before peak detection
low = 20 / (Fs/2);
[b, a] = butter(2, low, 'low');
Filt_data = filtfilt(b, a, input_data);

% Invert to detect negative peaks as positive peaks ("pits")
pits = -Filt_data;

% Peak detection parameters:
% - MinPeakHeight / Prominence: reject small deflections
% - MinPeakDistance: avoid double-counting within one step
[all_pits, locs_peaks] = findpeaks( ...
    pits, ...
    "MinPeakHeight",      200, ...
    "MinPeakProminence",  010, ...
    "MinPeakDistance",    35);

% Sort step candidates and remove edge artifacts (start/end buffer)
sorted_locs  = sort(locs_peaks);
start_thresh = 120;
end_thresh   = length(input_data) - 120;

sorted_locs  = sorted_locs(sorted_locs > start_thresh & sorted_locs < end_thresh);

% Output step indices (samples)
StepLocations = sorted_locs;

%% ------------------------------------------------------------------------
% Optional plots to visually inspect detected steps.
% Set q to:
%   1 = plots for accelerometer-based method (legacy / commented-out)
%   2 = plots for gyro-based method (current)
%   3 = no plots
% ------------------------------------------------------------------------
q = 3;

if q == 1
    % Plots for the accelerometer Z approach (legacy block)
    t = (1:length(Filt_data)) / Fs;

    figure()
    subplot(3,1,1)
    plot(t, input_data)
    xline(peaklocations, 'r')
    title('Raw signal with detected (small) peaks')

    subplot(3,1,2)
    plot(t, Filt_data)
    xline(locs_peaks/60, 'r')
    title('Filtered signal with detected (small) peaks')

    subplot(3,1,3)
    plot(t, input_data)
    xline(locs_peaks/60, 'g')
    title('Raw signal with detected (small) peaks (alt)')

elseif q == 2
    % Plots for the gyroscope approach (current block)
    t = (1:length(input_data));

    figure()
    subplot(2,1,1)
    plot(t, Filt_data)
    xline(sorted_locs/60, 'r')
    title('Detected steps over filtered gyroscope signal')

    subplot(2,1,2)
    plot(t, plotdata)
    xline(sorted_locs/60, 'r')
    title('Detected steps over secondary signal (e.g., Acc)')

elseif q == 3
    % No plots
end
end
