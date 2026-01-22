%% Movella DOT preprocessing + feature extraction (Masters thesis)
% Reads Movella DOT CSV files per participant and extracts:
% - Step events (per foot)
% - Pelvis mediolateral sway (Euler/roll)
% - Gravity-removed acceleration (quat-based)
% - Nonlinear metrics: Sample Entropy, Harmonic Ratio, RQA (CRP toolbox + TISEAN AMI)
% - Linear metrics: peaks, RMS, step frequency variability
%
% Notes:
% - This script assumes a specific folder layout: base_path/P{n}/<Sensor>_<Trial>.csv
% - External dependencies: lagcorrect(), FindSteps(), apply_step_corrections(), combine_trials(),
%   removeGravityWithQuat(), sampen(), CRP toolbox (Marwan), TISEAN mutual.exe, etc.
% - The computational logic is intentionally left unchanged; only comments + spacing were improved.

clc; clear; close all;

%% ========================= Participant-level loop =========================
for i_participant = 2:41  % Loop over all participants (P2..P41)
    base_path  = 'C:\Users\boris\OneDrive\Documenten\M Human Movement Sciences\MRP\data\Movella Data';
    folderPath = fullfile(base_path, sprintf('P%d', i_participant));
    Fs         = 60;  % Sampling frequency (Hz)

    %% ============================= Load data ==============================
    participant_id = sprintf('P%d', i_participant);

    % Trial names as used in the CSV filename suffix
    all_trials = {'Calibratie_Start','Ascent1','Ascent2','Descent1','Descent2','Calibratie_end'};

    % Sensor prefixes as used in the CSV filenames
    sensors = {'Pelvis','LFoot','RFoot','LWrist','RWrist'};

    % Containers for raw tables and time vectors per sensor/trial
    rawdata = struct();
    times   = struct();

    % Read each available CSV and create a time vector based on Fs
    for i_sensor = 1:length(sensors)
        sensor = sensors{i_sensor};
        rawdata.(sensor) = struct();
        times.(sensor)   = struct();

        for i_trial = 1:length(all_trials)
            trial    = all_trials{i_trial};
            filename = sprintf('%s_%s.csv', sensor, trial);
            filePath = fullfile(folderPath, filename);

            if isfile(filePath)
                tbl = readtable(filePath, 'VariableNamingRule', 'preserve');
                rawdata.(sensor).(trial) = tbl;

                n_samples = size(tbl, 1);
                times.(sensor).(trial) = (0:n_samples-1) / Fs;
            else
                warning('File not found: %s', filePath);
            end
        end
    end

    %% ============ Cross-correlation sync to estimate inter-sensor lag ============
    % Determine alignment based on the calibration trials only (start + end).
    % Note: optimal_lag is overwritten per datatype loop; the downstream code uses
    % the last computed optimal_lag in this block (kept as-is to preserve behavior).
    datatypes = {'Acc_X','Acc_Y','Acc_Z'};

    for i = [1,6]  % Use Calibratie_Start and Calibratie_end to estimate lag
        trialname = all_trials{i};

        for j = 1:length(datatypes)
            datatypename = datatypes{j};

            signal = {
                rawdata.Pelvis.(trialname).(datatypename), ...
                rawdata.LFoot.(trialname).(datatypename), ...
                rawdata.RFoot.(trialname).(datatypename)
            };

            time_vector = {
                times.Pelvis.(trialname), ...
                times.LFoot.(trialname), ...
                times.RFoot.(trialname)
            };

            % Align signals and obtain optimal lag shift (samples)
            [aligned_signals, optimal_lag] = lagcorrect(signal, time_vector);
        end
    end

    %% ============================= Apply time-lag correction =============================
    % For each trial and datatype:
    % - shift sensor time vectors using optimal_lag
    % - interpolate foot signals onto pelvis time grid
    % - trim/align by construction to the pelvis reference vector

    if i_participant == 1
        datatypes = {'Acc_X','Acc_Y','Acc_Z','Gyr_X','Gyr_Y','Gyr_Z','Euler_X','Euler_Y','Euler_Z'};
    else
        datatypes = {'Quat_W','Quat_X','Quat_Y','Quat_Z','Acc_X','Acc_Y','Acc_Z', ...
                     'Gyr_X','Gyr_Y','Gyr_Z','Mag_X','Mag_Y','Mag_Z'};
    end

    for i = 1:length(all_trials)
        trialname = all_trials{i};

        for j = 1:length(datatypes)
            datatypename = datatypes{j};

            signal = {
                rawdata.Pelvis.(trialname).(datatypename), ...
                rawdata.LFoot.(trialname).(datatypename), ...
                rawdata.RFoot.(trialname).(datatypename)
            };

            time_vector = {
                times.Pelvis.(trialname), ...
                times.LFoot.(trialname), ...
                times.RFoot.(trialname)
            };

            reference_signal = signal{1};
            reference_time   = time_vector{1};

            % Only align LFoot and RFoot onto Pelvis (u=2..3)
            for u = 2:3
                currentsensor = sensors{u};

                % Convert lag (samples) to time shift (seconds)
                time_shift = optimal_lag(u) * (time_vector{u}(2) - time_vector{u}(1));

                % Shift time grid and interpolate onto pelvis time vector
                adj_time = time_vector{u} - time_shift;

                data.(participant_id).(currentsensor).(trialname).(datatypename) = ...
                    interp1(adj_time, signal{u}, reference_time, 'linear', 0); % out-of-bounds -> 0

                % Store pelvis as raw table (unaltered)
                data.(participant_id).Pelvis.(trialname) = rawdata.Pelvis.(trialname);
            end
        end

        % Store common time vector used for aligned data
        t.(participant_id).(trialname) = reference_time;
    end

    %% ============================= Step recognition =============================
    % Detect step events per trial for each foot using gyroscope + vertical accel.
    feet  = {'LFoot','RFoot'};
    trials = {'Ascent1','Ascent2','Descent1','Descent2'};

    for ii = 1:2
        foot = feet{ii};

        for jj = 1:4
            trial_step = trials{jj};

            StepLocations = FindSteps( ...
                data.(participant_id).(foot).(trial_step).Gyr_Y, ...
                Fs, ...
                data.(participant_id).(foot).(trial_step).Acc_Z);

            rawresults.(participant_id).steps.(foot).(trial_step) = StepLocations;
        end
    end

    %% =================== Optional: diagnostic plots for step detection ===================
    plots = 0;  % 0 = no plots, 1 = plot step overlays for visual inspection

    if plots == 1
        low = 20 / (Fs/2);
        [b, a] = butter(2, low, 'low');

        for jj = 1:length(trials)
            trial_step = trials{jj};

            gyr_L = filtfilt(b, a, data.(participant_id).LFoot.(trial_step).Gyr_Y);
            gyr_R = filtfilt(b, a, data.(participant_id).RFoot.(trial_step).Gyr_Y);
            acc_L = filtfilt(b, a, data.(participant_id).LFoot.(trial_step).Acc_Z);
            acc_R = filtfilt(b, a, data.(participant_id).RFoot.(trial_step).Acc_Z);

            % Sample index used as x-axis (kept as-is)
            times = (1:length(gyr_L));

            stepLocs_L = rawresults.(participant_id).steps.LFoot.(trial_step);
            stepLocs_R = rawresults.(participant_id).steps.RFoot.(trial_step);

            figure('Name', sprintf('%s - %s', participant_id, trial_step), ...
                   'NumberTitle','off', 'Units','normalized', 'OuterPosition',[0 0 1 1]);

            subplot(2,2,1)
            plot(times, gyr_L); hold on; xline(stepLocs_L, 'r');
            title(sprintf('%s - %s - LFoot Gyr_Y', participant_id, trial_step));
            ylabel('Gyr_Y'); xlabel('Time (samples)');

            subplot(2,2,2)
            plot(times, acc_L); hold on; xline(stepLocs_L, 'r');
            title(sprintf('%s - %s - LFoot Acc_Z', participant_id, trial_step));
            ylabel('Acc_Z'); xlabel('Time (samples)');

            subplot(2,2,3)
            plot(times, gyr_R); hold on; xline(stepLocs_R, 'r');
            title(sprintf('%s - %s - RFoot Gyr_Y', participant_id, trial_step));
            ylabel('Gyr_Y'); xlabel('Time (samples)');

            subplot(2,2,4)
            plot(times, acc_R); hold on; xline(stepLocs_R, 'r');
            title(sprintf('%s - %s - RFoot Acc_Z', participant_id, trial_step));
            ylabel('Acc_Z'); xlabel('Time (samples)');
        end
    end

    %% ============== Apply manual corrections to step events (if any) ==============
    new_results = apply_step_corrections2(rawresults, participant_id);
    results.(participant_id) = new_results.(participant_id);

    %% === Optional: plot foot step events over pelvis acceleration (sanity check) ===
    plotjes = 0;  % 0 = off, 1 = on

    if plotjes == 1
        figure('Name', sprintf('%s - %s', participant_id, trial), ...
               'NumberTitle','off', 'Units','normalized', 'OuterPosition',[0 0 1 1]);

        low = 5 / (Fs/2);
        [b, a] = butter(2, low, 'low');

        for z = 1:4
            trial = trials{z};

            subplot(2,2,z)
            plot(t.(participant_id).(trial), filtfilt(b, a, data.(participant_id).Pelvis.(trial).Acc_Y)); hold on;
            hL = xline(results.(participant_id).steps.LFoot.(trial)/Fs, 'g');
            hR = xline(results.(participant_id).steps.RFoot.(trial)/Fs, 'r');

            title(sprintf('%s - %s - Steps on Pelvis ACC-Y', participant_id, trial));
            ylabel('Acc (m/s^2'); xlabel('Time (s)');

            if z == 4
                legend([hL(1), hR(1)], {'LFoot','RFoot'}, 'Location','southeast');
            end
        end
    end
end

%% ============================= Euler angles (pelvis) =============================
% Convert pelvis quaternions to Euler angles (XYZ convention), then compute
% mediolateral sway as the filtered roll component.

for i = 2:41  % P1 had Euler but no quaternions (comment preserved)
    participant_id = sprintf('P%d', i);

    for j = 1:4
        trial = trials{j};

        if i == 1
            % Placeholder for P1 logic (kept as-is)
        else
            n_samples = size(data.(participant_id).Pelvis.(trial).Acc_X);

            % Convert quaternions to Euler angles (roll, pitch, yaw)
            [roll_X, pitch_Y, yaw_Z] = quat2angle( ...
                [data.(participant_id).Pelvis.(trial).Quat_W, ...
                 data.(participant_id).Pelvis.(trial).Quat_X, ...
                 data.(participant_id).Pelvis.(trial).Quat_Y, ...
                 data.(participant_id).Pelvis.(trial).Quat_Z], ...
                 'XYZ'); % orientation convention based on Cudejko 2022

            eulerAngles = [roll_X, pitch_Y, yaw_Z];

            data.(participant_id).Pelvis.(trial).Euler_X = roll_X;
            data.(participant_id).Pelvis.(trial).Euler_Y = pitch_Y;
            data.(participant_id).Pelvis.(trial).Euler_Z = yaw_Z;
        end

        % Low-pass filter Euler angles
        Fc = 6;
        Wn = Fc / (Fs / 2);
        [b, a] = butter(2, Wn, 'low');
        eulerFilt = filtfilt(b, a, eulerAngles);

        % Mediolateral sway is roll (X rotation) after filtering
        mediolateral_sway = eulerFilt(:, 1);
        results.(participant_id).mediolateral_sway.(trial) = mediolateral_sway;
    end
end

%% =================== Determine analysis window per trial (step-based) ===================
% For each trial, combine L/R step indices, determine which foot starts,
% and define a start/stop index based on selecting steps 3..18 (inclusive).

participants = cell(41,1);
for i = 1:41
    participants{i} = ['P' num2str(i)];
end

trials = {'Ascent1','Ascent2','Descent1','Descent2'};

for p = 2:numel(participants)
    pid = participants{p};

    for q = 1:numel(trials)
        trial = trials{q};

        L = results.(pid).steps.LFoot.(trial);
        R = results.(pid).steps.RFoot.(trial);

        % Merge + sort step events across both feet
        allSteps = sort([L(:); R(:)]);

        % Label steps as L/R based on membership in L (kept as original logic)
        footLabels = repmat('?', size(allSteps));
        footLabels( ismember(allSteps, L) ) = 'L';
        footLabels(~ismember(allSteps, L))  = 'R';

        % Determine which foot stepped first
        if footLabels(1) == 'L'
            first = 'LFoot';
        else
            first = 'RFoot';
        end

        steps.(pid).(trial).allSteps    = allSteps;
        steps.(pid).(trial).footLabels  = footLabels;
        steps.(pid).(trial).firstFoot   = first;

        % Define analysis window using step indices 3..19 (boundaries)
        steps.(pid).(trial).start = steps.(pid).(trial).allSteps(3);
        steps.(pid).(trial).stop  = steps.(pid).(trial).allSteps(19);
    end

    % Debug/example selection (kept as-is)
    start = steps.(pid).Ascent1.start;
    stop  = steps.(pid).Ascent1.stop;
end

%% ============================= Combine trials =============================
% Combine trials into larger segments (implementation inside combine_trials.m).
[data, steps] = combine_trials(data, steps, participants);  % no blending

%% ===== Global-to-local rotation + gravity removal using quaternions =====
% removeGravityWithQuat() creates Acc_noG_* fields, used below as Acc_*.
for i_p = 2:41
    for i_t = 1:4
        participant = participants{i_p};
        trial       = trials{i_t};

        sensor = 'Pelvis';
        data = removeGravityWithQuat(data, participant, sensor, trial, 'SensorToWorldQuat', true);

        sensor = 'LFoot';
        data = removeGravityWithQuat(data, participant, sensor, trial, 'SensorToWorldQuat', true);

        sensor = 'RFoot';
        data = removeGravityWithQuat(data, participant, sensor, trial, 'SensorToWorldQuat', true);
    end
end

%% ===== Replace Acc_* with gravity-removed acceleration (Acc_noG_*) =====
% From this point forward, Acc_X/Y/Z correspond to "no gravity" acceleration.
for p = 2:41
    pid = participants{p};

    for q = 1:numel(trials)
        trial = trials{q};

        for i_s = 1:3
            sensor = sensors{i_s};

            data.(pid).(sensor).(trial).Acc_X = data.(pid).(sensor).(trial).Acc_noG_X;
            data.(pid).(sensor).(trial).Acc_Y = data.(pid).(sensor).(trial).Acc_noG_Y;
            data.(pid).(sensor).(trial).Acc_Z = data.(pid).(sensor).(trial).Acc_noG_Z;
        end
    end
end

%% ============================= Start of analysis =============================
% Metrics computed below:
% - Sample Entropy (SampEn)
% - Harmonic Ratio (HR) using windowed FFT and step frequency estimate
% - Linear features (later)
% - RQA (later)

%% ------------------------ Sample Entropy filter ------------------------
Fc = 15;
Wn = Fc / (Fs / 2);
[b_se, a_se] = butter(2, Wn, 'low');

%% ------------------------ Harmonic Ratio filter ------------------------
Fc = 15;
Wn = Fc / (Fs / 2);
[b_hr, a_hr] = butter(2, Wn, 'low');

directions = {'x','y','z'};

for p = 2:41
    pid = participants{p};

    for q = 1:numel(trials)
        trial = trials{q};

        % Define the index window T for the analysis segment
        if ismember(q, [1 2 3 4])
            T = steps.(pid).(trial).start : steps.(pid).(trial).stop;
        elseif ismember(q, [5 6])
            % For combined trials, use full available step range (kept as-is)
            T = steps.(pid).(trial).allSteps(1) : steps.(pid).(trial).allSteps(end);
        end

        for i = 1:3
            Dir = directions{i};

            % Select pelvis acceleration component based on axis
            if i == 1
                input = data.(pid).Pelvis.(trial).Acc_X;
            elseif i == 2
                input = data.(pid).Pelvis.(trial).Acc_Y;
            elseif i == 3
                input = data.(pid).Pelvis.(trial).Acc_Z;
            end

            %% ===================== Sample Entropy (SampEn) =====================
            filt_input_se = filtfilt(b_se, a_se, input);

            % SampEn parameters
            m = 2;      % embedding dimension
            r = 0.2;    % tolerance (note: standardization handled inside sampen() in your setup)

            [SE, C, M, R] = sampen(filt_input_se(T), m, r, 1, 0);

            %% ====================== Harmonic Ratio (HR) ======================
            % HR is computed using the amplitude of harmonic components derived from FFT.
            % Step frequency is estimated from detected steps; fundamental is f0 = stepfreq/2.
            % For ML axis (y) you invert odd/even definition (kept as original).

            f_search_min = 0.5;  %#ok<NASGU> % reserved for alternative method (unused)
            f_search_max = 3.0;  %#ok<NASGU> % reserved for alternative method (unused)

            filt_input_hr = filtfilt(b_hr, a_hr, input);

            seg = filt_input_hr(T);
            numHarmonics = 10;

            N = length(seg);
            w = hann(N);

            seg   = seg(:) - mean(seg);
            seg_w = seg .* w;

            % Zero-padding for finer FFT resolution
            Nfft = 8 * N;

            Y = fft(seg_w, Nfft);
            A = abs(Y)/N;

            % Window amplitude compensation
            A = A / (sum(w)/N);

            % Single-sided spectrum
            halfIdx = floor(Nfft/2) + 1;
            A_single = A(1:halfIdx);
            A_single(2:end-(mod(Nfft,2)==0)) = 2 * A_single(2:end-(mod(Nfft,2)==0));

            f = (0:halfIdx-1) * (Fs / Nfft);

            % Step frequency estimate based on steps 3..19 (inclusive) for this trial
            trialsteps = steps.(pid).(trial).allSteps(3:19) / Fs;
            trialtime  = trialsteps(end) - trialsteps(1);
            f_step     = (length(trialsteps) - 1) / trialtime;
            f0         = f_step/2;

            results.(pid).Linear.(trial).stepfreq = f_step;

            % Extract harmonic amplitudes at multiples of f0
            harmonicFreqs = f0 * (1:numHarmonics);
            harmonicAmp   = nan(numHarmonics,1);

            for k = 1:numHarmonics
                if harmonicFreqs(k) >= f(end), continue; end
                [~, idx] = min(abs(f - harmonicFreqs(k)));
                harmonicAmp(k) = A_single(idx);
            end

            evenSum = sum(harmonicAmp(2:2:end));
            oddSum  = sum(harmonicAmp(1:2:end));

            if evenSum == 0 || oddSum == 0
                HR_evenodd = NaN;
            else
                % Axis-specific convention: y uses odd/even inversion (kept as original)
                if strcmpi(Dir,'y')
                    HR_evenodd = oddSum / evenSum;
                else
                    HR_evenodd = evenSum / oddSum;
                end
            end

            %% ============================= Store results =============================
            results.(pid).SampEn.(trial).(Dir) = SE;
            results.(pid).HR.(trial).(Dir)     = HR_evenodd;
        end
    end
end

%% ================== TISEAN + CRP Toolbox RQA block ==================
% Requirements:
% - TISEAN 3.x Windows binaries available locally (mutual.exe required here)
% - Marwan CRP Toolbox on MATLAB path
%
% Pipeline:
% 1) Estimate delay tau via AMI (mutual information) using TISEAN mutual.exe
% 2) Estimate embedding dimension m via CRP Toolbox fnn()
% 3) Choose epsilon based on fraction of max distance (or fixed RR if enabled)
% 4) Compute RQA metrics (RR, DET, LAM) + store RP matrix

crpRoot = 'C:\Users\boris\OneDrive\Documenten\M Human Movement Sciences\MRP\data\CRPtool-master';
assert(isfolder(crpRoot), 'CRP toolbox folder not found: %s', crpRoot);
addpath(genpath(crpRoot));
rehash toolboxcache;

% Parameter log: {pid, trial, axis, tau, m, epsilon, RR, DET, LAM}
paramRows = {};

tiseanBin = 'C:\Users\boris\OneDrive\Documenten\M Human Movement Sciences\MRP\data\Tisean_3.0.0\bin';
mutualExe = fullfile(tiseanBin,'mutual.exe');
assert(isfile(mutualExe), 'mutual.exe not found: %s', mutualExe);

% Quick check: executable responds
[stH1,~] = system(['"',mutualExe,'" -h']);
assert(stH1==0, 'mutual.exe is not executable');

outDir = fullfile(pwd, 'RQA_data');
if ~exist(outDir, 'dir'), mkdir(outDir); end

maxDelay = 40;                      % AMI: maximum delay (samples)
datatype = {'Acc_X','Acc_Y','Acc_Z'};
EPS_FRAC = 0.1;                     % epsilon = 10% of max pairwise distance

% CRP Toolbox settings
CRP_LMIN     = 2;                   %#ok<NASGU>
CRP_TTHEILER = 3;                   % overwritten later: (m-1)*tau
CRP_NORM     = 'euclidean';         %#ok<NASGU>

% FNN settings (CRP toolbox)
FNN_MAX_M   = 8;
FNN_R       = 15;
FNN_S       = 5;
FNN_MAXN    = 4000;
FNN_METHOD  = 'euclidean';
FNN_THRESH  = 0.01;                 % threshold for choosing m (kept as original)
FNN_SEED    = 42;

% FNN caching
FNN_CACHE_DIR = fullfile(outDir,'cache_fnn');
if ~exist(FNN_CACHE_DIR,'dir'), mkdir(FNN_CACHE_DIR); end

% Plot flags
MAKE_RP_PLOTS             = false;   % will create binary recurrence plots
MAKE_COLOR_DISTANCE_PLOTS = false;   % will create heatmap-like euclidian distance plots

Use_Fixed_tau = true;     % fixed tau is used in MRP, fixed recurrence rate is also interesting however.
fixedRR       = false;

RPpngDir = fullfile(outDir,'RP_png');
if MAKE_RP_PLOTS && ~exist(RPpngDir,'dir'), mkdir(RPpngDir); end

DistPngDir = fullfile(outDir,'Distance_png');
if MAKE_COLOR_DISTANCE_PLOTS && ~exist(DistPngDir,'dir'), mkdir(DistPngDir); end

for p = 2:10  % optionally extend to 2:41
    pid = participants{p};

    for q = 1:4
        trial = trials{q};

        if ismember(q, [1 2 3 4])
            T = steps.(pid).(trial).start : steps.(pid).(trial).stop;
        else
            T = steps.(pid).(trial).allSteps(1) : steps.(pid).(trial).allSteps(end);
        end

        for iD = 1:3
            axisName = datatype{iD};

            % Extract and crop data segment
            x = data.(pid).Pelvis.(trial).(axisName);
            x = x(T);
            x = x(~isnan(x));

            % Normalize x (z-score)
            mu  = mean(x);
            sig = std(x);
            if ~isfinite(sig) || sig == 0
                warning('%s/%s/%s: std=0 or NaN, skipping.', pid, trial, axisName);
                continue
            end
            x = (x - mu) ./ sig;

            if numel(x) < (maxDelay+5)
                warning('Too few samples for %s/%s/%s, skipping.', pid, trial, axisName);
                continue
            end

            % AMI I/O files
            base    = sprintf('%s_%s_%s', pid, trial, axisName);
            inFile  = fullfile(outDir, [base '_in.txt']);
            amiFile = fullfile(outDir, [base '_ami.txt']);
            writematrix(x, inFile, 'Delimiter','\t');

            % Run TISEAN mutual.exe to obtain AMI curve
            cmdMut = sprintf('"%s" -c 1 -D %d -o "%s" "%s"', mutualExe, maxDelay, amiFile, inFile);
            [st1, out1] = system(cmdMut);

            dA = dir(amiFile);
            if st1 ~= 0 || isempty(dA) || dA.bytes == 0
                error('mutual.exe failed for %s\nCMD: %s\nOUT:\n%s', base, cmdMut, out1);
            end

            % Choose tau as first local minimum in AMI
            tau = pick_tau_from_ami(amiFile);

            % Optional override: fixed tau for all conditions
            if Use_Fixed_tau
                tau = 6;
            end

            % FNN-based embedding dimension selection with caching
            xseg = x(:);
            rng(FNN_SEED);

            fnnKey = sprintf('%s_%s_%s_tau%d_mmax%d_r%d_s%g_n%d_%s.mat', ...
                pid, trial, axisName, tau, FNN_MAX_M, FNN_R, FNN_S, FNN_MAXN, FNN_METHOD);
            fnnCacheFile = fullfile(FNN_CACHE_DIR, fnnKey);

            if exist(fnnCacheFile,'file')
                S = load(fnnCacheFile);
                FNNvec = S.FNNvec;
                mOpt   = S.mOpt;
            else
                FNNvec = fnn(xseg, FNN_MAX_M, tau, FNN_R, FNN_S, FNN_MAXN, ...
                             FNN_METHOD, 'nogui','silent');
                mOpt = pick_m_from_crp_fnn(FNNvec, FNN_THRESH);
                save(fnnCacheFile, 'FNNvec','mOpt','tau','FNN_MAX_M','FNN_R','FNN_S','FNN_MAXN','FNN_METHOD','FNN_THRESH');
            end

            fprintf('%s | %s | %s  -> tau = %d (AMI), m = %d (FNN)\n', pid, trial, axisName, tau, mOpt);

            % Embed and compute distances to set epsilon
            N1 = numel(x) - (mOpt-1)*tau;
            if N1 < 2
                warning('%s/%s/%s: too short after embedding (N1=%d), skipping.', pid, trial, axisName, N1);
                continue
            end

            X = zeros(N1, mOpt);
            for ii = 1:mOpt
                X(:,ii) = x((1:N1) + (ii-1)*tau);
            end

            dist = pdist(X,'euclidean');
            if isempty(dist)
                warning('%s/%s/%s: no pairwise distances, skipping.', pid, trial, axisName);
                continue
            end

            epsFixed = EPS_FRAC * max(dist);

            if fixedRR
                TARGET_RR_PCT = 1;
                CRP_TTHEILER  = (mOpt - 1) * tau;
                epsFixed = eps_for_fixed_rr(xseg, mOpt, tau, TARGET_RR_PCT, CRP_TTHEILER, true);

                if isnan(epsFixed)
                    warning('%s/%s/%s: could not determine epsilon for fixed RR', pid, trial, axisName);
                    continue
                end
            end

            % Theiler window: exclude temporally close neighbors
            CRP_TTHEILER = (mOpt - 1) * tau;

            % Compute RQA metrics via CRQA
            Y = crqa(x, mOpt, tau, epsFixed, [], 1, 2, 2, CRP_TTHEILER, ...
                     'euclidean','nonormalize','nogui','silent');

            RR  = Y(1,1);
            DET = Y(1,2);
            LAM = Y(1,6);

            % Compute recurrence plot matrix (binary)
            R = crp(x, mOpt, tau, epsFixed, 'euclidean','nonormalize','nogui','silent');

            % Store results
            results.(pid).RQA.(trial).(axisName).tau = tau;
            results.(pid).RQA.(trial).(axisName).m   = mOpt;
            results.(pid).RQA.(trial).(axisName).RR  = RR;
            results.(pid).RQA.(trial).(axisName).DET = DET;
            results.(pid).RQA.(trial).(axisName).LAM = LAM;
            results.(pid).RQA.(trial).(axisName).RP  = R;

            results.(pid).RQA.(trial).(axisName).FNN = struct( ...
                'vec',    FNNvec, ...
                'tau',    tau, ...
                'mOpt',   mOpt, ...
                'params', struct('R',FNN_R,'S',FNN_S,'N',FNN_MAXN,'method',FNN_METHOD,'thresh',FNN_THRESH));

            paramRows(end+1,:) = {pid, trial, axisName, tau, mOpt, epsFixed, RR, DET, LAM};

            % Optional: export recurrence plot
            if MAKE_RP_PLOTS
                fig = figure('Visible','on');
                imagesc(R);
                axis image off;
                colormap(flipud(gray));
                set(gca, 'YDir', 'normal'); axis image off;

                title(sprintf('%s | %s | %s  (RR=%.2f%%, eps=%.3g, m=%d, \\tau=%d)', ...
                    pid, trial, axisName, 100*RR, epsFixed, mOpt, tau), 'Interpreter','none');

                exportgraphics(fig, fullfile(RPpngDir, sprintf('%s_%s_%s_RP.png', pid, trial, axisName)), 'Resolution',150);
            end

            % Optional: export distance heatmap (no thresholding)
            if MAKE_COLOR_DISTANCE_PLOTS
                D = squareform(dist);

                fig2 = figure('Visible','on');
                imagesc(D);
                axis image off;
                colormap(turbo);
                colorbar;
                set(gca, 'YDir', 'normal'); axis image off;

                title(sprintf('%s | %s | %s  Distance heatmap (m=%d, \\tau=%d)', ...
                    pid, trial, axisName, mOpt, tau), 'Interpreter','none');

                exportgraphics(fig2, fullfile(DistPngDir, sprintf('%s_%s_%s_DIST.png', pid, trial, axisName)), 'Resolution',150);
            end
        end
    end
end

% Summarize and save parameter table
ParamsTbl = cell2table(paramRows, ...
    'VariableNames', {'Participant','Trial','Axis','Tau','M','Epsilon','RR','DET','LAM'});
disp(ParamsTbl);
writetable(ParamsTbl, fullfile(outDir,'rqa_params_with_rr.csv'));

if exist('groupsummary','file') == 2
    SummaryTbl = groupsummary(ParamsTbl, {'Participant','Trial'}, 'median', {'Tau','M'});
    disp(SummaryTbl);
    writetable(SummaryTbl, 'tisean_params_summary.csv');
end

%% ===================== Linear measures =====================
% Step frequency is already stored in the HR block.
% Here we compute:
% - peak vertical foot acceleration (starting foot)
% - peak pelvis acceleration per axis
% - RMS pelvis acceleration per axis
% - step frequency variability (std)

low = 20 / (Fs/2);
[b, a] = butter(2, low, 'low');

for p = 2:41
    pid = participants{p};

    for q = 1:4
        trial = trials{q};

        % Use the detected starting foot for foot-based vertical peak
        sensor = steps.(pid).(trial).firstFoot;

        acc_voetvt = data.(pid).(sensor).(trial).Acc_Z;

        acc_pelvx = data.(pid).Pelvis.(trial).Acc_X;
        acc_pelvy = data.(pid).Pelvis.(trial).Acc_Y;
        acc_pelvz = data.(pid).Pelvis.(trial).Acc_Z;

        idx = steps.(pid).(trial).start : steps.(pid).(trial).stop;

        acc_voetvt_f = filtfilt(b, a, acc_voetvt(idx));
        acc_pelvx_f  = filtfilt(b, a, acc_pelvx(idx));
        acc_pelvy_f  = filtfilt(b, a, acc_pelvy(idx));
        acc_pelvz_f  = filtfilt(b, a, acc_pelvz(idx));

        peak_acc_footVT = max(abs(acc_voetvt_f));
        peak_acc_pelvx  = max(abs(acc_pelvx_f));
        peak_acc_pelvy  = max(abs(acc_pelvy_f));
        peak_acc_pelvz  = max(abs(acc_pelvz_f));

        RMSacc_pelvx = sqrt(mean(acc_pelvx_f.^2));
        RMSacc_pelvy = sqrt(mean(acc_pelvy_f.^2));
        RMSacc_pelvz = sqrt(mean(acc_pelvz_f.^2));

        % Step interval variability (based on steps 3..18)
        dt = diff(steps.(pid).(trial).allSteps(3:18));  % step intervals (samples)
        f  = 1 ./ dt;                                   % step frequency proxy (1/sample)

        mean_f = mean(f); %#ok<NASGU>
        std_f  = std(f);

        results.(pid).Linear.(trial).peakfootacc  = peak_acc_footVT;
        results.(pid).Linear.(trial).peakpelacc_x = peak_acc_pelvx;
        results.(pid).Linear.(trial).peakpelacc_y = peak_acc_pelvy;
        results.(pid).Linear.(trial).peakpelacc_z = peak_acc_pelvz;

        results.(pid).Linear.(trial).RMSpelacc_x  = RMSacc_pelvx;
        results.(pid).Linear.(trial).RMSpelacc_y  = RMSacc_pelvy;
        results.(pid).Linear.(trial).RMSpelacc_z  = RMSacc_pelvz;

        results.(pid).Linear.(trial).std_stepf    = std_f;
    end
end

%% ===== Add averaged trials: meanascent (A1+A2)/2, meandescent (D1+D2)/2 =====
participants = fieldnames(results);

axes_xyz = {'x','y','z'};
axes_rqa = {'Acc_X','Acc_Y','Acc_Z'};

nl_measures_xyz    = {'SampEn','HR'};                     % axis-based nonlinear measures
rqa_metrics        = {'RR','DET','LAM'};                  % RQA metrics per embedded axis
lin_noaxis         = {'peakfootacc','stepfreq','std_stepf'};
lin_axes_basenames = {'peakpelacc','RMSpelacc'};

avg2 = @(a,b) mean([a,b],'omitnan');  % safe 2-point mean

for p = 1:numel(participants)
    pid = participants{p};

    trAsc = {'Ascent1','Ascent2'};
    trDes = {'Descent1','Descent2'};

    % --- Nonlinear: SampEn and HR (x/y/z) ---
    for k = 1:numel(nl_measures_xyz)
        meas = nl_measures_xyz{k};

        for a = 1:numel(axes_xyz)
            ax = axes_xyz{a};

            a1 = tryget(results, {pid,meas,trAsc{1},ax}, NaN);
            a2 = tryget(results, {pid,meas,trAsc{2},ax}, NaN);
            if ~(isnan(a1) && isnan(a2))
                results.(pid).(meas).meanascent.(ax) = avg2(a1,a2);
            end

            d1 = tryget(results, {pid,meas,trDes{1},ax}, NaN);
            d2 = tryget(results, {pid,meas,trDes{2},ax}, NaN);
            if ~(isnan(d1) && isnan(d2))
                results.(pid).(meas).meandescent.(ax) = avg2(d1,d2);
            end
        end
    end

    % --- RQA: RR/DET/LAM per Acc_X/Y/Z ---
    for a = 1:numel(axes_rqa)
        axfld = axes_rqa{a};

        for r = 1:numel(rqa_metrics)
            mname = rqa_metrics{r};

            a1 = tryget(results, {pid,'RQA',trAsc{1},axfld,mname}, NaN);
            a2 = tryget(results, {pid,'RQA',trAsc{2},axfld,mname}, NaN);
            if ~(isnan(a1) && isnan(a2))
                results.(pid).RQA.meanascent.(axfld).(mname) = avg2(a1,a2);
            end

            d1 = tryget(results, {pid,'RQA',trDes{1},axfld,mname}, NaN);
            d2 = tryget(results, {pid,'RQA',trDes{2},axfld,mname}, NaN);
            if ~(isnan(d1) && isnan(d2))
                results.(pid).RQA.meandescent.(axfld).(mname) = avg2(d1,d2);
            end
        end
    end

    % --- Linear: non-axis ---
    for k = 1:numel(lin_noaxis)
        meas = lin_noaxis{k};

        a1 = tryget(results, {pid,'Linear',trAsc{1},meas}, NaN);
        a2 = tryget(results, {pid,'Linear',trAsc{2},meas}, NaN);
        if ~(isnan(a1) && isnan(a2))
            results.(pid).Linear.meanascent.(meas) = avg2(a1,a2);
        end

        d1 = tryget(results, {pid,'Linear',trDes{1},meas}, NaN);
        d2 = tryget(results, {pid,'Linear',trDes{2},meas}, NaN);
        if ~(isnan(d1) && isnan(d2))
            results.(pid).Linear.meandescent.(meas) = avg2(d1,d2);
        end
    end

    % --- Linear: axis-based fields stored as <basename>_<x|y|z> ---
    for b = 1:numel(lin_axes_basenames)
        base = lin_axes_basenames{b};

        for a = 1:numel(axes_xyz)
            ax  = axes_xyz{a};
            fld = sprintf('%s_%s', base, ax);

            a1 = tryget(results, {pid,'Linear',trAsc{1},fld}, NaN);
            a2 = tryget(results, {pid,'Linear',trAsc{2},fld}, NaN);
            if ~(isnan(a1) && isnan(a2))
                results.(pid).Linear.meanascent.(fld) = avg2(a1,a2);
            end

            d1 = tryget(results, {pid,'Linear',trDes{1},fld}, NaN);
            d2 = tryget(results, {pid,'Linear',trDes{2},fld}, NaN);
            if ~(isnan(d1) && isnan(d2))
                results.(pid).Linear.meandescent.(fld) = avg2(d1,d2);
            end
        end
    end
end

%% ===================== Build SPSS long-format table =====================
% Exports all measures into a uniform long table:
% Columns: Participant, AgeGroup, Trial, Axis, Measure, Value
%
% Axis:
% - 'x','y','z' for axis-based measures
% - 'No' for non-axis measures (e.g., stepfreq, peakfootacc)

participants = fieldnames(results);
save(results)  % NOTE: kept as-is, although typical MATLAB usage is save('results.mat','results')

ages_P2P41 = [25,23,70,74,79,24,26,23,75,67,68,69,20,67,22,23,25,75,26,71,67,24,77,68,25,70,77,75,32,24,25,71,75,23,66,69,75,24,26,66];
participants_P2P41 = arrayfun(@(x) sprintf('P%d',x), 2:41, 'UniformOutput', false);

young_group = participants_P2P41(ages_P2P41 < 35);
old_group   = participants_P2P41(ages_P2P41 >= 65);

trials_export = {'Ascent1','Descent1','Ascent2','Descent2','meanascent','meandescent'};

Rows = {};  % {Participant, AgeGroup, Trial, Axis, Measure, Value}

% Map RQA axis fieldnames to x/y/z naming for SPSS compatibility
axis_map_rqa = struct('Acc_X','x','Acc_Y','y','Acc_Z','z');

for p = 1:numel(participants_P2P41)
    pid = participants_P2P41{p};

    if ismember(pid, young_group)
        grp = "Young";
    elseif ismember(pid, old_group)
        grp = "Old";
    else
        continue
    end

    for it = 1:numel(trials_export)
        tr = trials_export{it};

        % --- SampEn & HR (x,y,z) ---
        for ax = ["x","y","z"]
            SE  = tryget(results, {pid,'SampEn',tr,char(ax)}, NaN);
            HRv = tryget(results, {pid,'HR',    tr,char(ax)}, NaN);

            Rows = add_row(Rows, pid, grp, tr, char(ax), 'SampEn', SE);
            Rows = add_row(Rows, pid, grp, tr, char(ax), 'HR',    HRv);
        end

        % --- Linear (axis-based pelvis metrics) ---
        for ax = ["x","y","z"]
            peakpel = tryget(results, {pid,'Linear',tr,['peakpelacc_' char(ax)]}, NaN);
            rmspel  = tryget(results, {pid,'Linear',tr,['RMSpelacc_'  char(ax)]}, NaN);

            Rows = add_row(Rows, pid, grp, tr, char(ax), 'peakpelacc', peakpel);
            Rows = add_row(Rows, pid, grp, tr, char(ax), 'RMSpelacc',  rmspel);
        end

        % --- Linear (non-axis measures) ---
        peakfoot = tryget(results, {pid,'Linear',tr,'peakfootacc'}, NaN);
        stepfq   = tryget(results, {pid,'Linear',tr,'stepfreq'},    NaN);
        stdfq    = tryget(results, {pid,'Linear',tr,'std_stepf'},   NaN);

        Rows = add_row(Rows, pid, grp, tr, 'No', 'peakfootacc', peakfoot);
        Rows = add_row(Rows, pid, grp, tr, 'No', 'stepfreq',    stepfq);
        Rows = add_row(Rows, pid, grp, tr, 'No', 'std_stepf',   stdfq);

        % --- RQA (RR/DET/LAM) ---
        for fld = fieldnames(axis_map_rqa)'
            dirFld = fld{1};                % 'Acc_X' etc.
            ax     = axis_map_rqa.(dirFld); % 'x'|'y'|'z'

            RR  = tryget(results, {pid,'RQA',tr,dirFld,'RR'},  NaN);
            DET = tryget(results, {pid,'RQA',tr,dirFld,'DET'}, NaN);
            LAM = tryget(results, {pid,'RQA',tr,dirFld,'LAM'}, NaN);

            Rows = add_row(Rows, pid, grp, tr, ax, 'RQA_RR',  RR);
            Rows = add_row(Rows, pid, grp, tr, ax, 'RQA_DET', DET);
            Rows = add_row(Rows, pid, grp, tr, ax, 'RQA_LAM', LAM);
        end
    end
end

SPSS_Long = cell2table(Rows, 'VariableNames', ...
    {'Participant','AgeGroup','Trial','Axis','Measure','Value'});

out_csv = fullfile(pwd, 'SPSS_alldata_uniform.csv');
writetable(SPSS_Long, out_csv);
disp(SPSS_Long(1:min(20,height(SPSS_Long)),:));
fprintf('CSV written: %s\n', out_csv);

%% ===================== Helper functions (must be at end of script) =====================
function val = tryget(S, flds, default)
%TRYGET Safely read nested struct fields without hard failures.
%   Example:
%     val = tryget(results, {pid,'SampEn',tr,'x'}, NaN)
%
% Returns:
%   - the found value (scalar if numeric)
%   - default if any field is missing or indexing fails

    val = default;
    try
        T = S;
        for k = 1:numel(flds)
            f = flds{k};
            if ~isfield(T, f), return; end
            T = T.(f);
        end

        if ~isempty(T)
            if isnumeric(T)
                val = T(1);
            else
                val = T;
            end
        end
    catch
        val = default;
    end
end

function Rows = add_row(Rows, pid, grp, tr, ax, me, val)
%ADD_ROW Append a single row to the long-format cell array.
%   Rows = add_row(Rows, pid, grp, tr, ax, me, val)

    Rows(end+1,:) = {pid, grp, tr, ax, me, val};
end

function v = fetch_val(R, pid, tr, meas, dir)
%FETCH_VAL Retrieve scalar values from results struct given measure type.
% Supported patterns:
% - SampEn / HR: results.(pid).(meas).(tr).(x|y|z)
% - RQA_*:       results.(pid).RQA.(tr).(Acc_X|Acc_Y|Acc_Z).(RR|DET|LAM)

    v = NaN;
    try
        if startsWith(meas,'RQA_')
            m = extractAfter(meas,'RQA_');

            if isfield(R.(pid),'RQA') && isfield(R.(pid).RQA,tr) && ...
               isfield(R.(pid).RQA.(tr),dir) && isfield(R.(pid).RQA.(tr).(dir), m)
                vv = R.(pid).RQA.(tr).(dir).(m);
                if ~isempty(vv), v = vv(1); end
            end
        else
            if isfield(R.(pid),meas) && isfield(R.(pid).(meas),tr) && ...
               isfield(R.(pid).(meas).(tr),dir)
                vv = R.(pid).(meas).(tr).(dir);
                if ~isempty(vv), v = vv(1); end
            end
        end
    catch
        v = NaN;
    end
end

%% ===================== Helper functions for RQA =====================
function tau = pick_tau_from_ami(amiFile)
%PICK_TAU_FROM_AMI Choose tau from an AMI curve (TISEAN output).
% Strategy:
% - Prefer the first local minimum of AMI
% - If none exists, use the global minimum (excluding delay=0)

    T = readtable(amiFile, 'FileType','text', 'ReadVariableNames',false);
    A = T{:,1:2};

    delays  = A(:,1);
    amiVals = A(:,2);

    dV = diff(amiVals);

    % Detect first local minimum by sign change in derivative
    locMin = find([false; dV(1:end-1) < 0 & dV(2:end) > 0; false], 1, 'first');

    if isempty(locMin)
        [~, locMin] = min(amiVals(2:end));
        locMin = locMin + 1;
    end

    tau = max(1, delays(locMin));
end

function mOpt = pick_m_from_crp_fnn(FNNvec, thresh)
%PICK_M_FROM_CRP_FNN Select embedding dimension m from a FNN curve.
% Rule:
% - choose the first m with FNN <= thresh
% - if never below thresh, choose m at minimum FNN
% - enforce m >= 2

    if nargin < 2 || isempty(thresh)
        thresh = 0.05;
    end

    idx = find(FNNvec <= thresh, 1, 'first');
    if isempty(idx)
        [~, idx] = min(FNNvec);
    end

    mOpt = max(2, idx);
end

function eps_opt = eps_for_fixed_rr(x, m, tau, targetRR_pct, theiler, doZscore)
%EPS_FOR_FIXED_RR Choose epsilon such that recurrence rate approximates targetRR_pct.
% This uses the distance distribution percentile after applying Theiler exclusion.

    if nargin < 6, doZscore = true; end
    x = x(:);

    N1 = numel(x) - (m-1)*tau;
    if N1 < 3, eps_opt = NaN; return; end

    % Embed
    Y = zeros(N1, m);
    for ii = 1:m
        Y(:,ii) = x((1:N1) + (ii-1)*tau);
    end

    if doZscore
        Y = (Y - mean(Y,1)) ./ std(Y,0,1);
    end

    % Pairwise distances and Theiler mask
    D = squareform(pdist(Y, 'euclidean'));

    mask = true(N1);
    mask(1:N1+1:end) = false;

    if theiler > 0
        for k = -theiler:theiler
            if k == 0, continue; end
            M = eye(N1, N1);
            M = circshift(M, k, 2);
            mask(M==1) = false;
        end
    end

    dvec = D(mask);
    if isempty(dvec), eps_opt = NaN; return; end

    eps_opt = quantile(dvec, targetRR_pct/100);
end

function eps_opt = eps_fixed_rr_bisect(x, m, tau, rr_target_pct, theiler)
%EPS_FIXED_RR_BISECT Alternative epsilon selection using bisection + crqa() calls.
% Kept as an unused helper for experimentation.

    x = (x(:)-mean(x))/std(x);

    N1 = numel(x)-(m-1)*tau;
    Y  = zeros(N1,m);
    for ii = 1:m
        Y(:,ii)=x((1:N1)+(ii-1)*tau);
    end

    d = pdist(Y,'euclidean');
    d = sort(d);

    lo = d(max(1, round(0.001*numel(d))));
    hi = d(min(numel(d), round(0.999*numel(d))));

    target = rr_target_pct;
    tol    = 0.05;

    for it = 1:25
        mid = 0.5*(lo+hi);

        M = crqa(x, m, tau, mid, [], 1, 2, 2, theiler, ...
                 'euclidean','normalize','nogui','silent');

        RR = M(1,1);

        if RR > target
            hi = mid;
        else
            lo = mid;
        end

        if abs(RR - target) <= tol
            break
        end
    end

    eps_opt = 0.5*(lo+hi);
end
