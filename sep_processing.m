%% SEP COMPLETE ANALYSIS - PROCESSING + BALANCED SELECTION + VISUALIZATION

clear all; close all; clc;

%% CONFIGURATION

PARTICIPANTS = {
    %'P001'
    %'P002'
    'P003'
    'P004'
    'P005'
    'P006'
    'P007'
    'P008'
    'P009'
    'P010'
};

REFERENCE_METHOD    = 'Cz';
N_EPOCHS_PER_SPEED  = 100;

%% File locations
CSV_FOLDER = 'data\csv_files';
MAT_FOLDER = 'data\mat_files';

if ~exist(CSV_FOLDER, 'dir') || ~exist(MAT_FOLDER, 'dir')
    error('Data folders not found');
end

%% Analysis parameters
FS = 1200;

% Trigger detection
REF_THRESHOLD          = 0.5;
MIC_PERCENTILE_LOW     = 10;
MIC_PERCENTILE_HIGH    = 90;
MIC_THRESHOLD_FRACTION = 0.3;
MIN_ISI_MS             = 500;
EXPECTED_ISI_MS        = 2000;
REF_SEARCH_BEFORE_MS   = 200;
REF_SEARCH_AFTER_MS    = 400;

% Filtering
FILTER_ORDER = 1000;
FILTER_LOW   = 2;
FILTER_HIGH  = 32;

% Epoching
EPOCH_START    = -0.2;
EPOCH_END      =  0.8;
BASELINE_START = -0.2;
BASELINE_END   =  0;

% Component windows (seconds)
P1_START  = 0.040;   P1_END  = 0.090;
N1_START  = 0.100;   N1_END  = 0.160;
P2_START  = 0.180;   P2_END  = 0.235;
P1N1_START = P1_START;  P1N1_END = N1_END;
N1P2_START = N1_START;  N1P2_END = P2_END;

% Artifact rejection
ARTIFACT_THRESHOLD = 20;   % µV

% Microphone delay
MICROPHONE_DELAY_MS      = 55;
MICROPHONE_DELAY_SAMPLES = round(MICROPHONE_DELAY_MS * FS / 1000);

 % Electrode list defined ONCE here, used everywhere below
ELECTRODES = {'FC3','C3','CP3','FC4','C4','CP4'};
N_ELECTRODES = length(ELECTRODES);

%% Speed-condition colours (used in all figures)
SPEED_COLORS = [
    0.00, 0.45, 0.74;   % Blue   – 1.0 m/s
    0.85, 0.33, 0.10;   % Orange – 3.0 m/s
    0.93, 0.69, 0.13;   % Yellow – 5.0 m/s
    0.49, 0.18, 0.56;   % Purple – 7.0 m/s
    0.47, 0.67, 0.19;   % Green  – 10.0 m/s
];

%% Derived parameters
epoch_samples = round(EPOCH_START*FS):round(EPOCH_END*FS);
n_samples     = length(epoch_samples);
time_vector   = epoch_samples / FS;
time_ms       = time_vector * 1000;   % milliseconds – defined once here
baseline_idx  = time_vector >= BASELINE_START & time_vector < BASELINE_END;

P1_idx   = time_vector >= P1_START   & time_vector <= P1_END;
N1_idx   = time_vector >= N1_START   & time_vector <= N1_END;
P2_idx   = time_vector >= P2_START   & time_vector <= P2_END;
P1N1_idx = time_vector >= P1N1_START & time_vector <= P1N1_END;
N1P2_idx = time_vector >= N1P2_START & time_vector <= N1P2_END;

b = fir1(FILTER_ORDER, [FILTER_LOW, FILTER_HIGH]/(FS/2), 'bandpass');

fprintf('  Participants:  %d\n',  length(PARTICIPANTS));
fprintf('  Reference:     %s\n',  REFERENCE_METHOD);
fprintf('  Sampling rate: %d Hz\n', FS);
fprintf('  Filter:        %d–%d Hz\n', FILTER_LOW, FILTER_HIGH);
fprintf('  Artifact thr:  %d µV\n', ARTIFACT_THRESHOLD);
fprintf('  Mic delay:     %d ms (%d samples)\n', MICROPHONE_DELAY_MS, MICROPHONE_DELAY_SAMPLES);
fprintf('  Balanced sel:  %d epochs/speed/participant\n\n', N_EPOCHS_PER_SPEED);

%% PART 1: PROCESS ALL PARTICIPANTS

GROUP_C3  = [];  GROUP_C4  = [];
GROUP_CP3 = [];  GROUP_CP4 = [];
GROUP_FC3 = [];  GROUP_FC4 = [];
GROUP_trial_info         = [];
GROUP_participant_labels = [];

PARTICIPANT_AVERAGES   = struct();
successful_participants = {};
failed_participants     = {};

for p = 1:length(PARTICIPANTS)
    PARTICIPANT_ID = PARTICIPANTS{p};

    fprintf('Processing %d/%d: %s\n', p, length(PARTICIPANTS), PARTICIPANT_ID);

    try
        %% Load CSV
        csv_file = fullfile(CSV_FOLDER, sprintf('%s_trials.csv', PARTICIPANT_ID));
        if ~exist(csv_file, 'file')
            warning('CSV file not found for %s', PARTICIPANT_ID);
            failed_participants{end+1} = PARTICIPANT_ID;
            continue;
        end
        trial_info = readtable(csv_file);

        %% Find MAT files
        participant_num     = str2double(PARTICIPANT_ID(2:end));
        participant_num_str = sprintf('%02d', participant_num);
        trial_files = dir(fullfile(MAT_FOLDER, sprintf('s%s_tr*_SS.mat', participant_num_str)));

        if isempty(trial_files)
            warning('No MAT files for %s', PARTICIPANT_ID);
            failed_participants{end+1} = PARTICIPANT_ID;
            continue;
        end
        n_blocks   = length(trial_files);
        n_csv_rows = height(trial_info);
        fprintf('  Blocks:    %d\n', n_blocks);
        fprintf('  CSV rows:  %d  (expected ~%d per block = ~%d total)\n', ...
                n_csv_rows, round(n_csv_rows/n_blocks), n_csv_rows);

        DATA_FILES = cell(1, length(trial_files));
        for i = 1:length(trial_files)
            DATA_FILES{i} = fullfile(MAT_FOLDER, trial_files(i).name);
        end

        %% Participant-level epoch storage
        all_epochs_C3  = [];  all_epochs_C4  = [];
        all_epochs_CP3 = [];  all_epochs_CP4 = [];
        all_epochs_FC3 = [];  all_epochs_FC4 = [];
        all_trial_info_clean = [];

 % raw_trial_counter tracks position in trial_info (CSV).
        %% CSV has exactly TRIALS_PER_BLOCK rows per block.
        %% We must always increment by this fixed number – NOT by how many triggers were detected.
        TRIALS_PER_BLOCK = round(height(trial_info) / length(DATA_FILES));
        raw_trial_counter = 0;
        fprintf('  CSV: %d rows / %d blocks = %d trials per block\n', ...
                height(trial_info), length(DATA_FILES), TRIALS_PER_BLOCK);

        for block = 1:length(DATA_FILES)
            %% Load block
            data   = load(DATA_FILES{block});
            fields = fieldnames(data);
            raw_data = data.(fields{1});

            %% Channel extraction
            ref_trigger = raw_data(2,  :);
            mic_trigger = raw_data(4,  :);
            C3          = raw_data(6,  :);
            C4          = raw_data(7,  :);
            CP3         = raw_data(8,  :);
            CP4         = raw_data(9,  :);
            FC3         = raw_data(10, :);
            FC4         = raw_data(11, :);
            right_ear   = raw_data(12, :);
            left_ear    = raw_data(13, :);
            Cz          = raw_data(14, :);

            %% Adaptive MIC-first with unique REF matching
            ref_norm = (ref_trigger - min(ref_trigger)) / ...
                       (max(ref_trigger) - min(ref_trigger));

            mic_baseline = prctile(mic_trigger, MIC_PERCENTILE_LOW);
            mic_high     = prctile(mic_trigger, MIC_PERCENTILE_HIGH);
            mic_range    = mic_high - mic_baseline;

            mic_norm = (mic_trigger - mic_baseline) / mic_range;
            mic_norm(mic_norm < 0) = 0;
            mic_norm(mic_norm > 1) = 1;

            mic_binary    = mic_norm > MIC_THRESHOLD_FRACTION;
            all_mic_rising = find(diff(mic_binary) == 1);

            MIN_ISI_SAMPLES = round(MIN_ISI_MS * FS / 1000);

            mic_triggers = [];
            last_trigger = -inf;
            for i = 1:length(all_mic_rising)
                if (all_mic_rising(i) - last_trigger) >= MIN_ISI_SAMPLES
                    mic_triggers(end+1) = all_mic_rising(i);
                    last_trigger = all_mic_rising(i);
                end
            end

            %% Remove trial starts (gaps > 3× expected ISI)
            GAP_THRESH = round(EXPECTED_ISI_MS * 3 * FS / 1000);
            if length(mic_triggers) > 1
                ISI          = diff(mic_triggers);
                trial_starts = [1, find(ISI > GAP_THRESH) + 1];
                keep_mask    = true(1, length(mic_triggers));
                keep_mask(trial_starts) = false;
                mic_triggers = mic_triggers(keep_mask);
            end

            %% REF unique matching
            ref_binary  = ref_norm > REF_THRESHOLD;
            ref_edges   = find(diff(ref_binary) == 1);

            %% REF INTEGRITY CHECK: reconstruct missing REF triggers as NaN placeholders.
            %% Each NaN means "stimulus was presented but REF trigger not recorded here."
            %% This ensures REF rank always maps correctly to CSV row number.
            n_ref_found = length(ref_edges);
            if n_ref_found > 0 && n_ref_found ~= TRIALS_PER_BLOCK
                % Step 1: detect gaps > 1.5× median ITI and insert NaN placeholders
                if n_ref_found > 1
                    median_iti   = median(diff(ref_edges));
                    ref_edges_full = ref_edges(1);
                    for ri = 2:n_ref_found
                        gap = ref_edges(ri) - ref_edges_full(end);
                        n_missing_mid = round(gap / median_iti) - 1;
                        for mi = 1:n_missing_mid
                            ref_edges_full(end+1) = NaN; %#ok
                        end
                        ref_edges_full(end+1) = ref_edges(ri); %#ok
                    end
                else
                    ref_edges_full = ref_edges;
                end
                % Step 2: prepend remaining NaNs at the beginning (missing first trigger)
                n_still_missing = TRIALS_PER_BLOCK - length(ref_edges_full);
                if n_still_missing > 0
                    ref_edges_full = [NaN(1, n_still_missing), ref_edges_full];
                end
                % Step 3: trim to TRIALS_PER_BLOCK if somehow over-filled
                if length(ref_edges_full) > TRIALS_PER_BLOCK
                    ref_edges_full = ref_edges_full(1:TRIALS_PER_BLOCK);
                end
                fprintf('    [REF] Block %d: %d/%d REF triggers found. Inserted %d NaN placeholder(s).\n', ...
                        block, n_ref_found, TRIALS_PER_BLOCK, TRIALS_PER_BLOCK - n_ref_found);
                ref_edges = ref_edges_full;
            end

            ref_search_before_samples = round(REF_SEARCH_BEFORE_MS * FS / 1000);
            ref_search_after_samples  = round(REF_SEARCH_AFTER_MS  * FS / 1000);

            mic_to_ref_candidates = cell(1, length(mic_triggers));
            mic_to_ref_distance   = cell(1, length(mic_triggers));

            for i = 1:length(mic_triggers)
                mic_sample  = mic_triggers(i);
                candidates  = ref_edges(ref_edges >= (mic_sample - ref_search_before_samples) & ...
                                        ref_edges <= (mic_sample + ref_search_after_samples));
                mic_to_ref_candidates{i} = candidates;
                mic_to_ref_distance{i}   = abs(candidates - mic_sample);
            end

            matched_ref = zeros(1, length(mic_triggers));
            used_refs   = false(1, length(ref_edges));

            for i = 1:length(mic_triggers)
                if isempty(mic_to_ref_candidates{i})
                    matched_ref(i) = NaN;
                    continue;
                end
                candidates       = mic_to_ref_candidates{i};
                distances        = mic_to_ref_distance{i};
                available_mask   = true(1, length(candidates));
                for j = 1:length(candidates)
                    ref_idx = find(ref_edges == candidates(j), 1);
                    if ~isempty(ref_idx) && used_refs(ref_idx)
                        available_mask(j) = false;
                    end
                end
                if any(available_mask)
                    av_cand = candidates(available_mask);
                    av_dist = distances(available_mask);
                    [~, ci] = min(av_dist);
                    best_ref = av_cand(ci);
                    matched_ref(i) = best_ref;
                    ref_idx = find(ref_edges == best_ref, 1);
                    used_refs(ref_idx) = true;
                else
                    matched_ref(i) = NaN;
                end
            end

            validated_mask       = ~isnan(matched_ref);
            triggers             = mic_triggers(validated_mask);
            matched_ref_valid    = matched_ref(validated_mask);  % REF sample for each trigger

            %% REF-RANK MAPPING:
            %% ref_edges may now contain NaN placeholders for missing triggers.
            %% For each validated trigger, find its matched REF sample in ref_edges
            %% (only among real/non-NaN entries). Its 1-based index = CSV row offset.
            real_ref_mask   = ~isnan(ref_edges);
            real_ref_samples = ref_edges(real_ref_mask);
            real_ref_ranks   = find(real_ref_mask);  % rank in full ref_edges = CSV offset

            trigger_csv_offsets = nan(1, length(triggers));
            for ti = 1:length(triggers)
                match_in_real = find(real_ref_samples == matched_ref_valid(ti), 1);
                if ~isempty(match_in_real)
                    trigger_csv_offsets(ti) = real_ref_ranks(match_in_real);
                end
            end

            %% block_start_raw = first CSV row of this block (fixed step = TRIALS_PER_BLOCK)
            block_start_raw   = raw_trial_counter + 1;
            n_mic_postISI     = length(mic_triggers);
            raw_trial_counter = raw_trial_counter + TRIALS_PER_BLOCK;

            %% SPEED-LABEL VERIFICATION via REF rank → CSV direct lookup.
            %% Each epoch's speed label comes directly from trial_info(block_start_raw + ref_rank - 1).
            %% Unmatched REF ranks = stimuli with no valid MIC detection (e.g. first trigger issue).
            if ~isempty(trial_info) && ~all(isnan(trigger_csv_offsets))
                valid_offsets = trigger_csv_offsets(~isnan(trigger_csv_offsets));
                mapped_rows   = block_start_raw + valid_offsets - 1;
                in_range      = mapped_rows >= 1 & mapped_rows <= height(trial_info);

                if any(in_range)
                    mapped_speeds  = trial_info.speed(mapped_rows(in_range));
                    matched_ranks  = valid_offsets(in_range);
                    all_ranks      = 1:TRIALS_PER_BLOCK;
                    unmatched_ranks = setdiff(all_ranks, matched_ranks);

                    %% Speed count string
                    u_spd = unique(mapped_speeds);
                    spd_str = '';
                    for si = 1:length(u_spd)
                        spd_str = [spd_str, sprintf('%.0fm/s:%d ', u_spd(si), sum(mapped_speeds==u_spd(si)))]; %#ok
                    end

                    %% Unmatched REF ranks and their CSV speeds (= missed stimuli)
                    unmatched_str = '';
                    if ~isempty(unmatched_ranks)
                        for ui = 1:min(5, length(unmatched_ranks))   % show max 5
                            ur  = unmatched_ranks(ui);
                            row = block_start_raw + ur - 1;
                            if row >= 1 && row <= height(trial_info)
                                unmatched_str = [unmatched_str, ...
                                    sprintf('rank%d(%.0fm/s) ', ur, trial_info.speed(row))]; %#ok
                            end
                        end
                        if length(unmatched_ranks) > 5
                            unmatched_str = [unmatched_str, sprintf('... +%d more', length(unmatched_ranks)-5)];
                        end
                    end

                    fprintf('    [SPD] Block %d: %d/%d REF matched → %s', ...
                            block, sum(in_range), TRIALS_PER_BLOCK, strtrim(spd_str));
                    if ~isempty(unmatched_ranks)
                        fprintf(' | unmatched: %s', strtrim(unmatched_str));
                    end
                end
            end

            %% MIC count info (echoes > 150 are harmless; < 150 = missed detections)
            if n_mic_postISI ~= TRIALS_PER_BLOCK
                fprintf('    [MIC] Block %d: %d mic triggers (expected %d, diff %+d)\n', ...
                        block, n_mic_postISI, TRIALS_PER_BLOCK, n_mic_postISI - TRIALS_PER_BLOCK);
            end

            if length(triggers) < 50
                continue;
            end

            %% Filtering
            C3_filt  = filtfilt(b, 1, C3  - mean(C3));
            C4_filt  = filtfilt(b, 1, C4  - mean(C4));
            CP3_filt = filtfilt(b, 1, CP3 - mean(CP3));
            CP4_filt = filtfilt(b, 1, CP4 - mean(CP4));
            FC3_filt = filtfilt(b, 1, FC3 - mean(FC3));
            FC4_filt = filtfilt(b, 1, FC4 - mean(FC4));
            Cz_filt  = filtfilt(b, 1, Cz  - mean(Cz));
            right_ear_filt = filtfilt(b, 1, right_ear - mean(right_ear));
            left_ear_filt  = filtfilt(b, 1, left_ear  - mean(left_ear));

            %% Re-referencing
            switch REFERENCE_METHOD
                case 'none'
                    C3_reref  = C3_filt;   C4_reref  = C4_filt;
                    CP3_reref = CP3_filt;  CP4_reref = CP4_filt;
                    FC3_reref = FC3_filt;  FC4_reref = FC4_filt;

                case 'earlobes'
                    earlobes_avg = (right_ear_filt + left_ear_filt) / 2;
                    C3_reref  = C3_filt  - earlobes_avg;
                    C4_reref  = C4_filt  - earlobes_avg;
                    CP3_reref = CP3_filt - earlobes_avg;
                    CP4_reref = CP4_filt - earlobes_avg;
                    FC3_reref = FC3_filt - earlobes_avg;
                    FC4_reref = FC4_filt - earlobes_avg;

                case 'left_ear'
                    C3_reref  = C3_filt  - left_ear_filt;
                    C4_reref  = C4_filt  - left_ear_filt;
                    CP3_reref = CP3_filt - left_ear_filt;
                    CP4_reref = CP4_filt - left_ear_filt;
                    FC3_reref = FC3_filt - left_ear_filt;
                    FC4_reref = FC4_filt - left_ear_filt;

                case 'right_ear'
                    C3_reref  = C3_filt  - right_ear_filt;
                    C4_reref  = C4_filt  - right_ear_filt;
                    CP3_reref = CP3_filt - right_ear_filt;
                    CP4_reref = CP4_filt - right_ear_filt;
                    FC3_reref = FC3_filt - right_ear_filt;
                    FC4_reref = FC4_filt - right_ear_filt;

                case 'Cz'
                    C3_reref  = C3_filt  - Cz_filt;
                    C4_reref  = C4_filt  - Cz_filt;
                    CP3_reref = CP3_filt - Cz_filt;
                    CP4_reref = CP4_filt - Cz_filt;
                    FC3_reref = FC3_filt - Cz_filt;
                    FC4_reref = FC4_filt - Cz_filt;

                case 'hemisphere'
                    ipsi_avg  = (C3_filt + CP3_filt + FC3_filt) / 3;
                    contra_avg = (C4_filt + CP4_filt + FC4_filt) / 3;
                    C4_reref  = C4_filt  - ipsi_avg;
                    CP4_reref = CP4_filt - ipsi_avg;
                    FC4_reref = FC4_filt - ipsi_avg;
                    C3_reref  = C3_filt  - contra_avg;
                    CP3_reref = CP3_filt - contra_avg;
                    FC3_reref = FC3_filt - contra_avg;

                case 'average'
                    common_avg = (C3_filt+C4_filt+CP3_filt+CP4_filt+FC3_filt+FC4_filt)/6;
                    C3_reref  = C3_filt  - common_avg;
                    C4_reref  = C4_filt  - common_avg;
                    CP3_reref = CP3_filt - common_avg;
                    CP4_reref = CP4_filt - common_avg;
                    FC3_reref = FC3_filt - common_avg;
                    FC4_reref = FC4_filt - common_avg;

                case 'laplacian'
                    C3_reref  = C3_filt  - (CP3_filt + FC3_filt) / 2;
                    C4_reref  = C4_filt  - (CP4_filt + FC4_filt) / 2;
                    CP3_reref = CP3_filt - (C3_filt  + FC3_filt) / 2;
                    CP4_reref = CP4_filt - (C4_filt  + FC4_filt) / 2;
                    FC3_reref = FC3_filt - (C3_filt  + CP3_filt) / 2;
                    FC4_reref = FC4_filt - (C4_filt  + CP4_filt) / 2;
            end

            %% Epoch extraction
            n_trials = length(triggers);
            block_epochs_C3  = zeros(n_trials, n_samples);
            block_epochs_C4  = zeros(n_trials, n_samples);
            block_epochs_CP3 = zeros(n_trials, n_samples);
            block_epochs_CP4 = zeros(n_trials, n_samples);
            block_epochs_FC3 = zeros(n_trials, n_samples);
            block_epochs_FC4 = zeros(n_trials, n_samples);

            valid_trials = [];
            for i = 1:n_trials
                trial_idx  = triggers(i) - MICROPHONE_DELAY_SAMPLES;
                epoch_idx  = trial_idx + epoch_samples;
                if min(epoch_idx) > 0 && max(epoch_idx) <= length(C3_reref)
                    block_epochs_C3(i,:) = C3_reref(epoch_idx);
                    block_epochs_C4(i,:) = C4_reref(epoch_idx);
                    block_epochs_CP3(i,:) = CP3_reref(epoch_idx);
                    block_epochs_CP4(i,:) = CP4_reref(epoch_idx);
                    block_epochs_FC3(i,:) = FC3_reref(epoch_idx);
                    block_epochs_FC4(i,:) = FC4_reref(epoch_idx);
                    valid_trials(end+1) = i;
                end
            end

            block_epochs_C3  = block_epochs_C3(valid_trials,:);
            block_epochs_C4  = block_epochs_C4(valid_trials,:);
            block_epochs_CP3 = block_epochs_CP3(valid_trials,:);
            block_epochs_CP4 = block_epochs_CP4(valid_trials,:);
            block_epochs_FC3 = block_epochs_FC3(valid_trials,:);
            block_epochs_FC4 = block_epochs_FC4(valid_trials,:);

            %% Baseline correction
            for i = 1:size(block_epochs_C3, 1)
                block_epochs_C3(i,:) = block_epochs_C3(i,:) - mean(block_epochs_C3 (i,baseline_idx));
                block_epochs_C4(i,:) = block_epochs_C4(i,:) - mean(block_epochs_C4 (i,baseline_idx));
                block_epochs_CP3(i,:) = block_epochs_CP3(i,:) - mean(block_epochs_CP3(i,baseline_idx));
                block_epochs_CP4(i,:) = block_epochs_CP4(i,:) - mean(block_epochs_CP4(i,baseline_idx));
                block_epochs_FC3(i,:) = block_epochs_FC3(i,:) - mean(block_epochs_FC3(i,baseline_idx));
                block_epochs_FC4(i,:) = block_epochs_FC4(i,:) - mean(block_epochs_FC4(i,baseline_idx));
            end

            %% Artifact rejection
            artifact_free = true(size(block_epochs_C3, 1), 1);
            for i = 1:size(block_epochs_C3, 1)
                if max(abs(block_epochs_C3(i,:))) > ARTIFACT_THRESHOLD || ...
                   max(abs(block_epochs_C4(i,:))) > ARTIFACT_THRESHOLD || ...
                   max(abs(block_epochs_CP3(i,:))) > ARTIFACT_THRESHOLD || ...
                   max(abs(block_epochs_CP4(i,:))) > ARTIFACT_THRESHOLD || ...
                   max(abs(block_epochs_FC3(i,:))) > ARTIFACT_THRESHOLD || ...
                   max(abs(block_epochs_FC4(i,:))) > ARTIFACT_THRESHOLD
                    artifact_free(i) = false;
                end
            end

 % Map trial_info using ORIGINAL trigger indices
            %% and keep epochs/trial_info in sync.
            %%
            %% valid_and_clean_local: trigger indices within block that are valid AND clean
            %% global_indices:        position in participant full trial_info table
            %% in_bounds:             which global_indices are within trial_info row count
            %%
            %% Only epochs WITH a matching trial_info row are kept.
            %% This guarantees size(GROUP_C3,1) == height(GROUP_trial_info).

            % Apply artifact rejection to epochs
            block_epochs_C3  = block_epochs_C3(artifact_free,:);
            block_epochs_C4  = block_epochs_C4(artifact_free,:);
            block_epochs_CP3 = block_epochs_CP3(artifact_free,:);
            block_epochs_CP4 = block_epochs_CP4(artifact_free,:);
            block_epochs_FC3 = block_epochs_FC3(artifact_free,:);
            block_epochs_FC4 = block_epochs_FC4(artifact_free,:);

            % Compute global trial_info indices for clean epochs
            valid_and_clean_local = valid_trials(artifact_free);
            %% Map each surviving epoch to its CSV row via REF rank.
            %% valid_and_clean_local = indices into the 'triggers' array.
            %% trigger_csv_offsets(i) = REF rank of trigger i = correct CSV offset.
            raw_offsets    = trigger_csv_offsets(valid_and_clean_local);
            has_offset     = ~isnan(raw_offsets);
            if any(~has_offset)
                fprintf('    [WARN] Block %d: %d epoch(s) have no REF rank – dropped.\n', ...
                        block, sum(~has_offset));
                % filter these out now
                block_epochs_C3  = block_epochs_C3(has_offset,:);
                block_epochs_C4  = block_epochs_C4(has_offset,:);
                block_epochs_CP3 = block_epochs_CP3(has_offset,:);
                block_epochs_CP4 = block_epochs_CP4(has_offset,:);
                block_epochs_FC3 = block_epochs_FC3(has_offset,:);
                block_epochs_FC4 = block_epochs_FC4(has_offset,:);
                raw_offsets      = raw_offsets(has_offset);
            end
            global_indices = block_start_raw + raw_offsets - 1;

            % Keep only epochs where trial_info row exists (safety check)
            %% With the corrected raw_trial_counter this should never drop anything.
            in_bounds  = global_indices <= height(trial_info);
            n_dropped  = sum(~in_bounds);
            if n_dropped > 0
                warning('Block %d: %d epoch(s) dropped – global_index exceeds trial_info rows (%d). Check trigger counting.', ...
                        block, n_dropped, height(trial_info));
            end
            block_epochs_C3  = block_epochs_C3(in_bounds,:);
            block_epochs_C4  = block_epochs_C4(in_bounds,:);
            block_epochs_CP3 = block_epochs_CP3(in_bounds,:);
            block_epochs_CP4 = block_epochs_CP4(in_bounds,:);
            block_epochs_FC3 = block_epochs_FC3(in_bounds,:);
            block_epochs_FC4 = block_epochs_FC4(in_bounds,:);

            % Collect matching trial_info rows
            valid_global = global_indices(in_bounds);
            for i = 1:length(valid_global)
                all_trial_info_clean = [all_trial_info_clean; trial_info(valid_global(i),:)];
            end

            %% Per-block diagnostic summary
            n_trig_blk   = length(triggers);          % validated triggers in block
            n_valid_blk  = length(valid_trials);       % survived epoch boundary check
            n_clean_blk  = sum(artifact_free);         % survived artifact rejection
            n_inbounds_blk = sum(in_bounds);           % survived trial_info bounds check
            fprintf('    Block %d: %3d triggers → %3d valid epochs → %3d clean → %3d in CSV → %3d kept\n', ...
                    block, n_trig_blk, n_valid_blk, n_clean_blk, n_inbounds_blk, n_inbounds_blk);

            %% Concatenate
            all_epochs_C3  = [all_epochs_C3;  block_epochs_C3 ];
            all_epochs_C4  = [all_epochs_C4;  block_epochs_C4 ];
            all_epochs_CP3 = [all_epochs_CP3; block_epochs_CP3];
            all_epochs_CP4 = [all_epochs_CP4; block_epochs_CP4];
            all_epochs_FC3 = [all_epochs_FC3; block_epochs_FC3];
            all_epochs_FC4 = [all_epochs_FC4; block_epochs_FC4];
        end  % end block loop

        n_epochs = size(all_epochs_C3, 1);
        n_info   = height(all_trial_info_clean);
        if n_epochs ~= n_info
            error('SYNC ERROR in %s: %d epochs but %d trial_info rows. Check trigger counting.', ...
                  PARTICIPANT_ID, n_epochs, n_info);
        end
        fprintf('  Clean epochs: %d  (trial_info rows: %d ✓)\n', n_epochs, n_info);

        if n_epochs == 0
            warning('No valid epochs for %s', PARTICIPANT_ID);
            failed_participants{end+1} = PARTICIPANT_ID;
            continue;
        end

        %% Add to group
        GROUP_C3  = [GROUP_C3;  all_epochs_C3 ];
        GROUP_C4  = [GROUP_C4;  all_epochs_C4 ];
        GROUP_CP3 = [GROUP_CP3; all_epochs_CP3];
        GROUP_CP4 = [GROUP_CP4; all_epochs_CP4];
        GROUP_FC3 = [GROUP_FC3; all_epochs_FC3];
        GROUP_FC4 = [GROUP_FC4; all_epochs_FC4];
        GROUP_trial_info         = [GROUP_trial_info; all_trial_info_clean];
        new_pid_labels = repmat({PARTICIPANT_ID}, n_epochs, 1);
        GROUP_participant_labels = [GROUP_participant_labels; new_pid_labels];

        successful_participants{end+1} = PARTICIPANT_ID;
        fprintf('  ✓ Added to group\n');

    catch ME
        fprintf('  ❌ ERROR: %s\n', ME.message);
        failed_participants{end+1} = PARTICIPANT_ID;
    end
end  % end participant loop

fprintf('  Successful: %d participants\n', length(successful_participants));
fprintf('  Total epochs: %d\n', size(GROUP_C3, 1));

unique_speeds = unique(GROUP_trial_info.speed);
n_speeds      = length(unique_speeds);
fprintf('  Speeds: %d levels (', n_speeds);
fprintf('%.1f  ', unique_speeds);
fprintf('m/s)\n\n');

%% Save full dataset (no PARTICIPANT_AVERAGES yet – computed after balancing)
group_file = sprintf('GROUP_SEP_Data_%s.mat', REFERENCE_METHOD);
save(group_file, ...
     'GROUP_C3','GROUP_C4','GROUP_CP3','GROUP_CP4','GROUP_FC3','GROUP_FC4', ...
     'GROUP_trial_info','GROUP_participant_labels', ...
     'successful_participants','failed_participants', ...
     'time_vector','time_ms','P1_idx','N1_idx','P2_idx', ...
     'unique_speeds','REFERENCE_METHOD','FS', '-v7.3');
fprintf('  Saved: %s\n\n', group_file);

%% PART 2: BALANCED EPOCH SELECTION

fprintf('Target: %d epochs per speed per participant\n\n', N_EPOCHS_PER_SPEED);

%% Epoch counts before selection
fprintf('Epoch counts BEFORE balanced selection:\n\n');
fprintf('Participant | ');
for s = 1:n_speeds; fprintf('%6.1f m/s | ', unique_speeds(s)); end
fprintf('Total\n');
fprintf('------------|');
for s = 1:n_speeds; fprintf('----------|'); end
fprintf('------\n');

%% Check which participants have enough epochs per speed
eligible_participants = {};
excluded_participants = {};

for p = 1:length(successful_participants)
    pid          = successful_participants{p};
    p_mask       = strcmp(GROUP_participant_labels, pid);
    p_trial_info = GROUP_trial_info(p_mask, :);
    fprintf('%11s | ', pid);
    total     = 0;
    min_count = Inf;
    for s = 1:n_speeds
        cnt = sum([p_trial_info.speed] == unique_speeds(s));
        fprintf('%9d | ', cnt);
        total     = total + cnt;
        min_count = min(min_count, cnt);
    end
    if min_count >= N_EPOCHS_PER_SPEED
        fprintf('%5d  ✓\n', total);
        eligible_participants{end+1} = pid;
    else
        fprintf('%5d  ✗ EXCLUDED (min=%d < %d)\n', total, min_count, N_EPOCHS_PER_SPEED);
        excluded_participants{end+1} = pid;
    end
end

fprintf('\n  Eligible:  %d participants\n', length(eligible_participants));
if ~isempty(excluded_participants)
    fprintf('  Excluded:  ');
    for i = 1:length(excluded_participants)
        fprintf('%s ', excluded_participants{i});
    end
    fprintf('(< %d epochs in at least one speed condition)\n', N_EPOCHS_PER_SPEED);
end

if isempty(eligible_participants)
    error('No participants have %d epochs per speed! Lower N_EPOCHS_PER_SPEED.', N_EPOCHS_PER_SPEED);
end

%% Select balanced epochs (eligible participants only)
BALANCED_C3  = [];  BALANCED_C4  = [];
BALANCED_CP3 = [];  BALANCED_CP4 = [];
BALANCED_FC3 = [];  BALANCED_FC4 = [];
BALANCED_trial_info         = [];
BALANCED_participant_labels = [];

for p = 1:length(eligible_participants)
    pid              = eligible_participants{p};
    p_mask           = strcmp(GROUP_participant_labels, pid);
    p_indices        = find(p_mask);
    p_trial_info     = GROUP_trial_info(p_mask, :);

    for s = 1:n_speeds
        speed         = unique_speeds(s);
        s_local       = find([p_trial_info.speed] == speed);
        sel_local     = s_local(1:N_EPOCHS_PER_SPEED);
        sel_global    = p_indices(sel_local);

        BALANCED_C3  = [BALANCED_C3;  GROUP_C3(sel_global,:)];
        BALANCED_C4  = [BALANCED_C4;  GROUP_C4(sel_global,:)];
        BALANCED_CP3 = [BALANCED_CP3; GROUP_CP3(sel_global,:)];
        BALANCED_CP4 = [BALANCED_CP4; GROUP_CP4(sel_global,:)];
        BALANCED_FC3 = [BALANCED_FC3; GROUP_FC3(sel_global,:)];
        BALANCED_FC4 = [BALANCED_FC4; GROUP_FC4(sel_global,:)];
        new_info = GROUP_trial_info(sel_global,:);
        BALANCED_trial_info = [BALANCED_trial_info; new_info];
        new_labels = GROUP_participant_labels(sel_global);
        BALANCED_participant_labels = [BALANCED_participant_labels; new_labels];
    end
    fprintf('  ✓ %s: selected %d epochs (%d per speed)\n', ...
            pid, N_EPOCHS_PER_SPEED*n_speeds, N_EPOCHS_PER_SPEED);
end

%% Override successful_participants for downstream analysis
successful_participants = eligible_participants;
n_participants = length(successful_participants);
expected_total = N_EPOCHS_PER_SPEED * n_speeds * n_participants;
fprintf('\n  Total balanced: %d  (expected %d)  %s\n\n', ...
        size(BALANCED_C3,1), expected_total, ...
        iif(size(BALANCED_C3,1)==expected_total, '✓', '✗'));

 % Compute PARTICIPANT_AVERAGES from BALANCED data
fprintf('Computing participant averages from balanced data...\n');
PARTICIPANT_AVERAGES = struct();
ELECTRODE_DATA = struct();
ELECTRODE_DATA.C3  = BALANCED_C3;
ELECTRODE_DATA.C4  = BALANCED_C4;
ELECTRODE_DATA.CP3 = BALANCED_CP3;
ELECTRODE_DATA.CP4 = BALANCED_CP4;
ELECTRODE_DATA.FC3 = BALANCED_FC3;
ELECTRODE_DATA.FC4 = BALANCED_FC4;
for p = 1:n_participants
    pid    = successful_participants{p};
    p_mask = strcmp(BALANCED_participant_labels, pid);
    for e = 1:N_ELECTRODES
        el = ELECTRODES{e};
        PARTICIPANT_AVERAGES.(pid).(el) = mean(ELECTRODE_DATA.(el)(p_mask,:), 1);
    end
end

%% Save balanced dataset
balanced_file = sprintf('BALANCED_SEP_Data_%deps_%s.mat', N_EPOCHS_PER_SPEED, REFERENCE_METHOD);
save(balanced_file, ...
     'BALANCED_C3','BALANCED_C4','BALANCED_CP3','BALANCED_CP4', ...
     'BALANCED_FC3','BALANCED_FC4', ...
     'BALANCED_trial_info','BALANCED_participant_labels', ...
     'PARTICIPANT_AVERAGES','successful_participants', ...
     'time_vector','time_ms','P1_idx','N1_idx','P2_idx','P1N1_idx','N1P2_idx', ...
     'FS','unique_speeds','REFERENCE_METHOD','N_EPOCHS_PER_SPEED','-v7.3');
fprintf('  Saved: %s\n\n', balanced_file);

%% PART 3: VISUALIZATIONS

participant_colors = lines(n_participants);

%% Pre-compute speed-stratified grand averages per electrode
 % SEM computed over participant means (correct N=n_participants)
GA  = struct();  % Grand average:  n_speeds × n_samples
SEM = struct();  % SEM:            n_speeds × n_samples

for e = 1:N_ELECTRODES
    el        = ELECTRODES{e};
    GA.(el)   = zeros(n_speeds, n_samples);
    SEM.(el)  = zeros(n_speeds, n_samples);

    for s = 1:n_speeds
        speed  = unique_speeds(s);
        vp_avg = zeros(n_participants, n_samples);

        for pp = 1:n_participants
            pid      = successful_participants{pp};
            pp_mask  = strcmp(BALANCED_participant_labels, pid);
            sp_mask  = [BALANCED_trial_info.speed] == speed;
            combined = pp_mask & sp_mask;
            vp_avg(pp,:) = mean(ELECTRODE_DATA.(el)(combined,:), 1);
        end

        GA.(el)(s,:)  = mean(vp_avg, 1);
        SEM.(el)(s,:) = std(vp_avg, 0, 1) / sqrt(n_participants);
    end
end

%% Figure 1: Grand Average – all speeds per electrode
fprintf('Creating Figure 1: Grand Average...\n');
fig1 = figure('Position',[50,50,1800,1200],'Color','white');
for e = 1:N_ELECTRODES
    el = ELECTRODES{e};
    subplot(2, 3, e); hold on;
    for s = 1:n_speeds
        fill([time_ms, fliplr(time_ms)], ...
             [GA.(el)(s,:)+SEM.(el)(s,:), fliplr(GA.(el)(s,:)-SEM.(el)(s,:))], ...
             SPEED_COLORS(s,:), 'FaceAlpha', 0.2, 'EdgeColor', 'none');
    end
    for s = 1:n_speeds
        plot(time_ms, GA.(el)(s,:), 'Color', SPEED_COLORS(s,:), 'LineWidth', 2.5, ...
             'DisplayName', sprintf('%.1f m/s', unique_speeds(s)));
    end
    xline(0,'--k','LineWidth',1.5,'Alpha',0.7);
    yline(0,'-k','LineWidth',0.5,'Alpha',0.5);
    yl = ylim;
    patch([P1_START P1_END P1_END P1_START]*1000, [yl(1) yl(1) yl(2) yl(2)], ...
          [0.9 0.9 0.9], 'EdgeColor','none','FaceAlpha',0.3);
    xlim([time_ms(1) time_ms(end)]);
    xlabel('Time (ms)','FontSize',12,'FontWeight','bold');
    ylabel('Amplitude (µV)','FontSize',12,'FontWeight','bold');
    title(el,'FontSize',14,'FontWeight','bold');
    if e==1; legend('Location','best','FontSize',9); end
    grid on; set(gca,'FontSize',11); hold off;
end
sgtitle(sprintf('Grand Average SEP (Mean ± SEM, N=%d) | Ref: %s', n_participants, REFERENCE_METHOD), ...
        'FontSize',16,'FontWeight','bold');
fname = sprintf('BALANCED_GrandAverage_%deps_%s.png', N_EPOCHS_PER_SPEED, REFERENCE_METHOD);
saveas(fig1, fname);  fprintf('  Saved: %s\n', fname);

%% Figure 2: Individual overlay 
fprintf('Creating Figure 2: Individual overlay...\n');
fig2 = figure('Position',[100,100,1800,1200],'Color','white');
for e = 1:N_ELECTRODES
    el = ELECTRODES{e};
    subplot(2, 3, e); hold on;
    for pp = 1:n_participants
        pid = successful_participants{pp};
        plot(time_ms, PARTICIPANT_AVERAGES.(pid).(el), '-', ...
             'LineWidth',1, 'Color',[participant_colors(pp,:) 0.5], ...
             'DisplayName', pid);
    end
    plot(time_ms, mean(ELECTRODE_DATA.(el),1), 'k-', ...
         'LineWidth',4, 'DisplayName','Grand Mean');
    xline(0,'--k','LineWidth',1.5,'Alpha',0.7);
    yline(0,'-k','LineWidth',0.5,'Alpha',0.5);
    xlim([time_ms(1) time_ms(end)]);
    xlabel('Time (ms)','FontSize',12,'FontWeight','bold');
    ylabel('Amplitude (µV)','FontSize',12,'FontWeight','bold');
    title(el,'FontSize',14,'FontWeight','bold');
    if e==1; legend('Location','best','FontSize',8); end
    grid on; set(gca,'FontSize',11); hold off;
end
sgtitle(sprintf('Individual Participants + Grand Average (N=%d) | Ref: %s', ...
        n_participants, REFERENCE_METHOD), 'FontSize',16,'FontWeight','bold');
fname = sprintf('BALANCED_IndividualOverlay_%deps_%s.png', N_EPOCHS_PER_SPEED, REFERENCE_METHOD);
saveas(fig2, fname);  fprintf('  Saved: %s\n', fname);

%% Figure 3: Speed comparison grid
fprintf('Creating Figure 3: Speed comparison grid...\n');
fig3 = figure('Position',[200,200,2000,1200],'Color','white');
for e = 1:N_ELECTRODES
    el = ELECTRODES{e};
    for s = 1:n_speeds
        subplot(N_ELECTRODES, n_speeds, (e-1)*n_speeds + s);
        plot(time_ms, GA.(el)(s,:), '-', 'LineWidth',2.5, 'Color',SPEED_COLORS(s,:));
        hold on;
        fill([time_ms, fliplr(time_ms)], ...
             [GA.(el)(s,:)+SEM.(el)(s,:), fliplr(GA.(el)(s,:)-SEM.(el)(s,:))], ...
             SPEED_COLORS(s,:),'FaceAlpha',0.2,'EdgeColor','none');
        xline(0,'--k','LineWidth',1); yline(0,'-k','LineWidth',0.5);
        hold off;
        xlim([time_ms(1) time_ms(end)]);
        if s==1; ylabel(el,'FontSize',11,'FontWeight','bold'); end
        if e==1; title(sprintf('%.1f m/s',unique_speeds(s)),'FontSize',12,'FontWeight','bold'); end
        if e==N_ELECTRODES; xlabel('Time (ms)','FontSize',10); end
        if s>1; set(gca,'YTickLabel',[]); end
        grid on; set(gca,'FontSize',9);
    end
end
sgtitle(sprintf('Speed Comparison Grid (N=%d) | Ref: %s', n_participants, REFERENCE_METHOD), ...
        'FontSize',16,'FontWeight','bold');
fname = sprintf('BALANCED_SpeedGrid_%deps_%s.png', N_EPOCHS_PER_SPEED, REFERENCE_METHOD);
saveas(fig3, fname);  fprintf('  Saved: %s\n', fname);

%% Figure 4: Right Hemisphere (contralateral) – horizontal
fprintf('Creating Figure 4: Right hemisphere...\n');
fig4 = figure('Position',[50,100,1800,500],'Color','white');
right_els = {'FC4','C4','CP4'};
for e = 1:3
    el = right_els{e};
    subplot(3,1,e); hold on;
    for s = 1:n_speeds
        fill([time_ms, fliplr(time_ms)], ...
             [GA.(el)(s,:)+SEM.(el)(s,:), fliplr(GA.(el)(s,:)-SEM.(el)(s,:))], ...
             SPEED_COLORS(s,:),'FaceAlpha',0.2,'EdgeColor','none');
    end
    for s = 1:n_speeds
        plot(time_ms, GA.(el)(s,:), 'Color',SPEED_COLORS(s,:),'LineWidth',2.5, ...
             'DisplayName',sprintf('%.1f m/s',unique_speeds(s)));
    end
    xline(0,'--k','LineWidth',1.5,'Alpha',0.7);
    yline(0,'-k','LineWidth',0.5,'Alpha',0.5);
    xlim([time_ms(1) time_ms(end)]);
    xlabel('Time (ms)','FontSize',12,'FontWeight','bold');
    ylabel('Amplitude (µV)','FontSize',12,'FontWeight','bold');
    title(el,'FontSize',14,'FontWeight','bold');
    if e==1; legend('Location','best','FontSize',9); end
    grid on; set(gca,'FontSize',11); hold off;
end
sgtitle(sprintf('Grand Average – Right Hemisphere (Contralateral) | Ref: %s', REFERENCE_METHOD), ...
        'FontSize',16,'FontWeight','bold');
fname = sprintf('GrandAverage_RightHemisphere_%s.png', REFERENCE_METHOD);
saveas(fig4, fname);  fprintf('  Saved: %s\n', fname);

%% Figure 5: Left Hemisphere (ipsilateral) – horizontal
fprintf('Creating Figure 5: Left hemisphere...\n');
fig5 = figure('Position',[50,650,1800,500],'Color','white');
left_els = {'FC3','C3','CP3'};
for e = 1:3
    el = left_els{e};
    subplot(3,1,e); hold on;
    for s = 1:n_speeds
        fill([time_ms, fliplr(time_ms)], ...
             [GA.(el)(s,:)+SEM.(el)(s,:), fliplr(GA.(el)(s,:)-SEM.(el)(s,:))], ...
             SPEED_COLORS(s,:),'FaceAlpha',0.2,'EdgeColor','none');
    end
    for s = 1:n_speeds
        plot(time_ms, GA.(el)(s,:), 'Color',SPEED_COLORS(s,:),'LineWidth',2.5, ...
             'DisplayName',sprintf('%.1f m/s',unique_speeds(s)));
    end
    xline(0,'--k','LineWidth',1.5,'Alpha',0.7);
    yline(0,'-k','LineWidth',0.5,'Alpha',0.5);
    xlim([time_ms(1) time_ms(end)]);
    xlabel('Time (ms)','FontSize',12,'FontWeight','bold');
    ylabel('Amplitude (µV)','FontSize',12,'FontWeight','bold');
    title(el,'FontSize',14,'FontWeight','bold');
    if e==1; legend('Location','best','FontSize',9); end
    grid on; set(gca,'FontSize',11); hold off;
end
sgtitle(sprintf('Grand Average – Left Hemisphere (Ipsilateral) | Ref: %s', REFERENCE_METHOD), ...
        'FontSize',16,'FontWeight','bold');
fname = sprintf('GrandAverage_LeftHemisphere_%s.png', REFERENCE_METHOD);
saveas(fig5, fname);  fprintf('  Saved: %s\n', fname);

%% Figure 6: Electrode heatmaps (trials × time per speed)
fprintf('Creating Figure 6: Single-trial heatmaps...\n');
for e = 1:N_ELECTRODES
    el         = ELECTRODES{e};
    el_data    = ELECTRODE_DATA.(el);
    fig6       = figure('Position',[50+(e-1)*80,50,2200,800],'Color','white');
    ax_hm      = gobjects(1, n_speeds);

    %% Panel geometry (leave 15% at top for sgtitle, 8% at bottom for xlabel)
    pw   = 0.1;   % war automatisch ~0.17, jetzt fester kleinerer Wert
    hgap = 0.05;   % Abstand zwischen Panels (war 0.02) 
    ph = 0.8;    % panel height – shorter than default to give sgtitle room
    yb = 0.09;    % bottom edge of panels

    for s = 1:n_speeds
        sp_mask = [BALANCED_trial_info.speed] == unique_speeds(s);
        sp_data = el_data(sp_mask,:);
        total_content = n_speeds*pw + (n_speeds-1)*hgap + 0.008 + 0.012;  % panels + colorbar
        x_start = (1 - total_content) / 2;  % zentriert
        xl = x_start + (s-1)*(pw+hgap);        %xl = 0.06 + (s-1)*(pw+0.02);   % left edge of this panel
        ax_hm(s) = axes('Position', [xl, yb, pw, ph]); %#ok<LAXES>
        imagesc(time_ms, 1:size(sp_data,1), sp_data);
        xlim([-200 400]);
        % caxis([-10 10]): ±10 µV covers the typical single-trial SEP range
        colormap(jet);  caxis([-10 10]);
        xlabel('Time (ms)','FontSize',11,'FontWeight','bold');
        if s==1; ylabel('Trial','FontSize',11,'FontWeight','bold'); end
        title(sprintf('%.1f m/s (n=%d)',unique_speeds(s),size(sp_data,1)), ...
              'FontSize',12,'FontWeight','bold');
        set(gca,'FontSize',10);
    end

    %% Single shared colorbar to the right of all panels
    p_last = get(ax_hm(n_speeds),'Position');
    cb = colorbar('Position',[p_last(1)+p_last(3)+0.008, yb, 0.012, ph]);
    ylabel(cb,'Amplitude (µV)','FontSize',11);

    %% sgtitle sits well above the subplot titles now
    sgtitle(sprintf('%s – Trials × Time | Ref: %s', el, REFERENCE_METHOD), ...
            'FontSize',14,'FontWeight','bold');
    fname = sprintf('BALANCED_Heatmap_%s_%deps_%s.png', el, N_EPOCHS_PER_SPEED, REFERENCE_METHOD);
    saveas(fig6, fname);  fprintf('  Saved: %s\n', fname);
end

%% PART 4: COMPREHENSIVE PEAK ANALYSIS (P1, N1, P2, COMPLEXES)

PEAK_RESULTS = struct();

for e = 1:N_ELECTRODES
    el      = ELECTRODES{e};
    el_data = ELECTRODE_DATA.(el);

    fprintf('Processing: %s\n', el);
    for s = 1:n_speeds
        speed     = unique_speeds(s);
        sp_label  = sprintf('speed_%d', round(speed*10));
        sp_mask   = [BALANCED_trial_info.speed] == speed;
        avg_wf    = mean(el_data(sp_mask,:), 1);

        [p1_amp, i] = max(avg_wf(P1_idx));  p1_lat = time_ms(P1_idx); p1_lat = p1_lat(i);
        [n1_amp, i] = min(avg_wf(N1_idx));  n1_lat = time_ms(N1_idx); n1_lat = n1_lat(i);
        seg_p2 = avg_wf(P2_idx);  t_p2 = time_ms(P2_idx);
        %% Find most prominent local maximum (version-compatible, no MinProminence needed)
        is_peak_p2 = false(size(seg_p2));
        for pk_i = 2:length(seg_p2)-1
            is_peak_p2(pk_i) = seg_p2(pk_i) > seg_p2(pk_i-1) && seg_p2(pk_i) > seg_p2(pk_i+1);
        end
        locs_p2 = find(is_peak_p2);
        pks_p2  = seg_p2(locs_p2);
        prom_p2 = zeros(size(pks_p2));
        for pk_i = 1:length(locs_p2)
            left_min  = min(seg_p2(1:locs_p2(pk_i)));
            right_min = min(seg_p2(locs_p2(pk_i):end));
            prom_p2(pk_i) = pks_p2(pk_i) - max(left_min, right_min);
        end
        keep_p2 = prom_p2 >= 0.05;
        locs_p2 = locs_p2(keep_p2);  pks_p2 = pks_p2(keep_p2);  prom_p2 = prom_p2(keep_p2);
        if ~isempty(pks_p2)
            [~, best_p2] = max(prom_p2);
            p2_amp = pks_p2(best_p2);
            p2_lat = t_p2(locs_p2(best_p2));
        else
            [p2_amp, pk_i] = max(seg_p2);
            p2_lat = t_p2(pk_i);
        end

        PEAK_RESULTS.(el).(sp_label).P1_amplitude  = p1_amp;
        PEAK_RESULTS.(el).(sp_label).P1_latency    = p1_lat;
        PEAK_RESULTS.(el).(sp_label).N1_amplitude  = n1_amp;
        PEAK_RESULTS.(el).(sp_label).N1_latency    = n1_lat;
        PEAK_RESULTS.(el).(sp_label).P2_amplitude  = p2_amp;
        PEAK_RESULTS.(el).(sp_label).P2_latency    = p2_lat;
        PEAK_RESULTS.(el).(sp_label).P1N1_amplitude = p1_amp - n1_amp;
        PEAK_RESULTS.(el).(sp_label).N1P2_amplitude = p2_amp - n1_amp;
        PEAK_RESULTS.(el).(sp_label).n_trials       = sum(sp_mask);

        fprintf('  %.1f m/s: P1=%.2f@%.0fms  N1=%.2f@%.0fms  P2=%.2f@%.0fms\n', ...
                speed, p1_amp, p1_lat, n1_amp, n1_lat, p2_amp, p2_lat);
    end
end

PEAK_RESULTS.metadata.electrodes    = ELECTRODES;
PEAK_RESULTS.metadata.unique_speeds = unique_speeds;
PEAK_RESULTS.metadata.P1_window     = [P1_START P1_END]*1000;
PEAK_RESULTS.metadata.N1_window     = [N1_START N1_END]*1000;
PEAK_RESULTS.metadata.P2_window     = [P2_START P2_END]*1000;

peak_file = sprintf('PEAK_ANALYSIS_Results_%s.mat', REFERENCE_METHOD);
save(peak_file, 'PEAK_RESULTS', '-v7.3');
fprintf('✓ Saved: %s\n\n', peak_file);

%% C4 summary table
fprintf('C4 Summary:\n');
for s = 1:n_speeds
    sl = sprintf('speed_%d', round(unique_speeds(s)*10));
            unique_speeds(s), ...
            PEAK_RESULTS.C4.(sl).P1_amplitude, PEAK_RESULTS.C4.(sl).N1_amplitude, ...
            PEAK_RESULTS.C4.(sl).P2_amplitude, PEAK_RESULTS.C4.(sl).P1N1_amplitude, ...
            PEAK_RESULTS.C4.(sl).N1P2_amplitude);
end

%% PART 5: WPSS ANALYSIS

WPSS_FS             = FS;
WPSS_FREQ_LIMITS    = [1 30];
WPSS_VOICES_PER_OCT = 20;
WPSS_GAMMA          = 3;
WPSS_BETA           = 3 * 3.29;

fprintf('Creating CWT filter bank...\n');
cwtObj = cwtfilterbank('SignalLength',     n_samples, ...
                       'SamplingFrequency', WPSS_FS, ...
                       'FrequencyLimits',   WPSS_FREQ_LIMITS, ...
                       'Wavelet',           'morse', ...
                       'VoicesPerOctave',   WPSS_VOICES_PER_OCT, ...
                       'WaveletParameters', [WPSS_GAMMA, WPSS_BETA]);

center_freqs = centerFrequencies(cwtObj);
n_freqs      = length(center_freqs);
fprintf('  Frequency bins: %d (%.2f – %.2f Hz)\n\n', ...
        n_freqs, min(center_freqs), max(center_freqs));

WPSS_results = struct();

for e = 1:N_ELECTRODES
    el      = ELECTRODES{e};
    el_data = ELECTRODE_DATA.(el);
    fprintf('Processing WPSS: %s\n', el);

    for s = 1:n_speeds
        speed    = unique_speeds(s);
        sp_label = sprintf('speed_%d', round(speed*10));
        sp_mask  = [BALANCED_trial_info.speed] == speed;
        data     = el_data(sp_mask,:)';      % n_samples × n_trials
        nSweeps  = size(data, 2);
        data     = data - mean(data, 1);     % remove DC per trial

        coeffs = complex(zeros(n_freqs, n_samples, nSweeps));
        for n = 1:nSweeps
            coeffs(:,:,n) = wt(cwtObj, data(:,n));
        end

        phi = angle(coeffs);
        WPSS_results.(el).(sp_label).wpss     = abs(mean(exp(1i.*phi), 3));
        WPSS_results.(el).(sp_label).powerST  = mean(abs(coeffs).^2, 3);
        WPSS_results.(el).(sp_label).powerAVG = abs(mean(coeffs, 3)).^2;
        WPSS_results.(el).(sp_label).n_trials = nSweeps;
        fprintf('  %.1f m/s: %d trials\n', speed, nSweeps);
    end
end

WPSS_results.metadata.center_frequencies = center_freqs;
WPSS_results.metadata.time_ms            = time_ms;
WPSS_results.metadata.unique_speeds      = unique_speeds;
WPSS_results.metadata.electrodes         = ELECTRODES;
WPSS_results.metadata.sampling_rate      = WPSS_FS;
WPSS_results.metadata.frequency_limits   = WPSS_FREQ_LIMITS;

wpss_file = sprintf('WPSS_Results_%s.mat', REFERENCE_METHOD);
save(wpss_file, 'WPSS_results', '-v7.3');
fprintf('✓ Saved: %s\n\n', wpss_file);

%% WPSS component summary (C4)
%  Frequency band: 5–15 Hz (theta–alpha).
%  Rationale: The WPSS spectrograms show that stimulus-evoked phase
%  synchronisation is concentrated in the ~5–20 Hz range, with the
%  strongest activity between 5 and 15 Hz.  This band captures the
%  dominant theta–alpha synchronisation without including higher
%  frequencies (>15 Hz) where WPSS values are near baseline.
%  Using a band average rather than a single frequency is more robust
%  against bin-level fluctuations and avoids the ambiguity of selecting
%  a single "dominant" frequency.
WPSS_BAND_LOW  = 5;
WPSS_BAND_HIGH = 15;

p1_tmask = time_ms >= 45  & time_ms <= 90;
n1_tmask = time_ms >= 100 & time_ms <= 160;
p2_tmask = time_ms >= 200 & time_ms <= 225;
comp_fmask = center_freqs >= WPSS_BAND_LOW & center_freqs <= WPSS_BAND_HIGH;

fprintf('WPSS analysis band: %d–%d Hz  (%d frequency bins)\n', ...
        WPSS_BAND_LOW, WPSS_BAND_HIGH, sum(comp_fmask));

wpss_p1 = zeros(n_speeds,1);
wpss_n1 = zeros(n_speeds,1);
wpss_p2 = zeros(n_speeds,1);

for s = 1:n_speeds
    sl = sprintf('speed_%d', round(unique_speeds(s)*10));
    w  = WPSS_results.C4.(sl).wpss;
    wpss_p1(s) = mean(mean(w(comp_fmask, p1_tmask)));
    wpss_n1(s) = mean(mean(w(comp_fmask, n1_tmask)));
    wpss_p2(s) = mean(mean(w(comp_fmask, p2_tmask)));
end

fprintf('WPSS Summary (C4, %d–%d Hz band average):\n', WPSS_BAND_LOW, WPSS_BAND_HIGH);
for s = 1:n_speeds
            unique_speeds(s), wpss_p1(s), wpss_n1(s), wpss_p2(s));
end
        mean(wpss_p1), mean(wpss_n1), mean(wpss_p2));

%% WPSS Grand Average 1D Plot – Band-averaged (5–15 Hz)
%  Instead of selecting a single "dominant" frequency, the WPSS time
%  course is computed as the mean across the theta–alpha band (5–15 Hz).
%  This approach:
%    (a) is robust against single-bin fluctuations,
%    (b) reflects the broadband synchronisation visible in the WPSS
%        spectrograms (Fig. WPSS),
%    (c) is consistent with the band used for the component table above.
%
%  Note: Lehser et al. (2018) used a single scale parameter (a=20,
%  ~15.36 Hz) because AM stimulation produces a narrowband response.
%  STM stimulation produces broadband synchronisation (see spectrograms),
%  making a band average the more appropriate summary measure.

fprintf('\nWPSS 1D plot: band-averaged over %d–%d Hz\n', ...
        WPSS_BAND_LOW, WPSS_BAND_HIGH);

% Also report peak frequency per speed for completeness
fprintf('  Peak frequencies per speed (within %d–%d Hz, full SEP window):\n', ...
        WPSS_BAND_LOW, WPSS_BAND_HIGH);
full_sep_tmask = time_ms >= 40 & time_ms <= 235;
for s = 1:n_speeds
    sl = sprintf('speed_%d', round(unique_speeds(s)*10));
    w  = WPSS_results.C4.(sl).wpss;
    freq_profile = mean(w(comp_fmask, full_sep_tmask), 2);
    band_freqs   = center_freqs(comp_fmask);
    [peak_val, peak_idx] = max(freq_profile);
    fprintf('    %.1f m/s: peak at %.2f Hz (WPSS = %.3f)\n', ...
            unique_speeds(s), band_freqs(peak_idx), peak_val);
end

fig_wpss1d = figure('Position', [100, 100, 900, 450], 'Color', 'white');
hold on;

for s = 1:n_speeds
    sl    = sprintf('speed_%d', round(unique_speeds(s)*10));
    w     = WPSS_results.C4.(sl).wpss;

    % Band-averaged WPSS: mean across all frequency bins in 5–15 Hz
    wpss_trace = mean(w(comp_fmask, :), 1);   % 1 × n_samples

    % Smooth slightly for readability (moving average, 5 samples ≈ 4ms)
    wpss_smooth = movmean(wpss_trace, 5);

    plot(time_ms, wpss_smooth, ...
         'Color',     SPEED_COLORS(s,:), ...
         'LineWidth', 2.0, ...
         'DisplayName', sprintf('%.1f m/s', unique_speeds(s)));
end

% Stimulus onset line
xline(0, '--k', 'LineWidth', 1.5, 'HandleVisibility', 'off');

% Formatting
xlim([-200 600]);
ylim([0 0.35]);
xlabel('Time (ms)', 'FontSize', 13, 'FontWeight', 'bold');
ylabel('WPSS', 'FontSize', 13, 'FontWeight', 'bold');
title(sprintf('Grand Average WPSS (%d–%d Hz) – C4 | Ref: %s', ...
              WPSS_BAND_LOW, WPSS_BAND_HIGH, REFERENCE_METHOD), ...
      'FontSize', 14, 'FontWeight', 'bold');
legend('Location', 'northeast', 'FontSize', 11);
grid on;
set(gca, 'FontSize', 12);
hold off;

fname_1d = sprintf('WPSS_GrandAverage_Band_%d_%dHz_%s.png', ...
                   WPSS_BAND_LOW, WPSS_BAND_HIGH, REFERENCE_METHOD);
saveas(fig_wpss1d, fname_1d);
fprintf('✓ Saved: %s\n', fname_1d);

%% WPSS Figures (C4, all speeds) – staggered layout: 1+3 top, 5+7+10 bottom
fig_w = figure('Position',[100,100,1800,900],'Color','white');

%% Panel geometry (normalised figure coordinates)
pw   = 0.17;    % panel width
ph   = 0.34;    % panel height
cb_w = 0.015;   % colorbar width
cb_g = 0.012;   % gap panel → colorbar
grpW = pw + cb_g + cb_w;  % total width per speed group
hgap = 0.055;   % horizontal gap between groups
vgap = 0.15;    % vertical gap between rows

%% Bottom row x-positions  (5.0 | 7.0 | 10.0 m/s)
% Centre the 3 groups in the figure
bot_left = (1 - (3*grpW + 2*hgap)) / 2;
xB = [bot_left, bot_left + grpW + hgap, bot_left + 2*(grpW + hgap)];
yB = 0.08;

%% Top row x-positions  (1.0 | 3.0 m/s, centred between bottom pairs)
% 1.0 m/s panel centre = midpoint between 5.0 and 7.0 panel centres
cx_1 = (xB(1)+pw/2 + xB(2)+pw/2) / 2;
cx_3 = (xB(2)+pw/2 + xB(3)+pw/2) / 2;
xT   = [cx_1 - pw/2, cx_3 - pw/2];
yT   = yB + ph + vgap;

%% Speed → row assignment
%  speeds(1)=1.0, speeds(2)=3.0  → top row   (index 1,2 in xT)
%  speeds(3)=5.0, speeds(4)=7.0, speeds(5)=10.0 → bottom row (index 1,2,3 in xB)
for s = 1:n_speeds
    sl = sprintf('speed_%d', round(unique_speeds(s)*10));

    if s <= 2        % top row
        xP = xT(s);  yP = yT;
    else             % bottom row
        xP = xB(s-2); yP = yB;
    end

    ax = axes('Position', [xP, yP, pw, ph]); 
    imagesc(time_ms, center_freqs, WPSS_results.C4.(sl).wpss);
    axis xy; colormap('jet'); caxis([0 0.35]);
    colorbar('Position', [xP+pw+cb_g, yP, cb_w, ph]);
    xline(0,'--w','LineWidth',2);
    xlim([-200 800]);
    xlabel('Time (ms)','FontSize',11,'FontWeight','bold');
    ylabel('Frequency (Hz)','FontSize',11,'FontWeight','bold');
    title(sprintf('%.1f m/s (n=%d)', unique_speeds(s), WPSS_results.C4.(sl).n_trials), ...
          'FontSize',12,'FontWeight','bold');
    set(gca,'FontSize',10);
end

sgtitle(sprintf('WPSS – C4 | Ref: %s', REFERENCE_METHOD),'FontSize',15,'FontWeight','bold');
fname = sprintf('WPSS_C4_AllSpeeds_%s.png', REFERENCE_METHOD);
saveas(fig_w, fname);  fprintf('✓ Saved: %s\n', fname);

%% WPSS Component Comparison Figure
fig_wpss_comparison = figure('Position', [100, 100, 1600, 600], 'Color', 'white');

% Subplot 1: Mean WPSS by speed (grouped by component)
subplot(1, 3, 1);
bar_data = [wpss_p1, wpss_n1, wpss_p2];
b = bar(unique_speeds, bar_data, 'grouped');
b(1).FaceColor = [0.2 0.6 0.9];  % P1 blue
b(2).FaceColor = [0.9 0.3 0.2];  % N1 red
b(3).FaceColor = [0.3 0.8 0.3];  % P2 green
xlabel('Speed (m/s)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Mean WPSS', 'FontSize', 12, 'FontWeight', 'bold');
title('WPSS by Component and Speed', 'FontSize', 13, 'FontWeight', 'bold');
legend('P1 (45-90ms)', 'N1 (100-160ms)', 'P2 (200-225ms)', 'Location', 'best', 'FontSize', 10);
ylim([0 0.3]);
grid on;
set(gca, 'FontSize', 11);

% Subplot 2: Mean WPSS by component (averaged across speeds)
subplot(1, 3, 2);
mean_wpss = [mean(wpss_p1), mean(wpss_n1), mean(wpss_p2)];
std_wpss  = [std(wpss_p1),  std(wpss_n1),  std(wpss_p2)];
bar_handles = bar(1:3, mean_wpss, 'FaceColor', 'flat');
bar_handles.CData(1,:) = [0.2 0.6 0.9];  % P1
bar_handles.CData(2,:) = [0.9 0.3 0.2];  % N1
bar_handles.CData(3,:) = [0.3 0.8 0.3];  % P2
hold on;
errorbar(1:3, mean_wpss, std_wpss, 'k.', 'LineWidth', 2, 'MarkerSize', 15);
hold off;
set(gca, 'XTick', 1:3, 'XTickLabel', {'P1', 'N1', 'P2'});
ylabel('Mean WPSS', 'FontSize', 12, 'FontWeight', 'bold');
title('Average WPSS by Component', 'FontSize', 13, 'FontWeight', 'bold');
ylim([0 0.3]);
grid on;
set(gca, 'FontSize', 11);

% Subplot 3: Heatmap (Component x Speed)
subplot(1, 3, 3);
wpss_matrix = [wpss_p1'; wpss_n1'; wpss_p2'];
imagesc(unique_speeds, 1:3, wpss_matrix);
colorbar; colormap('hot'); caxis([0 0.3]);
set(gca, 'YTick', 1:3, 'YTickLabel', {'P1', 'N1', 'P2'});
xlabel('Speed (m/s)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Component', 'FontSize', 12, 'FontWeight', 'bold');
title('WPSS Heatmap', 'FontSize', 13, 'FontWeight', 'bold');
% Annotate cells
for c = 1:3
    for s = 1:n_speeds
        text(unique_speeds(s), c, sprintf('%.2f', wpss_matrix(c, s)), ...
             'HorizontalAlignment', 'center', 'FontSize', 10, ...
             'FontWeight', 'bold', 'Color', 'white');
    end
end
set(gca, 'FontSize', 11);

sgtitle(sprintf('WPSS Component Comparison - C4 | Ref: %s', REFERENCE_METHOD), ...
        'FontSize', 16, 'FontWeight', 'bold');

wpss_comparison_filename = sprintf('WPSS_Component_Comparison_%s.png', REFERENCE_METHOD);
saveas(fig_wpss_comparison, wpss_comparison_filename);
fprintf('✓ Saved: %s\n\n', wpss_comparison_filename);

fprintf('✓ COMPLETE ANALYSIS FINISHED\n');
fprintf('  Participants: %d\n', n_participants);
for i = 1:n_participants; fprintf('    - %s\n', successful_participants{i}); end
fprintf('  Data files:   %s  |  %s\n', group_file, balanced_file);
fprintf('  Peak file:    %s\n', peak_file);
fprintf('  WPSS file:    %s\n\n', wpss_file);

%% HELPER FUNCTION

function result = iif(condition, true_val, false_val)
    if condition; result = true_val; else; result = false_val; end
end