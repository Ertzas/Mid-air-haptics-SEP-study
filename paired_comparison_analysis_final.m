%% paired_comparison_analysis.m
%  Analysis of paired comparison data (Intensity & Valence)
%  CSV format: trial,block,participant,speed_first,speed_second,
%              choice,chosen_speed,rt_ms,timestamp
%  choice: 1=Stim1, 2=Stim2, -1=No opinion

clear; clc; close all;

%% ========================================================================
%  CONFIGURATION
%  ========================================================================
dataDir = fullfile(pwd, 'data', 'csv_files_subjective');
speeds  = [1.0, 3.0, 5.0, 7.0, 10.0];
nSpeeds = length(speeds);
speedLabels = arrayfun(@(s) sprintf('%.1f', s), speeds, 'UniformOutput', false);

%% ========================================================================
%  1. LOAD FILES
%  ========================================================================
fprintf('========================================================\n');
fprintf('  PAIRED COMPARISON ANALYSIS\n');
fprintf('========================================================\n\n');

intFiles = dir(fullfile(dataDir, '*_intensity_*.csv'));
valFiles = dir(fullfile(dataDir, '*_valence_*.csv'));

fprintf('  Found: %d intensity files, %d valence files\n\n', ...
    length(intFiles), length(valFiles));

[intData, intPIDs] = loadAllCSV(intFiles, dataDir, 'Intensity');
[valData, valPIDs] = loadAllCSV(valFiles, dataDir, 'Valence');

allPIDs = unique([intPIDs; valPIDs]);
fprintf('  Total participants: %d\n', length(allPIDs));
fprintf('  Intensity: %s\n', strjoin(intPIDs, ', '));
fprintf('  Valence:   %s\n\n', strjoin(valPIDs, ', '));

%% ========================================================================
%  2. INTENSITY ANALYSIS
%  ========================================================================
intResults = [];
if ~isempty(intData)
    fprintf('========================================================\n');
    fprintf('  INTENSITY ANALYSIS\n');
    fprintf('========================================================\n\n');
    intResults = analyzeData(intData, intPIDs, speeds, 'Intensity');
end

%% ========================================================================
%  3. VALENCE ANALYSIS
%  ========================================================================
valResults = [];
if ~isempty(valData)
    fprintf('\n========================================================\n');
    fprintf('  VALENCE ANALYSIS\n');
    fprintf('========================================================\n\n');
    valResults = analyzeData(valData, valPIDs, speeds, 'Valence');
end

%% ========================================================================
%  4. COMPARISON PLOTS
%  ========================================================================
if ~isempty(intResults) && ~isempty(valResults)
    fprintf('\n========================================================\n');
    fprintf('  COMPARISON: INTENSITY vs VALENCE\n');
    fprintf('========================================================\n\n');
    plotComparison(intResults, valResults, speeds, speedLabels);
end

fprintf('\n  Analysis complete.\n');

%% ========================================================================
%  FUNCTIONS
%  ========================================================================

function [allData, pids] = loadAllCSV(files, dataDir, label)
    allData = [];
    pids = {};

    for f = 1:length(files)
        fpath = fullfile(dataDir, files(f).name);
        fprintf('  Loading: %s\n', files(f).name);

        T = readtable(fpath, 'Delimiter', ',', 'TextType', 'string');

        % Extract PID from filename
        tokens = regexp(files(f).name, '(P\d+)', 'tokens');
        if ~isempty(tokens)
            pid = string(tokens{1}{1});
        else
            pid = string(T.participant(1));
        end

        % PID column as string
        T.pid = repmat(pid, height(T), 1);

        if isempty(allData)
            allData = T;
        else
            allData = [allData; T];
        end

        if ~ismember(pid, pids)
            pids{end+1} = char(pid);
        end
    end

    pids = pids(:);

    if ~isempty(allData)
        allData.pid = string(allData.pid);

        nTotal = height(allData);
        nNoOp  = sum(allData.choice == -1);
        fprintf('  %s: %d trials loaded, %d no opinion (%.1f%%)\n\n', ...
            label, nTotal, nNoOp, 100*nNoOp/nTotal);
    end
end

function results = analyzeData(data, pids, speeds, label)
    nSpeeds = length(speeds);
    speedLabels = arrayfun(@(s) sprintf('%.1f', s), speeds, 'UniformOutput', false);

    data.pid = string(data.pid);

    % Filter "no opinion"
    validData = data(data.choice ~= -1, :);
    noOpData  = data(data.choice == -1, :);
    nValid    = height(validData);
    nNoOp     = height(noOpData);

    fprintf('  Valid trials:  %d (%.1f%%)\n', nValid, 100*nValid/height(data));
    fprintf('  No opinion:    %d (%.1f%%)\n\n', nNoOp, 100*nNoOp/height(data));

    % -------------------------------------------------------------------
    %  WIN MATRIX (overall)
    % -------------------------------------------------------------------
    W = buildWinMatrix(validData, speeds);

    fprintf('  Win matrix (overall, row beats column):\n');
    printMatrix(W, speedLabels);

    % -------------------------------------------------------------------
    %  BRADLEY-TERRY MODEL
    % -------------------------------------------------------------------
    [bt, btSE] = fitBradleyTerry(W);

    fprintf('\n  Bradley-Terry scale values:\n');
    for i = 1:nSpeeds
        fprintf('    %5s m/s:  %+.3f (SE: %.3f)\n', speedLabels{i}, bt(i), btSE(i));
    end

    % -------------------------------------------------------------------
    %  WIN RATE PER SPEED
    % -------------------------------------------------------------------
    fprintf('\n  Win rate per speed:\n');
    winRates = zeros(nSpeeds, 1);
    for i = 1:nSpeeds
        wins   = sum(W(i,:));
        losses = sum(W(:,i));
        total  = wins + losses;
        winRates(i) = wins / max(total, 1);
        fprintf('    %5s m/s: %3d/%3d (%.1f%%)\n', speedLabels{i}, wins, total, 100*winRates(i));
    end

    % -------------------------------------------------------------------
    %  PER PARTICIPANT
    % -------------------------------------------------------------------
    fprintf('\n  Per-participant analysis:\n');
    btPerSubj = zeros(length(pids), nSpeeds);

    for p = 1:length(pids)
        pidStr = string(pids{p});
        pValid = validData(validData.pid == pidStr, :);
        pAll   = data(data.pid == pidStr, :);

        if isempty(pValid)
            fprintf('    %s: no valid trials\n', pids{p});
            continue;
        end

        Wp = buildWinMatrix(pValid, speeds);
        [btp, ~] = fitBradleyTerry(Wp);
        btPerSubj(p,:) = btp;

        % Circular triads
        ct = circularTriads(Wp);
        maxCT = (nSpeeds^3 - nSpeeds) / 24;  % n=5 -> max=5

        % No opinion
        pNoOp  = sum(pAll.choice == -1);
        pTotal = height(pAll);

        fprintf('    %s: %d trials, %d no opinion (%.1f%%), circ. triads: %d/%d (%.1f%%)\n', ...
            pids{p}, height(pValid), pNoOp, 100*pNoOp/max(pTotal,1), ...
            ct, maxCT, 100*ct/max(maxCT,1));
    end

    % -------------------------------------------------------------------
    %  "NO OPINION" PER PAIR
    % -------------------------------------------------------------------
    if nNoOp > 0
        fprintf('\n  "No opinion" per pair:\n');
        fprintf('    %-15s  NoOp  Total  Rate\n', 'Pair');
        for i = 1:nSpeeds
            for j = i+1:nSpeeds
                mask1 = abs(data.speed_first - speeds(i)) < 0.01 & ...
                        abs(data.speed_second - speeds(j)) < 0.01;
                mask2 = abs(data.speed_first - speeds(j)) < 0.01 & ...
                        abs(data.speed_second - speeds(i)) < 0.01;
                tot = sum(mask1) + sum(mask2);
                nop = sum((mask1 | mask2) & data.choice == -1);
                if tot > 0
                    fprintf('    %5s vs %5s:  %3d / %3d  (%.1f%%)\n', ...
                        speedLabels{i}, speedLabels{j}, nop, tot, 100*nop/tot);
                end
            end
        end
    end

    % -------------------------------------------------------------------
    %  POSITION BIAS
    % -------------------------------------------------------------------
    nStim1 = sum(validData.choice == 1);
    nStim2 = sum(validData.choice == 2);

    z = (nStim1 - nValid*0.5) / sqrt(nValid * 0.25);
    pBinom = 2 * (1 - normcdf(abs(z)));
    fprintf('\n  Position bias: Stim1=%d (%.1f%%) vs Stim2=%d (%.1f%%), p=%.4f %s\n', ...
        nStim1, 100*nStim1/nValid, nStim2, 100*nStim2/nValid, pBinom, ...
        ternary(pBinom < 0.05, '***', 'n.s.'));

    % -------------------------------------------------------------------
    %  PLOTS
    % -------------------------------------------------------------------
    plotResults(bt, btSE, btPerSubj, pids, W, winRates, speeds, speedLabels, label);

    % -------------------------------------------------------------------
    %  RESULTS STRUCT
    % -------------------------------------------------------------------
    results.bt        = bt;
    results.btSE      = btSE;
    results.btPerSubj = btPerSubj;
    results.W         = W;
    results.winRates  = winRates;
    results.pids      = pids;
    results.nNoOp     = nNoOp;
    results.nValid    = nValid;
end

function W = buildWinMatrix(data, speeds)
    n = length(speeds);
    W = zeros(n);

    for row = 1:height(data)
        sf = data.speed_first(row);
        ss = data.speed_second(row);
        ch = data.choice(row);

        if ch == -1, continue; end

        idxF = find(abs(speeds - sf) < 0.01, 1);
        idxS = find(abs(speeds - ss) < 0.01, 1);

        if isempty(idxF) || isempty(idxS), continue; end

        if ch == 1
            W(idxF, idxS) = W(idxF, idxS) + 1;
        else
            W(idxS, idxF) = W(idxS, idxF) + 1;
        end
    end
end

function [params, se] = fitBradleyTerry(W)
    n = size(W, 1);
    params = ones(n, 1);

    for iter = 1:500
        old = params;
        for i = 1:n
            num = 0; den = 0;
            for j = 1:n
                if i == j, continue; end
                nij = W(i,j) + W(j,i);
                if nij == 0, continue; end
                num = num + W(i,j);
                den = den + nij / (params(i) + params(j));
            end
            if den > 0
                params(i) = num / den;
            end
        end
        params = params / sum(params);
        if max(abs(params - old)) < 1e-10, break; end
    end

    logP = log(params);
    logP = logP - mean(logP);

    se = zeros(n, 1);
    for i = 1:n
        info = 0;
        for j = 1:n
            if i == j, continue; end
            nij = W(i,j) + W(j,i);
            if nij == 0, continue; end
            pij = params(i) / (params(i) + params(j));
            info = info + nij * pij * (1 - pij);
        end
        if info > 0
            se(i) = 1 / sqrt(info);
        end
    end

    params = logP;
end

function ct = circularTriads(W)
    n = size(W, 1);
    ct = 0;
    D = zeros(n);
    for i = 1:n
        for j = 1:n
            if i ~= j
                D(i,j) = double(W(i,j) > W(j,i));
            end
        end
    end
    for i = 1:n
        for j = i+1:n
            for k = j+1:n
                if (D(i,j) && D(j,k) && D(k,i)) || ...
                   (D(j,i) && D(k,j) && D(i,k))
                    ct = ct + 1;
                end
            end
        end
    end
end

function s = ternary(cond, a, b)
    if cond, s = a; else, s = b; end
end

function printMatrix(W, labels)
    n = size(W, 1);
    header = '           ';
    for j = 1:n
        header = [header, sprintf('%6s', labels{j})];
    end
    fprintf('  %s\n', header);
    for i = 1:n
        row = sprintf('    %5s  ', labels{i});
        for j = 1:n
            if i == j
                row = [row, '     -'];
            else
                row = [row, sprintf('%6d', W(i,j))];
            end
        end
        fprintf('  %s\n', row);
    end
end

function plotResults(bt, btSE, btPerSubj, pids, W, winRates, speeds, labels, titleStr)
    nS = length(speeds);

    figure('Name', titleStr, 'Position', [100 100 1200 500]);

    % 1. Bradley-Terry scale values
    subplot(1,3,1);
    errorbar(1:nS, bt, btSE, 'ko-', 'LineWidth', 1.5, ...
        'MarkerFaceColor', [0.2 0.4 0.8], 'MarkerSize', 8);
    set(gca, 'XTick', 1:nS, 'XTickLabel', labels);
    xlabel('Drawing Speed (m/s)');
    ylabel('BT Scale Value (log)');
    title(sprintf('%s: Bradley-Terry', titleStr));
    grid on;

    % 2. Win rate
    subplot(1,3,2);
    bar(winRates, 'FaceColor', [0.3 0.6 0.3]);
    set(gca, 'XTick', 1:nS, 'XTickLabel', labels);
    xlabel('Drawing Speed (m/s)');
    ylabel('Win Rate');
    title(sprintf('%s: Win Rate', titleStr));
    ylim([0 1]); yline(0.5, '--', 'Color', [0.5 0.5 0.5]);
    grid on;

    % 3. Per-participant BT values (individual data only, no overall mean)
    subplot(1,3,3);
    hold on;
    colors = lines(length(pids));
    for p = 1:length(pids)
        if any(btPerSubj(p,:) ~= 0)
            plot(1:nS, btPerSubj(p,:), 'o-', 'Color', colors(p,:), ...
                'LineWidth', 1, 'MarkerSize', 5, 'DisplayName', pids{p});
        end
    end
    hold off;
    set(gca, 'XTick', 1:nS, 'XTickLabel', labels);
    xlabel('Drawing Speed (m/s)'); ylabel('BT Scale Value');
    title(sprintf('%s: Per Participant', titleStr));
    legend('Location', 'best', 'FontSize', 7);
    grid on;

    sgtitle(sprintf('%s — %d Participants', titleStr, length(pids)), ...
        'FontSize', 14, 'FontWeight', 'bold');
end

function plotComparison(intR, valR, speeds, labels)
    nS = length(speeds);
    figure('Name', 'Comparison', 'Position', [150 150 900 400]);

    % BT values
    subplot(1,2,1);
    hold on;
    errorbar((1:nS)-0.1, intR.bt, intR.btSE, 'rs-', 'LineWidth', 1.5, ...
        'MarkerFaceColor', [0.8 0.2 0.2], 'MarkerSize', 8, 'DisplayName', 'Intensity');
    errorbar((1:nS)+0.1, valR.bt, valR.btSE, 'go-', 'LineWidth', 1.5, ...
        'MarkerFaceColor', [0.2 0.7 0.2], 'MarkerSize', 8, 'DisplayName', 'Valence');
    hold off;
    set(gca, 'XTick', 1:nS, 'XTickLabel', labels);
    xlabel('Drawing Speed (m/s)'); ylabel('BT Scale Value (log)');
    title('Intensity vs Valence: Bradley-Terry');
    legend('Location', 'best'); grid on;

    % Win rates
    subplot(1,2,2);
    b = bar([intR.winRates, valR.winRates]);
    b(1).DisplayName = 'Intensity';
    b(2).DisplayName = 'Valence';
    set(gca, 'XTick', 1:nS, 'XTickLabel', labels);
    xlabel('Drawing Speed (m/s)'); ylabel('Win Rate');
    title('Intensity vs Valence: Win Rate');
    legend('Location', 'best');
    ylim([0 1]); yline(0.5, '--', 'Color', [0.5 0.5 0.5], 'HandleVisibility', 'off'); grid on;

    sgtitle('Comparison: Intensity vs Valence', 'FontSize', 14, 'FontWeight', 'bold');

    % Correlation
    [r, p] = corr(intR.bt, valR.bt);
    fprintf('  BT correlation (Intensity vs Valence): r=%.3f, p=%.4f\n', r, p);

    fprintf('\n  Summary:\n');
    fprintf('    %-10s  %-12s  %-12s\n', 'Speed', 'Int (BT)', 'Val (BT)');
    for i = 1:nS
        fprintf('    %5s m/s   %+.3f        %+.3f\n', labels{i}, intR.bt(i), valR.bt(i));
    end
end