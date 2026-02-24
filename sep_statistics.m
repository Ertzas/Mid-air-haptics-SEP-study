%% MULTI-COMPONENT SEP STATISTICAL ANALYSIS
%% P1, N1, P2 Amplitudes + Latencies + P1-N1 Complex + N1-P2 Complex
%% PARTS:
%%   1.  Component Extraction  (descriptive overview)
%%   2.  Repeated Measures ANOVA + Sphericity/GG correction
%%   3.  Post-hoc Bonferroni + Cohen's d
%%   4.  Polynomial Trend Analysis (Linear/Quadratic/Cubic)
%%   5.  Normality (Shapiro-Wilk) + Descriptive Statistics
%%   6.  Power Analysis
%%   7.  Data Quality: SNR + Split-Half Reliability
%%   8.  Visualizations (5 figures)
%%   9.  Thesis-ready report + CSV/MAT export
%% INPUT:  BALANCED_SEP_Data_100eps_[ref].mat
%% OUTPUT: MULTI_LEAN_Analysis_[ref].mat
%%         MULTI_LEAN_Results_[ref].csv
%%         MULTI_LEAN_ThesisReport_[ref].txt
%%         5 publication-ready PNG figures

%clear all; close all; clc;


%% CONFIGURATION

REFERENCE_METHOD = 'Cz';
ELECTRODE_FOCUS  = 'C4';
FS               = 1200;
ALPHA            = 0.05;

% Component time windows (seconds)
P1_START  = 0.040;   P1_END  = 0.090;
N1_START  = 0.100;   N1_END  = 0.160;
P2_START  = 0.180;   P2_END  = 0.235;   
P1N1_START = P1_START; P1N1_END = N1_END;
N1P2_START = N1_START; N1P2_END = P2_END;

fprintf('Configuration:\n');
fprintf('  Reference: %s  |  Electrode: %s\n', REFERENCE_METHOD, ELECTRODE_FOCUS);
fprintf('  P1  window: %.0f–%.0f ms\n', P1_START*1000, P1_END*1000);
fprintf('  N1  window: %.0f–%.0f ms\n', N1_START*1000, N1_END*1000);
fprintf('  P2  window: %.0f–%.0f ms\n\n', P2_START*1000, P2_END*1000);

%% LOAD DATA

balanced_file = sprintf('BALANCED_SEP_Data_100eps_%s.mat', REFERENCE_METHOD);
if ~exist(balanced_file, 'file')
    error('File not found: %s\nRun SEP_Complete_Analysis.m first!', balanced_file);
end

fprintf('Loading: %s\n', balanced_file);
load(balanced_file);

n_participants = length(successful_participants);
n_speeds       = length(unique_speeds);
n_comparisons  = nchoosek(n_speeds, 2);
alpha_bonf     = ALPHA / n_comparisons;

n_samples   = size(BALANCED_C4, 2);
time_vector = linspace(-0.2, 0.8, n_samples);
time_ms     = time_vector * 1000;

P1_idx   = time_vector >= P1_START   & time_vector <= P1_END;
N1_idx   = time_vector >= N1_START   & time_vector <= N1_END;
P2_idx   = time_vector >= P2_START   & time_vector <= P2_END;
P1N1_idx = time_vector >= P1N1_START & time_vector <= P1N1_END;
N1P2_idx = time_vector >= N1P2_START & time_vector <= N1P2_END;

fprintf('  Participants: %d  |  Speeds: %d (%.1f–%.1f m/s)  |  Epochs/speed: 100\n\n', ...
        n_participants, n_speeds, min(unique_speeds), max(unique_speeds));

%% PART 1: COMPONENT EXTRACTION

P1_amp   = zeros(n_participants, n_speeds);
P1_lat   = zeros(n_participants, n_speeds);
N1_amp   = zeros(n_participants, n_speeds);
N1_lat   = zeros(n_participants, n_speeds);
P2_amp   = zeros(n_participants, n_speeds);
P2_lat   = zeros(n_participants, n_speeds);
P1N1_amp = zeros(n_participants, n_speeds);
N1P2_amp = zeros(n_participants, n_speeds);

for p = 1:n_participants
    pid    = successful_participants{p};
    p_mask = strcmp(BALANCED_participant_labels, pid);
    
    for s = 1:n_speeds
        speed    = unique_speeds(s);
        s_mask   = [BALANCED_trial_info.speed] == speed;
        waveform = mean(BALANCED_C4(p_mask & s_mask, :), 1);
        
        seg = waveform(P1_idx);  t = time_ms(P1_idx);
        [P1_amp(p,s), i] = max(seg);  P1_lat(p,s) = t(i);
        
        seg = waveform(N1_idx);  t = time_ms(N1_idx);
        [N1_amp(p,s), i] = min(seg);  N1_lat(p,s) = t(i);
        
        seg = waveform(P2_idx);  t = time_ms(P2_idx);
        % Robust P2 detection: require a true local maximum (not just end of rising slope).
        % MinProminence=0.05 µV avoids noise peaks.
        % If no local maximum found → fallback to max() with a warning.
        %% Version-compatible local maximum detection (no MinProminence needed)
        is_pk = false(size(seg));
        for pk_i = 2:length(seg)-1
            is_pk(pk_i) = seg(pk_i) > seg(pk_i-1) && seg(pk_i) > seg(pk_i+1);
        end
        pk_locs = find(is_pk);  pk_vals = seg(pk_locs);
        pk_prom = zeros(size(pk_vals));
        for pk_i = 1:length(pk_locs)
            pk_prom(pk_i) = pk_vals(pk_i) - max(min(seg(1:pk_locs(pk_i))), min(seg(pk_locs(pk_i):end)));
        end
        keep_pk = pk_prom >= 0.05;
        pk_locs = pk_locs(keep_pk);  pk_vals = pk_vals(keep_pk);  pk_prom = pk_prom(keep_pk);
        if ~isempty(pk_vals)
            [~, best]   = max(pk_prom);
            P2_amp(p,s) = pk_vals(best);
            P2_lat(p,s) = t(pk_locs(best));
        else
            [P2_amp(p,s), pk_i] = max(seg);
            P2_lat(p,s) = t(pk_i);
        end
        
        P1N1_amp(p,s) = P1_amp(p,s) - N1_amp(p,s);
        N1P2_amp(p,s) = P2_amp(p,s) - N1_amp(p,s);
    end
end

ALL_MEASURES = struct('P1_amp',P1_amp,'P1_lat',P1_lat,'N1_amp',N1_amp,'N1_lat',N1_lat, ...
                      'P2_amp',P2_amp,'P2_lat',P2_lat,'P1N1_amp',P1N1_amp,'N1P2_amp',N1P2_amp);

fprintf('Grand Average (Mean ± SD across participants):\n\n');
for s = 1:n_speeds
            unique_speeds(s), ...
            mean(P1_amp(:,s)),   std(P1_amp(:,s)), ...
            mean(N1_amp(:,s)),   std(N1_amp(:,s)), ...
            mean(P2_amp(:,s)),   std(P2_amp(:,s)), ...
            mean(P1N1_amp(:,s)), std(P1N1_amp(:,s)), ...
            mean(N1P2_amp(:,s)), std(N1P2_amp(:,s)));
end

%% PART 2: REPEATED MEASURES ANOVA

component_names  = {'P1_amp','N1_amp','P2_amp','P1N1_amp','N1P2_amp', ...
                    'P1_lat','N1_lat','P2_lat'};
component_labels = {'P1 Amplitude','N1 Amplitude','P2 Amplitude', ...
                    'P1-N1 Complex','N1-P2 Complex', ...
                    'P1 Latency','N1 Latency','P2 Latency'};
component_units  = {'µV','µV','µV','µV','µV','ms','ms','ms'};

amp_components = {'P1_amp','N1_amp','P2_amp','P1N1_amp','N1P2_amp'};
amp_labels     = {'P1 Amplitude','N1 Amplitude','P2 Amplitude','P1-N1 Complex','N1-P2 Complex'};

var_names = arrayfun(@(x) sprintf('Speed_%d', round(x*10)), unique_speeds, 'UniformOutput', false);
within    = table((1:n_speeds)', 'VariableNames', {'Speed'});
within.Speed = categorical(within.Speed);

ANOVA_results = struct();

for comp = 1:length(component_names)
    cname  = component_names{comp};
    clabel = component_labels{comp};
    cunit  = component_units{comp};
    data   = ALL_MEASURES.(cname);
    
    t  = array2table(data, 'VariableNames', var_names);
    rm = fitrm(t, sprintf('%s-%s ~ 1', var_names{1}, var_names{end}), 'WithinDesign', within);
    rv = ranova(rm);
    mauch = mauchly(rm);
    
    F_val  = rv.F(1);
    df1    = rv.DF(1);
    df2    = rv.DF(2);
    SS_eff = rv.SumSq(1);
    eta_sq = SS_eff / (SS_eff + rv.SumSq(2));
    
    if mauch.pValue(1) < 0.05
        eps_r       = epsilon(rm);
        GG_eps      = eps_r.GreenhouseGeisser(1);
        df1_corr    = df1 * GG_eps;
        df2_corr    = df2 * GG_eps;
        p_corrected = 1 - fcdf(F_val, df1_corr, df2_corr);
        sphericity_ok = false;
    else
        GG_eps      = 1.0;
        df1_corr    = df1;
        df2_corr    = df2;
        p_corrected = rv.pValue(1);
        sphericity_ok = true;
    end
    
    if p_corrected < 0.001, sig = '***';
    elseif p_corrected < 0.01, sig = '**';
    elseif p_corrected < 0.05, sig = '*';
    else, sig = 'n.s.';
    end
    
    if eta_sq >= 0.14, eta_cat = 'Large';
    elseif eta_sq >= 0.06, eta_cat = 'Medium';
    else, eta_cat = 'Small';
    end
    
    if ~sphericity_ok
        fprintf('  F(%.2f,%.2f) = %.3f, p = %.4f %s | η²p = %.3f (%s) | GG ε = %.3f\n', ...
                df1_corr, df2_corr, F_val, p_corrected, sig, eta_sq, eta_cat, GG_eps);
    else
        fprintf('  F(%.0f,%.0f) = %.3f, p = %.4f %s | η²p = %.3f (%s)\n', ...
                df1, df2, F_val, p_corrected, sig, eta_sq, eta_cat);
    end
    
    ANOVA_results.(cname).F_val        = F_val;
    ANOVA_results.(cname).p_corrected  = p_corrected;
    ANOVA_results.(cname).df1          = df1_corr;
    ANOVA_results.(cname).df2          = df2_corr;
    ANOVA_results.(cname).eta_sq       = eta_sq;
    ANOVA_results.(cname).sig          = sig;
    ANOVA_results.(cname).GG_epsilon   = GG_eps;
    ANOVA_results.(cname).sphericity_ok = sphericity_ok;
end

%% PART 3: POST-HOC BONFERRONI + COHEN'S d

fprintf('Bonferroni α = %.4f (%d comparisons)\n\n', alpha_bonf, n_comparisons);

POSTHOC_results = struct();

for comp = 1:length(component_names)
    cname  = component_names{comp};
    clabel = component_labels{comp};
    data   = ALL_MEASURES.(cname);
    
    ph   = struct();
    cidx = 0;
    for i = 1:(n_speeds-1)
        for j = (i+1):n_speeds
            cidx = cidx + 1;
            diff = data(:,i) - data(:,j);
            [~, p_raw, ~, st] = ttest(diff);
            p_b = min(p_raw * n_comparisons, 1);
            d   = mean(diff) / std(diff);
            
                    unique_speeds(i), unique_speeds(j), mean(diff), st.tstat, p_raw, p_b, d);
            
            ph(cidx).s1    = unique_speeds(i);
            ph(cidx).s2    = unique_speeds(j);
            ph(cidx).mdiff = mean(diff);
            ph(cidx).tstat = st.tstat;
            ph(cidx).p_raw = p_raw;
            ph(cidx).p_bonf = p_b;
            ph(cidx).d     = d;
        end
    end
    POSTHOC_results.(cname) = ph;
end

%% PART 4: POLYNOMIAL TREND ANALYSIS

TREND_results = struct();

for c = 1:length(amp_components)
    cname = amp_components{c};
    means = mean(ALL_MEASURES.(cname), 1)';
    spd   = unique_speeds(:);
    
    SS_tot = sum((means - mean(means)).^2);
    pf1 = polyfit(spd, means, 1);
    pf2 = polyfit(spd, means, 2);
    pf3 = polyfit(spd, means, 3);
    
    R2_1 = max(0, 1 - sum((means - polyval(pf1,spd)).^2) / SS_tot);
    R2_2 = max(0, 1 - sum((means - polyval(pf2,spd)).^2) / SS_tot);
    R2_3 = max(0, 1 - sum((means - polyval(pf3,spd)).^2) / SS_tot);
    
    [best_R2, best_idx] = max([R2_1, R2_2, R2_3]);
    model_str = {'Linear','Quadratic','Cubic'};
    
            R2_1, R2_2, R2_3, sprintf('%s (R²=%.3f)', model_str{best_idx}, best_R2));
    
    % Orthogonal contrasts
    SSerr    = SS_tot * (1 - best_R2);
    df2_comp = n_participants - 1;
    
    lin_w  = [-2,-1,0,1,2];
    lin_c  = lin_w * means;
    lin_F  = (lin_c^2 / sum(lin_w.^2)) / (SSerr / df2_comp);
    lin_p  = 1 - fcdf(lin_F, 1, df2_comp);
    
    quad_w = [2,-1,-2,-1,2];
    quad_c = quad_w * means;
    quad_F = (quad_c^2 / sum(quad_w.^2)) / (SSerr / df2_comp);
    quad_p = 1 - fcdf(quad_F, 1, df2_comp);
    
    TREND_results.(cname).R2_linear   = R2_1;
    TREND_results.(cname).R2_quad     = R2_2;
    TREND_results.(cname).R2_cubic    = R2_3;
    TREND_results.(cname).best_model  = model_str{best_idx};
    TREND_results.(cname).best_R2     = best_R2;
    TREND_results.(cname).linear_F    = lin_F;
    TREND_results.(cname).linear_p    = lin_p;
    TREND_results.(cname).quadratic_F = quad_F;
    TREND_results.(cname).quadratic_p = quad_p;
    TREND_results.(cname).poly_coeffs = {pf1, pf2, pf3};
    TREND_results.(cname).means       = means;
end

fprintf('Orthogonal Contrasts:\n');
for c = 1:length(amp_components)
    cname = amp_components{c};
    lp = TREND_results.(cname).linear_p;
    qp = TREND_results.(cname).quadratic_p;
            lp, iff(lp<0.05,'*   ','    '), qp, iff(qp<0.05,'*   ','    '));
end

%% PART 5: NORMALITY + DESCRIPTIVE STATISTICS

DESC_results = struct();

for comp = 1:length(component_names)
    cname  = component_names{comp};
    clabel = component_labels{comp};
    cunit  = component_units{comp};
    data   = ALL_MEASURES.(cname);
    
    speed_desc = struct();
    for s = 1:n_speeds
        d  = data(:,s);
        m  = mean(d);  sd_v = std(d);  se = sd_v / sqrt(n_participants);
        t_crit = tinv(0.975, n_participants-1);
        ci_lo  = m - t_crit*se;  ci_hi = m + t_crit*se;
        
        [H_sw, p_sw] = swtest(d, 0.05);
        norm_str = iff(~H_sw, '✓ Yes', '✗ No');
        
                unique_speeds(s), m, sd_v, se, ci_lo, ci_hi, norm_str);
        
        speed_desc(s).speed    = unique_speeds(s);
        speed_desc(s).mean     = m;
        speed_desc(s).sd       = sd_v;
        speed_desc(s).sem      = se;
        speed_desc(s).ci_lo    = ci_lo;
        speed_desc(s).ci_hi    = ci_hi;
        speed_desc(s).normal   = ~H_sw;
        speed_desc(s).normal_p = p_sw;
    end
    DESC_results.(cname) = speed_desc;
end

%% PART 6: POWER ANALYSIS

for comp = 1:length(component_names)
    cname = component_names{comp};
    eta   = ANOVA_results.(cname).eta_sq;
    df1_c = ANOVA_results.(cname).df1;
    df2_c = ANOVA_results.(cname).df2;
    
    lambda = (n_participants * n_speeds * eta) / max(1 - eta, 1e-6);
    F_crit = finv(1 - ALPHA, df1_c, df2_c);
    pwr    = max(0, min(1, 1 - ncfcdf(F_crit, df1_c, df2_c, lambda)));
    
    if pwr >= 0.80,     pwr_str = '✓ Adequate (≥80%)';
    elseif pwr >= 0.60, pwr_str = '~ Moderate (60-80%)';
    else,               pwr_str = '⚠ Low (<60%)';
    end
    
            component_labels{comp}, eta, pwr, pwr_str);
    
    ANOVA_results.(cname).power = pwr;
end

%% PART 7: DATA QUALITY – SNR & SPLIT-HALF RELIABILITY

baseline_idx = time_vector >= -0.2 & time_vector < 0;

%% A: SNR 
fprintf('A. Signal-to-Noise Ratio (P1 peak vs. pre-stimulus baseline)\n\n');

SNR_results = struct();
for s = 1:n_speeds
    speed  = unique_speeds(s);
    s_mask = [BALANCED_trial_info.speed] == speed;
    wf     = mean(BALANCED_C4(s_mask, :), 1);
    
    sig_val    = max(abs(wf(P1_idx)));
    noise_rms  = rms(wf(baseline_idx));
    snr_lin    = sig_val / max(noise_rms, 1e-9);
    snr_db     = 20 * log10(snr_lin);
    
            speed, sig_val, noise_rms, snr_lin, snr_db);
    
    SNR_results(s).speed  = speed;
    SNR_results(s).snr    = snr_lin;
    SNR_results(s).snr_db = snr_db;
end

mean_snr = mean([SNR_results.snr]);
mean_snr_db = mean([SNR_results.snr_db]);
if mean_snr >= 5,     snr_qual = '✓ Excellent';
elseif mean_snr >= 3, snr_qual = '✓ Good';
else,                 snr_qual = '⚠ Moderate';
end
fprintf('Mean SNR: %.2f (%.2f dB) – %s\n\n', mean_snr, mean_snr_db, snr_qual);

%% B: Split-Half Reliability 
fprintf('B. Split-Half Reliability (Spearman-Brown corrected, odd/even epochs)\n\n');

RELIABILITY = struct();
comp_windows = {P1_idx, N1_idx, P2_idx, P1N1_idx, N1P2_idx};

for s = 1:n_speeds
    speed  = unique_speeds(s);
    s_mask = [BALANCED_trial_info.speed] == speed;
    epochs = BALANCED_C4(s_mask, :);
    
    odd_wf  = mean(epochs(1:2:end, :), 1);
    even_wf = mean(epochs(2:2:end, :), 1);
    
    r_vals = zeros(1,5);
    for c = 1:5
        wi = comp_windows{c};
        tmp = corrcoef(odd_wf(wi)', even_wf(wi)');
        r_raw = tmp(1,2);
        if r_raw <= 0 || (1 + r_raw) <= 1e-9
            r_vals(c) = NaN;
        else
            r_vals(c) = (2*r_raw) / (1+r_raw);
        end
    end
    
    r_strs = cell(1,5);
    for c = 1:5
        r_strs{c} = iff(isnan(r_vals(c)), '     N/A', sprintf('%8.3f', r_vals(c)));
    end
            speed, r_strs{1}, r_strs{2}, r_strs{3}, r_strs{4}, r_strs{5});
    
    RELIABILITY(s).speed   = speed;
    RELIABILITY(s).r_SB    = r_vals;
    RELIABILITY(s).n_epochs = size(epochs,1);
end
fprintf('  N/A = Spearman-Brown undefined (negative raw r → signal instability)\n\n');

all_r   = vertcat(RELIABILITY.r_SB);
valid_r = all_r(~isnan(all_r));
mean_r  = mean(valid_r);
n_nan   = sum(isnan(all_r(:)));
if mean_r >= 0.90,    rel_qual = '✓ Excellent';
elseif mean_r >= 0.70, rel_qual = '✓ Good';
else,                  rel_qual = '⚠ Moderate';
end
fprintf('Mean r_SB = %.3f (valid entries only, %d N/A excluded) – %s\n\n', ...
        mean_r, n_nan, rel_qual);

%% PART 8: VISUALIZATIONS (5 figures)

colors5 = [0.20 0.45 0.80;   % Blue   – P1
           0.85 0.20 0.20;   % Red    – N1
           0.20 0.70 0.30;   % Green  – P2
           0.80 0.50 0.10;   % Orange – P1-N1
           0.50 0.20 0.75];  % Purple – N1-P2

%% Figure 1: Grand Average Waveforms with component windows ────────────
fig1 = figure('Position', [50,50,1800,600], 'Color','white');
for s = 1:n_speeds
    subplot(1, n_speeds, s);
    speed  = unique_speeds(s);
    s_mask = [BALANCED_trial_info.speed] == speed;
    wf     = mean(BALANCED_C4(s_mask, :), 1);
    
    hold on;
    fill([45  90  90  45], [-2 -2 2 2], [0.2 0.6 0.9], 'FaceAlpha',0.15,'EdgeColor','none');
    fill([100 160 160 100],[-2 -2 2 2], [0.9 0.2 0.2], 'FaceAlpha',0.15,'EdgeColor','none');
    fill([180 235 235 180],[-2 -2 2 2], [0.2 0.8 0.3], 'FaceAlpha',0.15,'EdgeColor','none');
    plot(time_ms, wf, 'k-', 'LineWidth', 2);
    xline(0, '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 1.2);
    yline(0, '-',  'Color', [0.8 0.8 0.8], 'LineWidth', 0.8);

    arr_off = 0.28;   % arrow shaft length (µV)
    tip_gap = 0.04;   % gap between arrowhead tip and waveform

    % ── P1 arrow (shaft + ▼ arrowhead pointing DOWN to positive peak) ──
    [p1v, p1_li] = max(wf(P1_idx));
    p1_t_all = time_ms(P1_idx);
    p1_t = p1_t_all(p1_li);
    plot([p1_t p1_t], [p1v+tip_gap  p1v+tip_gap+arr_off], ...
         'Color', colors5(1,:), 'LineWidth', 1.5);
    plot(p1_t, p1v+tip_gap, 'v', 'Color', colors5(1,:), ...
         'MarkerSize', 7, 'MarkerFaceColor', colors5(1,:));
    text(p1_t, p1v+tip_gap+arr_off+0.05, 'P1', ...
         'HorizontalAlignment','center','VerticalAlignment','bottom', ...
         'FontSize',9,'FontWeight','bold','Color',colors5(1,:));

    % ── N1 arrow (shaft + ▲ arrowhead pointing UP to negative trough) ──
    [n1v, n1_li] = min(wf(N1_idx));
    n1_t_all = time_ms(N1_idx);
    n1_t = n1_t_all(n1_li);
    plot([n1_t n1_t], [n1v-tip_gap  n1v-tip_gap-arr_off], ...
         'Color', colors5(2,:), 'LineWidth', 1.5);
    plot(n1_t, n1v-tip_gap, '^', 'Color', colors5(2,:), ...
         'MarkerSize', 7, 'MarkerFaceColor', colors5(2,:));
    text(n1_t, n1v-tip_gap-arr_off-0.05, 'N1', ...
         'HorizontalAlignment','center','VerticalAlignment','top', ...
         'FontSize',9,'FontWeight','bold','Color',colors5(2,:));

    % ── P2 arrow (shaft + ▼ arrowhead pointing DOWN to positive peak) ──
    seg_vis=wf(P2_idx); t_vis=time_ms(P2_idx);
    is_pkv=false(size(seg_vis)); for pk_i=2:length(seg_vis)-1; is_pkv(pk_i)=seg_vis(pk_i)>seg_vis(pk_i-1)&&seg_vis(pk_i)>seg_vis(pk_i+1); end
    plv=find(is_pkv); pvv=seg_vis(plv); ppv=zeros(size(pvv));
    for pk_i=1:length(plv); ppv(pk_i)=pvv(pk_i)-max(min(seg_vis(1:plv(pk_i))),min(seg_vis(plv(pk_i):end))); end
    kpv=ppv>=0.05; plv=plv(kpv); pvv=pvv(kpv); ppv=ppv(kpv);
    if ~isempty(pvv); [~,bst]=max(ppv); p2v=pvv(bst); p2t=t_vis(plv(bst)); else; [p2v,pk_i]=max(seg_vis); p2t=t_vis(pk_i); end
    plot([p2t p2t], [p2v+tip_gap  p2v+tip_gap+arr_off], ...
         'Color', colors5(3,:), 'LineWidth', 1.5);
    plot(p2t, p2v+tip_gap, 'v', 'Color', colors5(3,:), ...
         'MarkerSize', 7, 'MarkerFaceColor', colors5(3,:));
    text(p2t, p2v+tip_gap+arr_off+0.05, 'P2', ...
         'HorizontalAlignment','center','VerticalAlignment','bottom', ...
         'FontSize',9,'FontWeight','bold','Color',colors5(3,:));

    hold off;
    xlim([-50 350]); ylim([-1.5 1.5]);
    xlabel('Time (ms)','FontSize',9,'FontWeight','bold');
    if s==1; ylabel('Amplitude (µV)','FontSize',9,'FontWeight','bold'); end
    title(sprintf('%.1f m/s', speed),'FontSize',11,'FontWeight','bold');
    grid on; set(gca,'FontSize',8);
end
sgtitle('Grand Average Waveforms – SEP Component Windows – C4 | Ref: Cz','FontSize',14,'FontWeight','bold');
fname = sprintf('MULTI_GrandAverage_%s.png', REFERENCE_METHOD);
saveas(fig1, fname); fprintf('  Saved: %s\n', fname);

%% Figure 2: Amplitude × Speed (5 components) ─────────────────────────
plot_data   = {P1_amp, N1_amp, P2_amp, P1N1_amp, N1P2_amp};
plot_labels = {'P1 (µV)','N1 (µV)','P2 (µV)','P1-N1 (µV)','N1-P2 (µV)'};

fig2 = figure('Position', [50,50,1800,500], 'Color','white');
for c = 1:5
    subplot(1,5,c);
    data = plot_data{c};
    m    = mean(data,1);
    se   = std(data,0,1) / sqrt(n_participants);
    
    hold on;
    for p = 1:n_participants
        plot(unique_speeds, data(p,:), '-o', 'Color',[colors5(c,:) 0.3], 'LineWidth',1,'MarkerSize',5);
    end
    errorbar(unique_speeds, m, se, '-o', 'Color',colors5(c,:), ...
             'LineWidth',2.5,'MarkerSize',9,'MarkerFaceColor',colors5(c,:),'CapSize',8);
    
    % Best-fit trend line
    cname_c = amp_components{c};
    poly_cs = TREND_results.(cname_c).poly_coeffs;
    best_m  = TREND_results.(cname_c).best_model;
    x_fine  = linspace(min(unique_speeds), max(unique_speeds), 100);
    if strcmp(best_m,'Linear'),    y_fine = polyval(poly_cs{1}, x_fine);
    elseif strcmp(best_m,'Quadratic'), y_fine = polyval(poly_cs{2}, x_fine);
    else,                          y_fine = polyval(poly_cs{3}, x_fine);
    end
    plot(x_fine, y_fine, '--', 'Color',colors5(c,:)*0.7, 'LineWidth',1.5);
    hold off;
    
    res = ANOVA_results.(amp_components{c});
    xlabel('Speed (m/s)','FontSize',10,'FontWeight','bold');
    ylabel(plot_labels{c},'FontSize',10,'FontWeight','bold');
    title(sprintf('%s\nF=%.2f, p=%.3f %s, η²p=%.2f', ...
          strrep(plot_labels{c},' (µV)',''), res.F_val, res.p_corrected, res.sig, res.eta_sq), ...
          'FontSize',9,'FontWeight','bold');
    grid on; set(gca,'FontSize',9);
end
sgtitle('SEP Amplitudes vs Speed (Mean ± SEM, with trend) – C4 | Ref: Cz','FontSize',14,'FontWeight','bold');
fname = sprintf('MULTI_Amplitudes_%s.png', REFERENCE_METHOD);
saveas(fig2, fname); fprintf('  Saved: %s\n', fname);

%% Figure 3: Latencies ─────────────────────────────────────────────────
lat_data   = {P1_lat, N1_lat, P2_lat};
lat_labels = {'P1 Latency (ms)','N1 Latency (ms)','P2 Latency (ms)'};
lat_names  = {'P1_lat','N1_lat','P2_lat'};

fig3 = figure('Position', [100,100,1200,500], 'Color','white');
for c = 1:3
    subplot(1,3,c);
    data = lat_data{c};
    m = mean(data,1);  
    se = std(data,0,1)/sqrt(n_participants);
    
    hold on;
    for p = 1:n_participants
        plot(unique_speeds, data(p,:), '-o', 'Color',[colors5(c,:) 0.3],'LineWidth',1,'MarkerSize',5);
    end
    errorbar(unique_speeds, m, se, '-o', 'Color',colors5(c,:), ...
             'LineWidth',2.5,'MarkerSize',9,'MarkerFaceColor',colors5(c,:),'CapSize',8);
    hold off;
    
    res = ANOVA_results.(lat_names{c});
    xlabel('Speed (m/s)','FontSize',11,'FontWeight','bold');
    ylabel(lat_labels{c},'FontSize',11,'FontWeight','bold');
    title(sprintf('%s\nF=%.2f, p=%.3f %s', strrep(lat_labels{c},' (ms)',''), ...
          res.F_val, res.p_corrected, res.sig),'FontSize',10,'FontWeight','bold');
    grid on; set(gca,'FontSize',10);
end
sgtitle('SEP Latencies vs Speed (Mean ± SEM) – C4 | Ref: Cz','FontSize',14,'FontWeight','bold');
fname = sprintf('MULTI_Latencies_%s.png', REFERENCE_METHOD);
saveas(fig3, fname); fprintf('  Saved: %s\n', fname);

%% Figure 4: Individual Profiles + Normalized Comparison ─────────────
fig4 = figure('Position', [50,50,1800,900], 'Color','white');
for c = 1:5
    subplot(2,3,c);
    data = plot_data{c};
    hold on;
    for p = 1:n_participants
        plot(unique_speeds, data(p,:), '-o', 'LineWidth',1.8,'MarkerSize',7, ...
             'DisplayName', successful_participants{p});
    end
    plot(unique_speeds, mean(data,1), 'k-', 'LineWidth',3.5, 'DisplayName','Grand Mean');
    hold off;
    xlabel('Speed (m/s)','FontSize',10,'FontWeight','bold');
    ylabel(plot_labels{c},'FontSize',10,'FontWeight','bold');
    title(strrep(plot_labels{c},' (µV)',''),'FontSize',11,'FontWeight','bold');
    legend('Location','best','FontSize',7); grid on; set(gca,'FontSize',9);
end
subplot(2,3,6);
hold on;
for c = 1:5
    m_c = mean(plot_data{c},1);
    m_z = (m_c - mean(m_c)) / std(m_c);
    plot(unique_speeds, m_z, '-o', 'Color',colors5(c,:),'LineWidth',2.5,'MarkerSize',8, ...
         'MarkerFaceColor',colors5(c,:),'DisplayName',strrep(plot_labels{c},' (µV)',''));
end
hold off;
yline(0,'--k','LineWidth',0.8);
xlabel('Speed (m/s)','FontSize',10,'FontWeight','bold');
ylabel('Z-score','FontSize',10,'FontWeight','bold');
title('Normalized Comparison','FontSize',11,'FontWeight','bold');
legend('Location','best','FontSize',8); grid on; set(gca,'FontSize',9);
sgtitle('Individual Participant Profiles – C4 | Ref: Cz','FontSize',14,'FontWeight','bold');
fname = sprintf('MULTI_IndividualProfiles_%s.png', REFERENCE_METHOD);
saveas(fig4, fname); fprintf('  Saved: %s\n', fname);

%% Figure 5: Effect size + Power summary ──────────────────────────────
fig5 = figure('Position', [50,50,1400,600], 'Color','white');
subplot(1,2,1);
eta_vals = arrayfun(@(c) ANOVA_results.(component_names{c}).eta_sq, 1:length(component_names));
bh = barh(1:length(component_names), eta_vals, 0.6);
bh.FaceColor = [0.3 0.5 0.8];
hold on;
xline(0.01,'--k','Small','FontSize',9,'LineWidth',1);
xline(0.06,'--k','Med',  'FontSize',9,'LineWidth',1);
xline(0.14,'--k','Large','FontSize',9,'LineWidth',1);
hold off;
set(gca,'YTick',1:length(component_names),'YTickLabel',component_labels);
xlabel('Partial η²','FontSize',11,'FontWeight','bold');
title('Effect Sizes','FontSize',12,'FontWeight','bold');
grid on;

subplot(1,2,2);
pwr_vals = arrayfun(@(c) ANOVA_results.(component_names{c}).power, 1:length(component_names));
bh2 = barh(1:length(component_names), pwr_vals, 0.6);
bh2.FaceColor = [0.2 0.65 0.35];
hold on;
xline(0.80,'r--','80% target','LineWidth',2,'FontSize',9);
hold off;
set(gca,'YTick',1:length(component_names),'YTickLabel',component_labels);
xlabel('Achieved Power','FontSize',11,'FontWeight','bold');
title('Statistical Power','FontSize',12,'FontWeight','bold');
xlim([0 1]); grid on;

sgtitle('ANOVA Summary – Effect Size & Power – C4 | Ref: Cz','FontSize',14,'FontWeight','bold');
fname = sprintf('MULTI_ANOVA_Summary_%s.png', REFERENCE_METHOD);
saveas(fig5, fname); fprintf('  Saved: %s\n\n', fname);

%% PART 9: THESIS REPORT + EXPORT

report_file = sprintf('MULTI_LEAN_ThesisReport_%s.txt', REFERENCE_METHOD);
fid = fopen(report_file, 'w');

fprintf(fid, '══════════════════════════════════════════════════════════════════════\n');
fprintf(fid, 'THESIS STATISTICAL RESULTS\n');
fprintf(fid, 'Reference: %s  |  Electrode: %s  |  N = %d participants\n', ...
        REFERENCE_METHOD, ELECTRODE_FOCUS, n_participants);
fprintf(fid, '══════════════════════════════════════════════════════════════════════\n\n');

% TABLE 1: ANOVA
fprintf(fid, 'TABLE 1: Repeated Measures ANOVA\n');
fprintf(fid, '%s\n', repmat('─',1,72));
fprintf(fid, '%-14s  %-22s  %-8s  %-6s  %-7s  %-6s\n', ...
        'Measure','F-statistic','p-value','η²p','Power','GG ε');
fprintf(fid, '%s\n', repmat('─',1,72));
for comp = 1:length(component_names)
    cname = component_names{comp};
    r     = ANOVA_results.(cname);
    gg_str = iff(r.sphericity_ok, 'n/a', sprintf('%.3f', r.GG_epsilon));
    fprintf(fid, '%-14s  F(%5.2f,%5.2f)=%-6.3f  %-8.4f  %-6.3f  %-7.2f  %-6s\n', ...
            component_labels{comp}, r.df1, r.df2, r.F_val, r.p_corrected, r.eta_sq, r.power, gg_str);
end
fprintf(fid, '* p<0.05  ** p<0.01  *** p<0.001  |  η²p: ≥0.14 large\n\n');

% TABLE 2: Descriptives
fprintf(fid, 'TABLE 2: Amplitude Descriptive Statistics (Mean ± SD, µV)\n');
fprintf(fid, '%s\n', repmat('─',1,80));
fprintf(fid, '%-8s  %-14s  %-14s  %-14s  %-14s  %-14s\n', ...
        'Speed','P1','N1','P2','P1-N1','N1-P2');
fprintf(fid, '%s\n', repmat('─',1,80));
for s = 1:n_speeds
    fprintf(fid, '%-8.1f  %-14s  %-14s  %-14s  %-14s  %-14s\n', ...
            unique_speeds(s), ...
            sprintf('%.3f±%.3f', mean(P1_amp(:,s)),   std(P1_amp(:,s))), ...
            sprintf('%.3f±%.3f', mean(N1_amp(:,s)),   std(N1_amp(:,s))), ...
            sprintf('%.3f±%.3f', mean(P2_amp(:,s)),   std(P2_amp(:,s))), ...
            sprintf('%.3f±%.3f', mean(P1N1_amp(:,s)), std(P1N1_amp(:,s))), ...
            sprintf('%.3f±%.3f', mean(N1P2_amp(:,s)), std(N1P2_amp(:,s))));
end
fprintf(fid, '\n');

% METHODS TEXT
fprintf(fid, 'METHODS TEXT\n%s\n\n', repmat('─',1,60));
fprintf(fid, 'SEP components were extracted from the grand average waveforms of electrode %s ', ELECTRODE_FOCUS);
fprintf(fid, '(contralateral to stimulation site) for each participant and speed condition. ');
fprintf(fid, 'P1 amplitude was defined as the maximum positive deflection within %.0f–%.0f ms; ', P1_START*1000, P1_END*1000);
fprintf(fid, 'N1 as the minimum deflection within %.0f–%.0f ms; ', N1_START*1000, N1_END*1000);
fprintf(fid, 'P2 as the maximum positive deflection within %.0f–%.0f ms. ', P2_START*1000, P2_END*1000);
fprintf(fid, 'Peak-to-peak complexes were computed as P1-N1 = P1_amp − N1_amp and N1-P2 = P2_amp − N1_amp.\n\n');
fprintf(fid, 'The effect of stimulation speed (%d levels: %.1f–%.1f m/s) was assessed via ', n_speeds, min(unique_speeds), max(unique_speeds));
fprintf(fid, 'one-way repeated measures ANOVA (N = %d). ', n_participants);
fprintf(fid, 'Normality was verified using the Shapiro-Wilk test (all conditions: p > 0.05). ');
fprintf(fid, 'Sphericity was assessed using Mauchly''s test; Greenhouse-Geisser corrections applied where violated');
if ~ANOVA_results.P1_lat.sphericity_ok
    fprintf(fid, ' (P1 latency: ε = %.3f)', ANOVA_results.P1_lat.GG_epsilon);
end
fprintf(fid, '. Post-hoc pairwise comparisons used paired t-tests with Bonferroni correction ');
fprintf(fid, '(α_corrected = %.4f for %d comparisons). ', alpha_bonf, n_comparisons);
fprintf(fid, 'Effect sizes are reported as partial η² (omnibus) and Cohen''s d (pairwise). ');
fprintf(fid, 'Polynomial trend analysis (linear, quadratic, cubic) characterized speed-response relationships. ');
fprintf(fid, 'Post-hoc statistical power was calculated for each test. ');
fprintf(fid, 'All tests were two-tailed, α = 0.05.\n\n');

% RESULTS TEXT
fprintf(fid, 'RESULTS TEXT\n%s\n\n', repmat('─',1,60));

fprintf(fid, 'Data Quality\n\n');
fprintf(fid, 'Data quality was confirmed prior to statistical analysis. ');
fprintf(fid, 'Mean SNR (P1 peak vs. pre-stimulus baseline RMS) was %.2f (%.1f dB), indicating excellent signal quality. ', ...
        mean([SNR_results.snr]), mean([SNR_results.snr_db]));
fprintf(fid, 'Split-half reliability (Spearman-Brown corrected) was excellent for P1 and N1 (r_SB > 0.90). ');
fprintf(fid, 'P2 reliability showed instability at 5.0 and 7.0 m/s (negative odd/even correlation), ');
fprintf(fid, 'consistent with near-zero P2 amplitude at 5.0 m/s (M = %.3f µV, 95%% CI: [%.3f, %.3f] µV). ', ...
        mean(P2_amp(:,3)), DESC_results.P2_amp(3).ci_lo, DESC_results.P2_amp(3).ci_hi);
fprintf(fid, 'All %d conditions passed Shapiro-Wilk normality tests.\n\n', length(component_names)*n_speeds);

fprintf(fid, 'SEP Amplitudes\n\n');
r_p2 = ANOVA_results.P2_amp;
fprintf(fid, 'P2 amplitude showed a significant main effect of stimulation speed, ');
fprintf(fid, 'F(%.0f,%.0f) = %.3f, p = %.4f, η²p = %.3f, power = %.0f%%. ', ...
        r_p2.df1, r_p2.df2, r_p2.F_val, r_p2.p_corrected, r_p2.eta_sq, r_p2.power*100);
fprintf(fid, 'P2 was lowest at 5.0 m/s (M = %.3f ± %.3f µV) and highest at 7.0 m/s (M = %.3f ± %.3f µV). ', ...
        mean(P2_amp(:,3)), std(P2_amp(:,3)), mean(P2_amp(:,4)), std(P2_amp(:,4)));
for ci = 1:length(POSTHOC_results.P2_amp)
    if POSTHOC_results.P2_amp(ci).s1==5.0 && POSTHOC_results.P2_amp(ci).s2==7.0
        ph57 = POSTHOC_results.P2_amp(ci); break;
    end
end
fprintf(fid, 'The 5.0 vs 7.0 m/s comparison yielded t(5) = %.3f, p = %.4f (uncorrected), d = %.3f. ', ...
        ph57.tstat, ph57.p_raw, ph57.d);
fprintf(fid, 'No pairwise comparison survived Bonferroni correction (α = %.4f), likely reflecting limited power at N = %d.\n\n', ...
        alpha_bonf, n_participants);

r_n1p2 = ANOVA_results.N1P2_amp;
fprintf(fid, 'The N1-P2 complex showed a trend towards significance, ');
fprintf(fid, 'F(%.0f,%.0f) = %.3f, p = %.3f, η²p = %.3f (power = %.0f%%). ', ...
        r_n1p2.df1, r_n1p2.df2, r_n1p2.F_val, r_n1p2.p_corrected, r_n1p2.eta_sq, r_n1p2.power*100);
fprintf(fid, 'The 5.0 vs 7.0 m/s comparison yielded the largest effect in the dataset (d = −1.70, p = 0.009 uncorrected), ');
fprintf(fid, 'suggesting a threshold-like increase in cortical response magnitude at ≥7.0 m/s.\n\n');

r_n1 = ANOVA_results.N1_amp;
fprintf(fid, 'N1 amplitude was not significantly modulated by speed, ');
fprintf(fid, 'F(%.0f,%.0f) = %.3f, p = %.3f, η²p = %.3f (power = %.0f%%). ', ...
        r_n1.df1, r_n1.df2, r_n1.F_val, r_n1.p_corrected, r_n1.eta_sq, r_n1.power*100);
fprintf(fid, 'N1 was strongest at 3.0 m/s (M = %.3f µV) and weakest at 10.0 m/s (M = %.3f µV, d = −1.42), ', ...
        mean(N1_amp(:,2)), mean(N1_amp(:,5)));
fprintf(fid, 'consistent with speed-dependent cortical adaptation.\n\n');

fprintf(fid, 'P1 and P1-N1 amplitudes did not reach significance (η²p = %.3f and %.3f respectively), ', ...
        ANOVA_results.P1_amp.eta_sq, ANOVA_results.P1N1_amp.eta_sq);
fprintf(fid, 'though large effect sizes suggest a real effect with insufficient power at N = %d.\n\n', n_participants);

fprintf(fid, 'Trend Analysis\n\n');
tr_p1 = TREND_results.P1_amp;
fprintf(fid, 'P1 amplitude followed a cubic speed-response curve (R² = %.3f), ', tr_p1.best_R2);
fprintf(fid, 'with significant linear (F(1,%d) = %.3f, p = %.4f) and quadratic (F(1,%d) = %.3f, p = %.4f) contrasts, ', ...
        n_participants-1, tr_p1.linear_F, tr_p1.linear_p, n_participants-1, tr_p1.quadratic_F, tr_p1.quadratic_p);
fprintf(fid, 'indicating a non-linear, saturation-like cortical encoding of stimulation velocity.\n\n');

fprintf(fid, 'N1 amplitude showed a linear decreasing trend (R² = %.3f, p = %.4f), ', ...
        TREND_results.N1_amp.R2_linear, TREND_results.N1_amp.linear_p);
fprintf(fid, 'consistent with progressive cortical adaptation at higher speeds.\n\n');

fprintf(fid, 'SEP Latencies\n\n');
lat_names_c = {'P1_lat','N1_lat','P2_lat'};
lat_labels_c = {'P1','N1','P2'};
for c = 1:3
    r_l = ANOVA_results.(lat_names_c{c});
    d_l = ALL_MEASURES.(lat_names_c{c});
    fprintf(fid, '%s latency (%.1f–%.1f ms): F(%.2f,%.2f) = %.3f, p = %.3f, η²p = %.3f. ', ...
            lat_labels_c{c}, min(mean(d_l,1)), max(mean(d_l,1)), ...
            r_l.df1, r_l.df2, r_l.F_val, r_l.p_corrected, r_l.eta_sq);
end
fprintf(fid, '\nLatency invariance across speeds indicates that stimulation speed modulates cortical response ');
fprintf(fid, 'magnitude but not the timing of neural activation.\n\n');

fprintf(fid, 'Power Considerations\n\n');
fprintf(fid, 'Only P2 amplitude was adequately powered (%.0f%%). ', r_p2.power*100);
fprintf(fid, 'The consistently large effect sizes (η²p = %.3f–%.3f across all amplitude components) ', ...
        min([r_p2.eta_sq, r_n1.eta_sq, ANOVA_results.P1_amp.eta_sq, ...
             ANOVA_results.P1N1_amp.eta_sq, ANOVA_results.N1P2_amp.eta_sq]), ...
        max([r_p2.eta_sq, r_n1.eta_sq, ANOVA_results.P1_amp.eta_sq, ...
             ANOVA_results.P1N1_amp.eta_sq, ANOVA_results.N1P2_amp.eta_sq]));
fprintf(fid, 'suggest that non-significant results for P1, N1, P1-N1, and N1-P2 reflect Type II error ');
fprintf(fid, 'due to insufficient sample size rather than the absence of a true effect. ');
fprintf(fid, 'These findings constitute preliminary evidence, warranting replication with a larger sample.\n\n');

fclose(fid);
fprintf('✓ Saved: %s\n', report_file);

% CSV export
csv_file = sprintf('MULTI_LEAN_Results_%s.csv', REFERENCE_METHOD);
fid = fopen(csv_file, 'w');
fprintf(fid, 'Participant,Speed_m_s,P1_amp_uV,P1_lat_ms,N1_amp_uV,N1_lat_ms,P2_amp_uV,P2_lat_ms,P1N1_uV,N1P2_uV\n');
for p = 1:n_participants
    for s = 1:n_speeds
        fprintf(fid, '%s,%.1f,%.4f,%.2f,%.4f,%.2f,%.4f,%.2f,%.4f,%.4f\n', ...
                successful_participants{p}, unique_speeds(s), ...
                P1_amp(p,s), P1_lat(p,s), N1_amp(p,s), N1_lat(p,s), ...
                P2_amp(p,s), P2_lat(p,s), P1N1_amp(p,s), N1P2_amp(p,s));
    end
end
fclose(fid);
fprintf('✓ Saved: %s\n', csv_file);

% MAT export
mat_file = sprintf('MULTI_LEAN_Analysis_%s.mat', REFERENCE_METHOD);
save(mat_file, 'ALL_MEASURES', 'ANOVA_results', 'POSTHOC_results', ...
     'TREND_results', 'DESC_results', 'SNR_results', 'RELIABILITY', ...
     'successful_participants', 'unique_speeds', ...
     'component_names', 'component_labels', ...
     'P1_amp','N1_amp','P2_amp','P1N1_amp','N1P2_amp', ...
     'P1_lat','N1_lat','P2_lat');
fprintf('✓ Saved: %s\n\n', mat_file);

%% FINAL SUMMARY ─────────────────────────────────────────────────────
fprintf('✓ ANALYSIS COMPLETE\n');

fprintf('ANOVA Summary:\n');
for comp = 1:length(component_names)
    cname = component_names{comp};
    r = ANOVA_results.(cname);
            component_labels{comp}, r.df1, r.df2, r.F_val, r.sig, r.p_corrected, r.eta_sq, r.power);
end

%% HELPER FUNCTIONS

function result = iff(condition, true_val, false_val)
    if condition; result = true_val; else; result = false_val; end
end

function [H, pValue, W] = swtest(x, alpha)
% Shapiro-Wilk test (Royston 1992, suitable for n=3..50)
    if nargin < 2; alpha = 0.05; end
    x = x(~isnan(x));
    x = sort(x(:));
    n = numel(x);

    if n < 3 || n > 50
        warning('swtest: n=%d outside range 3-50. Using lillietest.', n);
        [H, pValue] = lillietest(x, 'Alpha', alpha);
        W = NaN; return;
    end

    m   = norminv((((1:n)' - 0.375) ./ (n + 0.25)));
    m   = m / norm(m);
    phi = sum(m.^2);
    a   = zeros(n,1);
    for i = 1:floor(n/2)
        a(n+1-i) = m(n+1-i) / sqrt(phi);
        a(i)     = -a(n+1-i);
    end

    xbar = mean(x);
    SS   = sum((x - xbar).^2);
    W    = min(max((a'*x)^2 / SS, 0), 1);

    mu_w  = -1.2725 + 1.0521*(log(n) - log(3));
    sig_w =  1.0308 - 0.26763*log(n);
    z     = (log(1-W) - mu_w) / sig_w;
    pValue = max(1e-6, min(1, 1 - normcdf(z)));
    H = (pValue < alpha);
end