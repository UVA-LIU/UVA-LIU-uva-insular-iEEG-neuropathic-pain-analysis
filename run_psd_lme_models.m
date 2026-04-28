function Results = run_psd_lme_models(T, rootPath)
% =========================================================================
% run_psd_lme_quadrant_ipsi_contra_parallel
%
% PSD mixed-effects analysis with:
%   • Stage 1: Trait-like quadrant effects
%   • Stage 2A: Pain-only state models (Subject-Centered)
%   • Stage 2B: Unpleasantness-only state models (Subject-Centered)
%   • Stage 2C: Temporal Drift (Day_Index) extracted from Stage 2A models
%
% UPDATES:
%   1. Added SLOW INDEX support (Metric='log10', Band='slowindex').
%   2. Calculates 'Temporal Drift' (Day_Index coefficient).
%   3. Adds FDR q-values to all tables.
% =========================================================================

%% ---------- OPTION B: AUTO-LOAD TABLE ----------
if nargin < 1 || isempty(T)
    if nargin < 2 || isempty(rootPath)
        rootPath = 'E:\DBS_Data';
    end
    fprintf('▶ Loading PSD_Epoch_Level_LongTable.csv (Option B)\n');
    T = readtable(fullfile(rootPath,'MixedEffects','PSD_Epoch_Level_LongTable.csv'));
end

if nargin < 2 || isempty(rootPath)
    rootPath = 'E:\DBS_Data';
end

%% ---------- VERIFY INPUT ----------
reqVars = {'Subject','Region','PainSide','Band','Metric','Value','DayIndex', ...
           'Pain_Intensity','Unpleasantness'};
assert(all(ismember(reqVars, T.Properties.VariableNames)), ...
    'Missing required variables in input table.');

%% ---------- FORMAT ----------
T.Subject        = categorical(string(T.Subject));
T.Region         = categorical(string(T.Region));
T.PainSide       = categorical(string(T.PainSide));
% Convert Band/Metric to lowercase string for consistent matching
T.Band           = lower(string(T.Band));
T.Metric         = lower(string(T.Metric));
T.Value          = double(T.Value);
T.DayIndex       = double(T.DayIndex);
T.Pain_Intensity = double(T.Pain_Intensity);
T.Unpleasantness = double(T.Unpleasantness);

%% ---------- DERIVE LATERALITY & AXIS ----------
lat = strings(height(T),1);
ax  = strings(height(T),1);

% Identify Ipsilateral vs Contralateral based on PainSide
subjects = unique(T.Subject);
for i = 1:numel(subjects)
    s = subjects(i);
    mask = T.Subject == s;
    ps   = unique(T.PainSide(mask));
    if isempty(ps), continue; end
    ps = string(ps(1));
    
    if strcmpi(ps, "Left")
        lat(mask & contains(string(T.Region),"Left"))  = "ipsilateral";
        lat(mask & contains(string(T.Region),"Right")) = "contralateral";
    elseif strcmpi(ps, "Right")
        lat(mask & contains(string(T.Region),"Right")) = "ipsilateral";
        lat(mask & contains(string(T.Region),"Left"))  = "contralateral";
    end
end

ax(contains(string(T.Region),"Anterior"))  = "Anterior";
ax(contains(string(T.Region),"Posterior")) = "Posterior";

T.Laterality = categorical(lat, {'ipsilateral','contralateral'});
T.Axis       = categorical(ax,  {'Anterior','Posterior'});

% Create QuadRegion (e.g., 'ipsilateral_Anterior')
T.QuadRegion = categorical( ...
    strcat(string(T.Laterality),'_',string(T.Axis)), ...
    {'ipsilateral_Anterior','ipsilateral_Posterior', ...
     'contralateral_Anterior','contralateral_Posterior'} );

%% ---------- CENTER DAY (GLOBAL CONTROL) ----------
T.DayIndex_c = T.DayIndex - mean(T.DayIndex,'omitnan');

%% ---------- MODEL DEFINITIONS ----------
% UPDATED: Added 'log10' metric for 'slowindex' band
jobs = struct( ...
    'metric', {'db',     'rel_pct', 'pdf_hz', 'pdp_db', 'log10'}, ...
    'bands',  {{'delta','theta','alpha','beta','gamma'}, ...
               {'delta','theta','alpha','beta','gamma'}, ...
               {'peak'}, ...
               {'peak'}, ...
               {'slowindex'}} ); 

formula_stage1 = 'Value ~ 1 + QuadRegion + DayIndex_c + (1|Subject)';
formula_pain   = 'Value ~ 1 + DayIndex_c + Pain_Intensity_c + (1|Subject)';
formula_unpl   = 'Value ~ 1 + DayIndex_c + Unpleasantness_c + (1|Subject)';

%% ---------- OUTPUT ----------
ts = datestr(now,'yyyymmdd_HHMMSS');
outDir = fullfile(rootPath,'MixedEffects','PSD_Results');
if ~exist(outDir,'dir'), mkdir(outDir); end

Stage1       = table();
Stage2_PAIN  = table();
Stage2_UNPL  = table();
Stage2_DRIFT = table(); 

%% ============================================================
%% STAGE 1 — TRAIT-LIKE QUADRANT EFFECTS
%% ============================================================
fprintf('Running Stage 1 (Trait)...\n');
for j = 1:numel(jobs)
    for b = 1:numel(jobs(j).bands)
        % Filter Data
        mask = (T.Metric == jobs(j).metric) & (T.Band == jobs(j).bands{b});
        ds = T(mask, :);
        
        if isempty(ds)
            % Optional: Warn if missing (useful for debugging Slow Index)
            % fprintf('  Skipping %s %s (No data)\n', jobs(j).metric, jobs(j).bands{b});
            continue; 
        end
        
        try
            lme = fitlme(ds, formula_stage1, 'FitMethod','REML');
            coef = lme.Coefficients;
            for i = 1:height(coef)
                cname = string(coef.Name{i});
                if cname=="(Intercept)" || contains(cname,"DayIndex"), continue; end
                
                t  = coef.tStat(i);
                df = coef.DF(i);
                r2 = t^2/(t^2+df);
                
                Stage1 = [Stage1; {
                    jobs(j).metric, jobs(j).bands{b}, cname, ...
                    coef.Estimate(i), coef.SE(i), coef.pValue(i), r2, ...
                    height(ds), numel(unique(ds.Subject))
                }];
            end
        catch
        end
    end
end

if ~isempty(Stage1)
    Stage1.Properties.VariableNames = ...
     {'Metric','Band','Contrast','Estimate','SE','pValue','PartialR2','N','Subjects'};
end

%% ============================================================
%% STAGE 2A — PAIN & TEMPORAL DRIFT (WITHIN-SUBJECT CENTERED)
%% ============================================================
fprintf('Running Stage 2A (Pain & Drift)...\n');
for j = 1:numel(jobs)
    for b = 1:numel(jobs(j).bands)
        for q = categories(T.QuadRegion)'
            % Filter Data
            mask = (T.Metric == jobs(j).metric) & ...
                   (T.Band == jobs(j).bands{b}) & ...
                   (T.QuadRegion == string(q));
            ds = T(mask, :);
            
            if height(ds)<40, continue; end
            
            % Center Pain within Subject
            G = findgroups(ds.Subject);
            subjMean = splitapply(@mean, ds.Pain_Intensity, G);
            ds.Pain_Intensity_c = ds.Pain_Intensity - subjMean(G);
            
            try
                lme = fitlme(ds, formula_pain, 'FitMethod','REML');
                c   = lme.Coefficients;
                
                % --- 1. Extract PAIN Effect ---
                idx = strcmp(c.Name,'Pain_Intensity_c');
                if any(idx)
                    t  = c.tStat(idx);
                    df = c.DF(idx);
                    r2 = t^2/(t^2+df);
                    
                    Stage2_PAIN = [Stage2_PAIN; {
                        jobs(j).metric, jobs(j).bands{b}, char(q), 'Pain_Intensity_c', ...
                        c.Estimate(idx), c.SE(idx), c.pValue(idx), r2, height(ds)
                    }];
                end
                
                % --- 2. Extract DRIFT Effect (DayIndex_c) ---
                idxD = strcmp(c.Name,'DayIndex_c');
                if any(idxD)
                    t  = c.tStat(idxD);
                    df = c.DF(idxD);
                    r2 = t^2/(t^2+df);
                    
                    Stage2_DRIFT = [Stage2_DRIFT; {
                        jobs(j).metric, jobs(j).bands{b}, char(q), 'Day_Index', ...
                        c.Estimate(idxD), c.SE(idxD), c.pValue(idxD), r2, height(ds)
                    }];
                end
            catch
            end
        end
    end
end

if ~isempty(Stage2_PAIN)
    Stage2_PAIN.Properties.VariableNames = ...
     {'Metric','Band','Quadrant','Predictor','Estimate','SE','pValue','PartialR2','N'};
end
if ~isempty(Stage2_DRIFT)
    Stage2_DRIFT.Properties.VariableNames = ...
     {'Metric','Band','Quadrant','Predictor','Estimate','SE','pValue','PartialR2','N'};
end

%% ============================================================
%% STAGE 2B — UNPLEASANTNESS-ONLY (WITHIN-SUBJECT CENTERED)
%% ============================================================
fprintf('Running Stage 2B (Unpleasantness)...\n');
for j = 1:numel(jobs)
    for b = 1:numel(jobs(j).bands)
        for q = categories(T.QuadRegion)'
            % Filter Data
            mask = (T.Metric == jobs(j).metric) & ...
                   (T.Band == jobs(j).bands{b}) & ...
                   (T.QuadRegion == string(q));
            ds = T(mask, :);
            
            if height(ds)<40, continue; end
            
            % Center Unpleasantness
            G = findgroups(ds.Subject);
            subjMean = splitapply(@mean, ds.Unpleasantness, G);
            ds.Unpleasantness_c = ds.Unpleasantness - subjMean(G);
            
            try
                lme = fitlme(ds, formula_unpl, 'FitMethod','REML');
                c   = lme.Coefficients;
                idx = strcmp(c.Name,'Unpleasantness_c');
                
                if any(idx)
                    t  = c.tStat(idx);
                    df = c.DF(idx);
                    r2 = t^2/(t^2+df);
                    
                    Stage2_UNPL = [Stage2_UNPL; {
                        jobs(j).metric, jobs(j).bands{b}, char(q), 'Unpleasantness_c', ...
                        c.Estimate(idx), c.SE(idx), c.pValue(idx), r2, height(ds)
                    }];
                end
            catch
            end
        end
    end
end

if ~isempty(Stage2_UNPL)
    Stage2_UNPL.Properties.VariableNames = ...
     {'Metric','Band','Quadrant','Predictor','Estimate','SE','pValue','PartialR2','N'};
end

%% ============================================================
%% APPLY FDR CORRECTION (BENJAMINI-HOCHBERG)
%% ============================================================
fprintf('Applying FDR Correction (q-value)...\n');
if ~isempty(Stage1),       Stage1.qValue = fdr_bh(Stage1.pValue); end
if ~isempty(Stage2_PAIN),  Stage2_PAIN.qValue = fdr_bh(Stage2_PAIN.pValue); end
if ~isempty(Stage2_UNPL),  Stage2_UNPL.qValue = fdr_bh(Stage2_UNPL.pValue); end
if ~isempty(Stage2_DRIFT), Stage2_DRIFT.qValue = fdr_bh(Stage2_DRIFT.pValue); end

%% ---------- SAVE RESULTS (MAT + CSV) ----------
resultsMat = fullfile(outDir, ['LME_PSD_QUADRANT_PARALLEL_' ts '.mat']);
save(resultsMat, 'Stage1', 'Stage2_PAIN', 'Stage2_UNPL', 'Stage2_DRIFT', '-v7.3');

writetable(Stage1,       fullfile(outDir, ['LME_PSD_STAGE1_TRAIT_' ts '.csv']));
writetable(Stage2_PAIN,  fullfile(outDir, ['LME_PSD_STAGE2_PAIN_' ts '.csv']));
writetable(Stage2_UNPL,  fullfile(outDir, ['LME_PSD_STAGE2_UNPLEASANTNESS_' ts '.csv']));
writetable(Stage2_DRIFT, fullfile(outDir, ['LME_PSD_STAGE2_DRIFT_' ts '.csv']));

fprintf('✅ Saved PSD mixed-effects results with FDR q-values:\n');
fprintf('   MAT: %s\n', resultsMat);
fprintf('   CSV: Stage1 / Stage2_PAIN / Stage2_UNPL / Stage2_DRIFT\n');

end

%% ============================================================
%% HELPER FUNCTION: FDR (Benjamini-Hochberg)
%% ============================================================
function [q] = fdr_bh(p)
    m = length(p);
    [sorted_p, sort_idx] = sort(p);
    q = zeros(m,1);
    if m > 0
        q(m) = sorted_p(m);
        for i = m-1:-1:1
            q(i) = min(q(i+1), sorted_p(i) * (m/i));
        end
    end
    q(q > 1) = 1;
    original_q = zeros(m,1);
    original_q(sort_idx) = q;
    q = original_q;
end