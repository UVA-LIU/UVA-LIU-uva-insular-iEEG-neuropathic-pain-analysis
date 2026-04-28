function Results = run_connectivity_lme_models(T, rootPath)
% =========================================================================
% run_connectivity_lme_all
%
% Combined connectivity mixed-effects analysis for the manuscript.
%
% INPUT
%   T         : connectivity epoch-level long table from
%               build_conn_epoch_longtable.m
%   rootPath  : root data folder
%
% ANALYSES
%   1) Trait architecture
%      a. Posterior vs Anterior gradient
%      b. Contralateral vs Ipsilateral laterality
%
%   2) State modulation
%      a. Pain intensity
%      b. Unpleasantness
%
% PIPELINE
%   - aggregates epoch-level connectivity to pair-day means
%   - applies QC threshold: minimum 10 epochs per pair-day
%   - transforms connectivity:
%       COH -> Fisher z
%       PLV -> logit
%   - centers DayIndex within subject
%   - centers PainIntensity and Unpleasantness within subject
%   - fits LME models with:
%       (1|Subject) + (1|Subject:PairID)
%
% OUTPUT
%   Results.TraitGradient
%   Results.TraitLaterality
%   Results.State
%
% SAVES
%   <rootPath>/MixedEffects/CONN_Results/
% =========================================================================

if nargin < 2 || isempty(rootPath)
    rootPath = 'E:\DBS_Data';
end

if nargin < 1 || isempty(T)
    matFile = fullfile(rootPath, 'MixedEffects', 'CONN_Epoch_Level_LongTable.mat');
    csvFile = fullfile(rootPath, 'MixedEffects', 'CONN_Epoch_Level_LongTable.csv');

    if exist(matFile, 'file') == 2
        fprintf('Loading T from MAT: %s\n', matFile);
        S = load(matFile, 'T');
        T = S.T;
    elseif exist(csvFile, 'file') == 2
        fprintf('Loading T from CSV: %s\n', csvFile);
        T = readtable(csvFile);
    else
        error('Connectivity long table not found.');
    end
end

assert(~isempty(T), 'Input table T is empty.');

reqVars = {'PairID','Subject','PairType','Metric','Band','Value','DayIndex',...
           'Pain_Intensity','Unpleasantness'};
assert(all(ismember(reqVars, T.Properties.VariableNames)), ...
    'Input table missing required variables.');

%% ------------------------------------------------------------------------
% PREP
% -------------------------------------------------------------------------
fprintf('Preparing connectivity table...\n');

T.Subject        = categorical(string(T.Subject));
T.PairID         = categorical(string(T.PairID));
T.PairType       = upper(strtrim(string(T.PairType)));
T.Metric         = categorical(strtrim(lower(string(T.Metric))));
T.Band           = categorical(strtrim(lower(string(T.Band))));
T.DayIndex       = double(T.DayIndex);
T.Value          = double(T.Value);
T.Pain_Intensity = double(T.Pain_Intensity);
T.Unpleasantness = double(T.Unpleasantness);

% Keep only valid rows
keep = ~isnan(T.Value);
T = T(keep, :);

if isempty(T)
    error('No valid connectivity rows remain after filtering.');
end

%% ------------------------------------------------------------------------
% OUTPUT SETUP
% -------------------------------------------------------------------------
ts = datestr(now, 'yyyymmdd_HHMMSS');
outDir = fullfile(rootPath, 'MixedEffects', 'CONN_Results');
if ~exist(outDir, 'dir')
    mkdir(outDir);
end

logFile = fullfile(outDir, ['LOG_CONN_ALL_' ts '.txt']);
fid = fopen(logFile, 'w');
dualfprintf(fid, '=== CONNECTIVITY MIXED-EFFECTS ANALYSIS ===\n');
dualfprintf(fid, 'Trait + State combined script\n');
dualfprintf(fid, 'QC: pair-day mean, minimum 10 epochs\n\n');

traitCols = {'Metric','Band','Region','Contrast','Estimate','SE','tStat','DF',...
             'pValue','PartialR2','N_PairDays','Subjects','FDR'};
traitTypes = {'string','string','string','string','double','double','double',...
              'double','double','double','double','double','double'};

stateCols = {'Predictor','Metric','Band','Estimate','SE','tStat','DF',...
             'pValue','PartialR2','N_PairDays','Subjects','FDR'};
stateTypes = {'string','string','string','double','double','double',...
              'double','double','double','double','double','double'};

TraitGradient = table('Size',[0,numel(traitCols)], ...
    'VariableNames', traitCols, 'VariableTypes', traitTypes);

TraitLaterality = table('Size',[0,numel(traitCols)], ...
    'VariableNames', traitCols, 'VariableTypes', traitTypes);

StateResults = table('Size',[0,numel(stateCols)], ...
    'VariableNames', stateCols, 'VariableTypes', stateTypes);

jobs = [ ...
    struct('metric',"plv", 'bands',["delta","theta","alpha","beta","gamma"]); ...
    struct('metric',"coh", 'bands',["delta","theta","alpha","beta","gamma"]) ...
    ];

%% ------------------------------------------------------------------------
% MAIN LOOP
% -------------------------------------------------------------------------
for j = 1:numel(jobs)
    metricStr = jobs(j).metric;
    bandsStr  = jobs(j).bands;

    for b = 1:numel(bandsStr)
        bandStr = bandsStr(b);
        metricCat = categorical(strtrim(lower(string(metricStr))));
        bandCat   = categorical(strtrim(lower(string(bandStr))));

        ds_epoch = T(T.Metric == metricCat & T.Band == bandCat, :);
        if isempty(ds_epoch)
            continue;
        end

        % Aggregate epoch-level to pair-day means
        G = findgroups(ds_epoch.Subject, ds_epoch.DayIndex, ds_epoch.PairID);
        nEpochs = splitapply(@numel, ds_epoch.Value, G);
        MeanVal = splitapply(@nanmean, ds_epoch.Value, G);
        Subject = splitapply(@(x) x(1), ds_epoch.Subject, G);
        PairID  = splitapply(@(x) x(1), ds_epoch.PairID, G);
        PairType = splitapply(@(x) x(1), ds_epoch.PairType, G);
        DayIdx  = splitapply(@(x) x(1), ds_epoch.DayIndex, G);
        Pain    = splitapply(@nanmean, ds_epoch.Pain_Intensity, G);
        Unpl    = splitapply(@nanmean, ds_epoch.Unpleasantness, G);

        qcMask = nEpochs >= 10;
        if sum(qcMask) < 10
            continue;
        end

        % Connectivity transform
        raw = MeanVal(qcMask);
        epsv = 1e-6;
        raw = min(max(raw, epsv), 1-epsv);

        if metricStr == "coh"
            ConnT = atanh(raw);
        else
            ConnT = log(raw ./ (1 - raw));
        end

        ds_agg = table( ...
            Subject(qcMask), PairID(qcMask), PairType(qcMask), DayIdx(qcMask), ...
            Pain(qcMask), Unpl(qcMask), ConnT, ...
            'VariableNames', {'Subject','PairID','PairType','DayIndex',...
                              'PainIntensity','Unpleasantness','ConnT'});

        % Subject-centered covariates
        Gsub = findgroups(ds_agg.Subject);

        MeanDay = splitapply(@mean, ds_agg.DayIndex, Gsub);
        ds_agg.DayIndex_c = ds_agg.DayIndex - MeanDay(Gsub);

        MeanPain = splitapply(@mean, ds_agg.PainIntensity, Gsub);
        ds_agg.PainIntensity_c = ds_agg.PainIntensity - MeanPain(Gsub);

        MeanUnpl = splitapply(@mean, ds_agg.Unpleasantness, Gsub);
        ds_agg.Unpleasantness_c = ds_agg.Unpleasantness - MeanUnpl(Gsub);

        %% ----------------------------------------------------------------
        % TRAIT A: Posterior vs Anterior
        % -----------------------------------------------------------------
        maskGrad = contains(ds_agg.PairType, "WI_");
        ds_grad = ds_agg(maskGrad, :);

        if height(ds_grad) >= 10
            ds_grad.Axis = strings(height(ds_grad),1);
            ds_grad.Axis(contains(ds_grad.PairType, "_A_")) = "Anterior";
            ds_grad.Axis(contains(ds_grad.PairType, "_P_")) = "Posterior";
            ds_grad.Axis = categorical(ds_grad.Axis, {'Anterior','Posterior'});

            if numel(unique(ds_grad.Axis)) == 2
                try
                    lme = fitlme(ds_grad, ...
                        'ConnT ~ Axis + DayIndex_c + (1|Subject) + (1|Subject:PairID)');
                    coef = lme.Coefficients;
                    idx = find(strcmp(coef.Name, 'Axis_Posterior'), 1);

                    if ~isempty(idx)
                        row = extract_trait_row(coef, idx, metricStr, bandStr, ...
                            "Global", "Posterior>Anterior", ...
                            height(ds_grad), numel(unique(ds_grad.Subject)));
                        TraitGradient = [TraitGradient; row]; %#ok<AGROW>

                        dualfprintf(fid, '[%s %s] Trait gradient: β=%.3f, p=%.3e\n', ...
                            metricStr, bandStr, coef.Estimate(idx), coef.pValue(idx));
                    end
                catch ME
                    dualfprintf(fid, '[%s %s] Trait gradient failed: %s\n', ...
                        metricStr, bandStr, ME.message);
                end
            end
        end

        %% ----------------------------------------------------------------
        % TRAIT B: Contralateral vs Ipsilateral
        % -----------------------------------------------------------------
        families = {'WI_A','WI_P','AP'};
        for f = 1:numel(families)
            fam = families{f};

            maskFam = contains(ds_agg.PairType, fam) & ~contains(ds_agg.PairType, "XH");
            ds_fam = ds_agg(maskFam, :);

            if height(ds_fam) < 10
                continue;
            end

            ds_fam.Lat = strings(height(ds_fam),1);
            ds_fam.Lat(contains(ds_fam.PairType, "IPSI"))   = "Ipsilateral";
            ds_fam.Lat(contains(ds_fam.PairType, "CONTRA")) = "Contralateral";
            ds_fam.Lat = categorical(ds_fam.Lat, {'Ipsilateral','Contralateral'});

            if numel(unique(ds_fam.Lat)) < 2
                continue;
            end

            try
                lme = fitlme(ds_fam, ...
                    'ConnT ~ Lat + DayIndex_c + (1|Subject) + (1|Subject:PairID)');
                coef = lme.Coefficients;
                idx = find(strcmp(coef.Name, 'Lat_Contralateral'), 1);

                if ~isempty(idx)
                    row = extract_trait_row(coef, idx, metricStr, bandStr, ...
                        string(fam), "Contra>Ipsi", ...
                        height(ds_fam), numel(unique(ds_fam.Subject)));
                    TraitLaterality = [TraitLaterality; row]; %#ok<AGROW>
                end
            catch ME
                dualfprintf(fid, '[%s %s %s] Trait laterality failed: %s\n', ...
                    metricStr, bandStr, fam, ME.message);
            end
        end

        %% ----------------------------------------------------------------
        % STATE: Pain / Unpleasantness
        % -----------------------------------------------------------------
        predictors = {'PainIntensity_c','Unpleasantness_c'};
        predictorLabels = {'PainIntensity','Unpleasantness'};

        for p = 1:numel(predictors)
            predVar = predictors{p};
            predLabel = predictorLabels{p};

            try
                formula = sprintf('ConnT ~ %s + DayIndex_c + (1|Subject) + (1|Subject:PairID)', predVar);
                lme = fitlme(ds_agg, formula);
                coef = lme.Coefficients;
                idx = find(strcmp(coef.Name, predVar), 1);

                if ~isempty(idx)
                    row = extract_state_row(coef, idx, predLabel, metricStr, bandStr, ...
                        height(ds_agg), numel(unique(ds_agg.Subject)));
                    StateResults = [StateResults; row]; %#ok<AGROW>

                    if coef.pValue(idx) < 0.05
                        dualfprintf(fid, '[%s %s] %s: β=%.3f, p=%.3e\n', ...
                            metricStr, bandStr, predLabel, coef.Estimate(idx), coef.pValue(idx));
                    end
                end
            catch ME
                dualfprintf(fid, '[%s %s] State %s failed: %s\n', ...
                    metricStr, bandStr, predLabel, ME.message);
            end
        end
    end
end

%% ------------------------------------------------------------------------
% FDR
% -------------------------------------------------------------------------
if ~isempty(TraitGradient)
    TraitGradient.FDR = calc_fdr_by_metric(TraitGradient);
end
if ~isempty(TraitLaterality)
    TraitLaterality.FDR = calc_fdr_by_metric(TraitLaterality);
end
if ~isempty(StateResults)
    StateResults.FDR = calc_fdr_by_predictor_metric(StateResults);
end

%% ------------------------------------------------------------------------
% SAVE
% -------------------------------------------------------------------------
writetable(TraitGradient,   fullfile(outDir, ['CONN_Trait_Gradient_' ts '.csv']));
writetable(TraitLaterality, fullfile(outDir, ['CONN_Trait_Laterality_' ts '.csv']));
writetable(StateResults,    fullfile(outDir, ['CONN_State_' ts '.csv']));

save(fullfile(outDir, ['CONN_All_Results_' ts '.mat']), ...
    'TraitGradient', 'TraitLaterality', 'StateResults', '-v7.3');

dualfprintf(fid, '\nSaved results to: %s\n', outDir);
fclose(fid);

fprintf('✅ Connectivity mixed-effects analysis complete.\n');

Results = struct();
Results.TraitGradient   = TraitGradient;
Results.TraitLaterality = TraitLaterality;
Results.State           = StateResults;
end

%% =========================================================================
% HELPERS
% =========================================================================
function row = extract_trait_row(coef, idx, metricStr, bandStr, regionStr, contrastStr, N, nSub)
est = coef.Estimate(idx);
se  = coef.SE(idx);
t   = coef.tStat(idx);
df  = coef.DF(idx);
p   = coef.pValue(idx);
r2  = t^2 / (t^2 + df);

row = {string(metricStr), string(bandStr), string(regionStr), string(contrastStr), ...
       est, se, t, df, p, r2, N, nSub, nan};
end

function row = extract_state_row(coef, idx, predictorStr, metricStr, bandStr, N, nSub)
est = coef.Estimate(idx);
se  = coef.SE(idx);
t   = coef.tStat(idx);
df  = coef.DF(idx);
p   = coef.pValue(idx);
r2  = t^2 / (t^2 + df);

row = {string(predictorStr), string(metricStr), string(bandStr), ...
       est, se, t, df, p, r2, N, nSub, nan};
end

function fdr_col = calc_fdr_by_metric(Tin)
if isempty(Tin)
    fdr_col = [];
    return;
end
fdr_col = ones(height(Tin),1);
metrics = unique(Tin.Metric);
for i = 1:numel(metrics)
    idx = Tin.Metric == metrics(i);
    fdr_col(idx) = fdr_bh_local(Tin.pValue(idx));
end
end

function fdr_col = calc_fdr_by_predictor_metric(Tin)
if isempty(Tin)
    fdr_col = [];
    return;
end
fdr_col = ones(height(Tin),1);
preds = unique(Tin.Predictor);
metrics = unique(Tin.Metric);

for p = 1:numel(preds)
    for m = 1:numel(metrics)
        idx = Tin.Predictor == preds(p) & Tin.Metric == metrics(m);
        if any(idx)
            fdr_col(idx) = fdr_bh_local(Tin.pValue(idx));
        end
    end
end
end

function adj_p = fdr_bh_local(pvals)
pvals = pvals(:);
m = length(pvals);
if m == 0
    adj_p = [];
    return;
end

[p_sorted, sort_ids] = sort(pvals);
adj_sorted = zeros(m,1);
for i = 1:m
    adj_sorted(i) = p_sorted(i) * m / i;
end
for i = m-1:-1:1
    adj_sorted(i) = min(adj_sorted(i), adj_sorted(i+1));
end
adj_sorted(adj_sorted > 1) = 1;

adj_p = zeros(m,1);
adj_p(sort_ids) = adj_sorted;
end

function dualfprintf(fid, varargin)
fprintf(varargin{:});
if fid > 0
    fprintf(fid, varargin{:});
end
end