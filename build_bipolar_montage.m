function status = build_bipolar_montage(rootPath, subjectIDs, regionDef, opts)
% build_bipolar_montage_from_raw
%
% Build region-constrained bipolar montages from <subj>_RawEEG.mat.
%
% Rule:
%   Bipolar pairs are created only between consecutive contacts within the
%   same anatomically defined region. No cross-region pairs are created.
%
% Required inputs
%   rootPath    : root data directory
%   subjectIDs  : cell array of subject IDs
%
% Optional inputs
%   regionDef   : region definition structure. May be:
%                 1) global region map used for all subjects, or
%                 2) struct with subject-level maps, e.g. regionDef.DBS001
%   opts        : optional structure:
%                 .writeLogs       logical, default true
%                 .writeSummary     logical, default true
%                 .outputSuffix     char, default '_BipolarEEG.mat'
%
% Expected input
%   <root>/<subj>/<subj>_setfiles/<subj>_RawEEG.mat
%
% Output
%   <root>/<subj>/<subj>_setfiles/<subj>_BipolarEEG.mat

if nargin < 1 || isempty(rootPath)
    error('rootPath is required.');
end
if nargin < 2 || isempty(subjectIDs)
    error('subjectIDs is required.');
end
if nargin < 3
    regionDef = [];
end
if nargin < 4 || isempty(opts)
    opts = struct();
end

if ~isfield(opts, 'writeLogs'),   opts.writeLogs = true; end
if ~isfield(opts, 'writeSummary'), opts.writeSummary = true; end
if ~isfield(opts, 'outputSuffix'), opts.outputSuffix = '_BipolarEEG.mat'; end

if ~exist(rootPath, 'dir')
    error('Root path not found: %s', rootPath);
end

status = struct('Subject',{}, 'Days',{}, 'TotalBipolar',{}, 'Warnings',{});

for s = 1:numel(subjectIDs)
    subj = char(subjectIDs{s});
    inDir  = fullfile(rootPath, subj, [subj '_setfiles']);
    inFile = fullfile(inDir, [subj '_RawEEG.mat']);

    if exist(inFile, 'file') ~= 2
        warning('Missing RawEEG file for %s: %s', subj, inFile);
        status(end+1) = struct('Subject', subj, 'Days', 0, ...
            'TotalBipolar', 0, 'Warnings', "missing RawEEG"); %#ok<AGROW>
        continue
    end

    L = load(inFile);
    if ~isfield(L, 'Mono') || isempty(L.Mono) || ~isfield(L, 'Meta')
        warning('RawEEG file for %s does not contain Mono/Meta.', subj);
        status(end+1) = struct('Subject', subj, 'Days', 0, ...
            'TotalBipolar', 0, 'Warnings', "missing Mono/Meta"); %#ok<AGROW>
        continue
    end

    Mono = L.Mono;
    Meta = L.Meta;
    PainSide = "";
    if isfield(L, 'PainSide')
        PainSide = string(L.PainSide);
    end

    subjRegionDef = get_subject_region_def(regionDef, subj);
    hasRegionMap = ~isempty(subjRegionDef);

    EEG_bipolar_list = cell(numel(Mono), 1);
    RegionByDay_idx  = cell(numel(Mono), 1);
    BipPairsByDay    = cell(numel(Mono), 1);
    BipLabelsByDay   = cell(numel(Mono), 1);

    fid = -1;
    if opts.writeLogs
        logPath = fullfile(inDir, [subj '_Bipolar_Build_Log.txt']);
        fid = fopen(logPath, 'w');
        if fid > 0
            fprintf(fid, 'Bipolar montage construction log\n');
            fprintf(fid, 'Subject: %s\n', subj);
            fprintf(fid, 'Strict within-region pairing: enabled\n');
            fprintf(fid, 'Region map provided: %d\n\n', hasRegionMap);
        end
    end

    totalBipolar = 0;

    for d = 1:numel(Mono)
        if fid > 0
            fprintf(fid, '\nDay index %d\n', d);
        end

        EEGd = Mono(d).EEG;
        if isempty(EEGd) || ~isfield(EEGd, 'data') || isempty(EEGd.data)
            if fid > 0, fprintf(fid, '  Empty data, skipped.\n'); end
            EEG_bipolar_list{d} = [];
            RegionByDay_idx{d}  = empty_region_index_struct();
            BipPairsByDay{d}    = [];
            BipLabelsByDay{d}   = {};
            continue
        end

        Xmono  = double(EEGd.data);
        nbchan = size(Xmono, 1);
        srate  = EEGd.srate;

        dayLabs = get_day_labels(Meta, d, nbchan);

        [regionIdx, regionNames, regionSides, regionAxes] = resolve_regions_for_day( ...
            Meta, d, dayLabs, subjRegionDef, hasRegionMap, fid);

        allPairs = [];
        allLabels = {};
        allRegions = {};
        allSides = {};
        allAxes = {};
        allLaterality = {};

        for r = 1:numel(regionNames)
            regName = regionNames{r};
            idx = regionIdx.(regName)(:).';

            if numel(idx) < 2
                if fid > 0
                    fprintf(fid, '  %-16s : fewer than 2 contacts, none built.\n', regName);
                end
                continue
            end

            pairs = [idx(1:end-1).' idx(2:end).'];

            for k = 1:size(pairs, 1)
                label1 = char(dayLabs(pairs(k,1)));
                label2 = char(dayLabs(pairs(k,2)));

                if isempty(label1), label1 = sprintf('Ch%d', pairs(k,1)); end
                if isempty(label2), label2 = sprintf('Ch%d', pairs(k,2)); end

                allPairs(end+1,:) = pairs(k,:); %#ok<AGROW>
                allLabels{end+1}  = sprintf('%s-%s', label1, label2); %#ok<AGROW>
                allRegions{end+1} = regName; %#ok<AGROW>
                allSides{end+1}   = regionSides{r}; %#ok<AGROW>
                allAxes{end+1}    = regionAxes{r}; %#ok<AGROW>
                allLaterality{end+1} = laterality_relative_to_pain(regionSides{r}, PainSide); %#ok<AGROW>
            end
        end

        if isempty(allPairs)
            if fid > 0, fprintf(fid, '  No bipolar pairs constructed.\n'); end
            EEG_bipolar_list{d} = [];
            RegionByDay_idx{d}  = regionIdx;
            BipPairsByDay{d}    = [];
            BipLabelsByDay{d}   = {};
            continue
        end

        nPairs = size(allPairs, 1);
        Xbip = zeros(nPairs, size(Xmono, 2), 'like', Xmono);

        for p = 1:nPairs
            Xbip(p,:) = Xmono(allPairs(p,2),:) - Xmono(allPairs(p,1),:);
        end

        EEG_bip = struct();
        EEG_bip.data       = Xbip;
        EEG_bip.srate      = srate;
        EEG_bip.chanlabels = allLabels(:);
        EEG_bip.pairs      = allPairs;
        EEG_bip.region     = allRegions(:);
        EEG_bip.side       = allSides(:);
        EEG_bip.axis       = allAxes(:);
        EEG_bip.laterality = allLaterality(:);

        EEG_bipolar_list{d} = EEG_bip;
        RegionByDay_idx{d}  = regionIdx;
        BipPairsByDay{d}    = allPairs;
        BipLabelsByDay{d}   = EEG_bip.chanlabels;

        totalBipolar = totalBipolar + nPairs;

        if fid > 0
            fprintf(fid, '  Built %d bipolar channels, sampling rate = %g Hz\n', nPairs, srate);
            for r = 1:numel(regionNames)
                regName = regionNames{r};
                mask = strcmp(allRegions, regName);
                fprintf(fid, '    %-16s : %d\n', regName, nnz(mask));
            end
        end
    end

    if fid > 0
        fclose(fid);
    end

    Meta.RegionByDay        = RegionByDay_idx;
    Meta.BipolarPairsByDay  = BipPairsByDay;
    Meta.BipolarLabelsByDay = BipLabelsByDay;

    outFile = fullfile(inDir, [subj opts.outputSuffix]);
    save(outFile, 'EEG_bipolar_list', 'Meta', 'PainSide', '-v7.3');

    if opts.writeSummary
        try
            DayIndex = (1:numel(Mono))';
            nBipolar = zeros(numel(Mono),1);
            for dd = 1:numel(Mono)
                if ~isempty(BipPairsByDay{dd})
                    nBipolar(dd) = size(BipPairsByDay{dd}, 1);
                end
            end
            Tsum = table(DayIndex, nBipolar, ...
                'VariableNames', {'DayIndex','nBipolar'});
            save(fullfile(inDir, [subj '_Bipolar_Summary.mat']), 'Tsum');
        catch ME
            warning('Failed to save bipolar summary for %s: %s', subj, ME.message);
        end
    end

    fprintf('Saved bipolar montage for %s: %s\n', subj, outFile);

    status(end+1) = struct('Subject', subj, 'Days', numel(Mono), ...
        'TotalBipolar', totalBipolar, 'Warnings', ""); %#ok<AGROW>
end
end

% =========================================================================
% helpers
% =========================================================================

function S = empty_region_index_struct()
S = struct('LeftAnterior', [], 'LeftPosterior', [], ...
           'RightAnterior', [], 'RightPosterior', []);
end

function subjRegionDef = get_subject_region_def(regionDef, subj)
subjRegionDef = [];

if isempty(regionDef)
    return
end

if isstruct(regionDef) && isfield(regionDef, subj)
    subjRegionDef = normalize_region_struct(regionDef.(subj));
else
    subjRegionDef = normalize_region_struct(regionDef);
end
end

function dayLabs = get_day_labels(Meta, d, nbchan)
if isfield(Meta, 'LabelsByDay') && numel(Meta.LabelsByDay) >= d && ~isempty(Meta.LabelsByDay{d})
    labels = string(Meta.LabelsByDay{d});
    labels = arrayfun(@(x) canonicalize_label(x), labels);

    if numel(labels) < nbchan
        tmp = labels;
        labels = strings(nbchan, 1);
        labels(1:numel(tmp)) = tmp;
    else
        labels = labels(1:nbchan);
    end
    dayLabs = labels(:);
else
    dayLabs = compose("Ch%d", 1:nbchan).';
end
end

function [regionIdx, regionNames, regionSides, regionAxes] = resolve_regions_for_day( ...
    Meta, d, dayLabs, subjRegionDef, hasRegionMap, fid)

regionNames = {'LeftAnterior','LeftPosterior','RightAnterior','RightPosterior'};
regionSides = {'Left','Left','Right','Right'};
regionAxes  = {'Anterior','Posterior','Anterior','Posterior'};
regionIdx = empty_region_index_struct();

if hasRegionMap && ~isempty(subjRegionDef)
    for j = 1:numel(regionNames)
        name = regionNames{j};
        labset = string(subjRegionDef.(name)(:));
        labset = arrayfun(@(x) canonicalize_label(x), labset);

        idxs = zeros(1, numel(labset));

        for k = 1:numel(labset)
            pos = find(dayLabs == labset(k), 1, 'first');
            if isempty(pos)
                if fid > 0
                    fprintf(fid, '  Region label not found today: %s, region: %s\n', labset(k), name);
                end
                idxs(k) = 0;
            else
                idxs(k) = pos;
            end
        end

        regionIdx.(name) = idxs(idxs >= 1 & idxs <= numel(dayLabs));
    end
    return
end

if isfield(Meta, 'RegionByDay') && numel(Meta.RegionByDay) >= d && ~isempty(Meta.RegionByDay{d})
    R = Meta.RegionByDay{d};

    if iscell(R) && numel(R) == numel(dayLabs)
        for i = 1:numel(dayLabs)
            regionLabel = string(R{i});
            isLeft = startsWith(dayLabs(i), "L");

            if regionLabel == "Anterior"
                if isLeft
                    regionIdx.LeftAnterior(end+1) = i;
                else
                    regionIdx.RightAnterior(end+1) = i;
                end
            elseif regionLabel == "Posterior"
                if isLeft
                    regionIdx.LeftPosterior(end+1) = i;
                else
                    regionIdx.RightPosterior(end+1) = i;
                end
            elseif any(regionLabel == string(regionNames))
                regionIdx.(char(regionLabel))(end+1) = i;
            end
        end
    elseif isstruct(R)
        R = normalize_region_struct(R);
        regionIdx = R;
    end
    return
end

if fid > 0
    fprintf(fid, '  No region information found for this day. No bipolar pairs built.\n');
end
end

function out = canonicalize_label(inLab)
t = upper(string(inLab));
t = regexprep(t, '[\s\-\_\(\)]', '');
t = regexprep(t, '^LI', 'L');
t = regexprep(t, '^RI', 'R');

if all(isstrprop(t,'digit'))
    n = str2double(t);
    if n >= 1 && n <= 9
        out = "L" + n; return
    elseif n >= 10 && n <= 18
        out = "R" + (n - 9); return
    end
end

mL = regexp(t, '^L(\d+)$', 'tokens', 'once');
if ~isempty(mL)
    k = str2double(mL{1});
    if k >= 1 && k <= 9
        out = "L" + k; return
    end
end

mR = regexp(t, '^R(\d+)$', 'tokens', 'once');
if ~isempty(mR)
    k = str2double(mR{1});
    if k >= 1 && k <= 9
        out = "R" + k; return
    elseif k >= 10 && k <= 18
        out = "R" + (k - 9); return
    end
end

out = t;
end

function out = normalize_region_struct(S)
out = struct('LeftAnterior', {{}}, 'LeftPosterior', {{}}, ...
             'RightAnterior', {{}}, 'RightPosterior', {{}});

if ~isstruct(S)
    return
end

fields = fieldnames(S);
for i = 1:numel(fields)
    key = lower(regexprep(fields{i}, '[^a-z0-9]', ''));
    values = S.(fields{i});

    if isstring(values), values = cellstr(values); end
    if ~iscell(values), values = cellstr(string(values)); end

    switch key
        case {'leftanterior','leftant','la','lant'}
            out.LeftAnterior = values;
        case {'leftposterior','leftpost','lp','lpost'}
            out.LeftPosterior = values;
        case {'rightanterior','rightant','ra','rant'}
            out.RightAnterior = values;
        case {'rightposterior','rightpost','rp','rpost'}
            out.RightPosterior = values;
    end
end
end

function lat = laterality_relative_to_pain(sideLabel, painSide)
sideLabel = lower(string(sideLabel));
painSide = lower(string(painSide));

if painSide == "" || sideLabel == ""
    lat = "";
elseif startsWith(sideLabel, "left") && any(painSide == ["l","left"])
    lat = "ipsi";
elseif startsWith(sideLabel, "right") && any(painSide == ["r","right"])
    lat = "ipsi";
else
    lat = "contra";
end
end