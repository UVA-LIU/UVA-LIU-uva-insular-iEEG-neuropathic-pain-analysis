function status = make_raweeg_from_setfiles(rootPath, subjectInfo, subjectRatings, regionDef, opts)
% load_raw_eeg_from_setfiles
%
% Load EEGLAB .set files, canonicalize monopolar labels, attach day-level
% ratings, derive left/right and ipsi/contra labels, and save per-subject
% RawEEG structures for downstream analysis.
%
% Required inputs
%   rootPath      : root data directory
%   subjectInfo   : struct with one field per subject. Each field contains a
%                   cell array of day tags in chronological order.
%
% Optional inputs
%   subjectRatings: struct or table-based ratings input
%   regionDef     : optional anatomical region definition
%   opts          : struct with optional fields:
%       .painSideMap          struct mapping subject IDs to 'L' or 'R'
%       .writeLogs            logical, default true
%       .writeSummaries       logical, default true
%       .writeMergedRatings   logical, default true
%       .logFileName          char, default 'eeg_load_log.txt'
%       .painSideCsvNames     cellstr, default {'Pain_Side.csv','pain_side.csv'}
%
% Supported root-level ratings files
%   Pain_intensity.csv
%   Unpleasantness.csv
%   Pain_Side.csv
%
% Accepted ratings formats
%   1. long format by day index: Subject, Day, Value
%   2. long format by date: Subject, Date, Value
%   3. wide format by date
%
% Outputs per subject
%   <subj>_RawEEG.mat
%   <subj>_Ratings_Merged.csv      if opts.writeMergedRatings = true
%   <subj>_RawEEG_Summary.txt      if opts.writeSummaries = true
%
% Notes
%   This public version does not contain cohort-specific defaults.
%   Subject metadata such as day tags and pain side should be provided
%   externally through subjectInfo and opts.painSideMap.

% -------------------------------------------------------------------------
% defaults
% -------------------------------------------------------------------------
if nargin < 1 || isempty(rootPath)
    error('rootPath is required.');
end
if nargin < 2 || isempty(subjectInfo)
    error('subjectInfo is required and must be a struct.');
end
if nargin < 3 || isempty(subjectRatings)
    subjectRatings = struct();
end
if nargin < 4
    regionDef = [];
end
if nargin < 5 || isempty(opts)
    opts = struct();
end

if ~isfield(opts, 'painSideMap'),        opts.painSideMap = struct(); end
if ~isfield(opts, 'writeLogs'),          opts.writeLogs = true; end
if ~isfield(opts, 'writeSummaries'),     opts.writeSummaries = true; end
if ~isfield(opts, 'writeMergedRatings'), opts.writeMergedRatings = true; end
if ~isfield(opts, 'logFileName'),        opts.logFileName = 'eeg_load_log.txt'; end
if ~isfield(opts, 'painSideCsvNames')
    opts.painSideCsvNames = {'Pain_Side.csv','pain_side.csv'};
end

% -------------------------------------------------------------------------
% preconditions
% -------------------------------------------------------------------------
if exist('pop_loadset','file') ~= 2
    error('EEGLAB function pop_loadset not found. Add EEGLAB to the MATLAB path.');
end
if ~exist(rootPath,'dir')
    error('Root path not found: %s', rootPath);
end

painSideMap = opts.painSideMap;

% optional pain side override from CSV
for i = 1:numel(opts.painSideCsvNames)
    sideCsv = fullfile(rootPath, opts.painSideCsvNames{i});
    if exist(sideCsv, 'file') == 2
        try
            Tside = readtable(sideCsv, 'TextType', 'string');
            v = lower(string(Tside.Properties.VariableNames));
            iS = find(v=="subject" | v=="subj" | v=="id", 1);
            iP = find(v=="side" | v=="pain_side" | v=="painside", 1);
            if ~isempty(iS) && ~isempty(iP)
                for r = 1:height(Tside)
                    sid = upper(strtrim(string(Tside{r,iS})));
                    if startsWith(sid, "DBS")
                        painSideMap.(char(sid)) = upper(strtrim(string(Tside{r,iP})));
                    end
                end
            end
        catch ME
            warning('Could not read pain side file %s: %s', sideCsv, ME.message);
        end
        break
    end
end

% -------------------------------------------------------------------------
% logging
% -------------------------------------------------------------------------
fid = -1;
if opts.writeLogs
    logFile = fullfile(rootPath, opts.logFileName);
    [fid,msg] = fopen(logFile,'w');
    if fid < 0
        error('Cannot open log file: %s', msg);
    end
    fprintf(fid, 'EEG load log\n');
    fprintf(fid, '===========================\n');
end

subjectIDs = fieldnames(subjectInfo);
status = struct('Subject',{},'DaysLoaded',{},'DaysMissing',{},'DaysFailed',{}, ...
                'Warnings',{},'SavePath',{});

% -------------------------------------------------------------------------
% per subject
% -------------------------------------------------------------------------
for s = 1:numel(subjectIDs)
    subjID = subjectIDs{s};
    dayTags = subjectInfo.(subjID);

    if ~iscell(dayTags)
        error('subjectInfo.%s must be a cell array of day tags.', subjID);
    end

    setFolder = fullfile(rootPath, subjID, [subjID '_setfiles']);

    if fid > 0
        fprintf(fid, '\nSubject: %s\n', subjID);
    end

    if ~exist(setFolder,'dir')
        warning('Missing setfiles folder: %s', setFolder);
        if fid > 0
            fprintf(fid, 'Missing setfiles folder.\n');
        end
        status(end+1) = struct( ...
            'Subject', subjID, ...
            'DaysLoaded', 0, ...
            'DaysMissing', numel(dayTags), ...
            'DaysFailed', 0, ...
            'Warnings', "Missing setfiles folder", ...
            'SavePath', "" ); %#ok<AGROW>
        continue
    end

    Mono = struct('EEG',{},'DayStr',{});
    Meta = struct( ...
        'Subject', subjID, ...
        'DayTags', {dayTags}, ...
        'DayIndex', [], ...
        'LabelsByDay', {cell(numel(dayTags),1)}, ...
        'OrigLabelsByDay', {cell(numel(dayTags),1)}, ...
        'NumericLabelsByDay', {cell(numel(dayTags),1)}, ...
        'ChannelSideByDay', {cell(numel(dayTags),1)}, ...
        'LateralityRelPainByDay', {cell(numel(dayTags),1)}, ...
        'SideLR_NumByDay', {cell(numel(dayTags),1)}, ...
        'IpsiContra_NumByDay', {cell(numel(dayTags),1)}, ...
        'RegionByDay', {cell(numel(dayTags),1)}, ...
        'srate', [], ...
        'nbchan', [], ...
        'ToolInfo', struct( ...
            'MatlabVersion', version, ...
            'FunctionName', mfilename, ...
            'RunTimestamp', char(datetime('now')) ) );

    if isfield(painSideMap, subjID)
        PainSide = upper(string(painSideMap.(subjID)));
    else
        PainSide = "";
    end

    Ratings = table(strings(0,1), zeros(0,1), nan(0,1), nan(0,1), ...
        'VariableNames', {'DayStr','DayIndex','Pain_Intensity','Unpleasantness'});

    missingDays = {};
    failedLoads = {};

    for d = 1:numel(dayTags)
        dayTag = string(dayTags{d});

        [fullFilePath, filename] = find_set_file(setFolder, subjID, dayTag);

        if fullFilePath == ""
            warning('Missing .set file for %s %s', subjID, dayTag);
            if fid > 0
                fprintf(fid, 'Missing .set file for day index %d.\n', d);
            end
            missingDays{end+1} = char(dayTag); %#ok<AGROW>
            Ratings = [Ratings; {dayTag, d, NaN, NaN}]; %#ok<AGROW>
            continue
        end

        try
            EEG = pop_loadset('filename', filename, 'filepath', setFolder);
        catch ME
            warning('Load failed for %s: %s', fullFilePath, ME.message);
            if fid > 0
                fprintf(fid, 'Load failed for day index %d: %s\n', d, ME.message);
            end
            failedLoads{end+1} = char(dayTag); %#ok<AGROW>
            Ratings = [Ratings; {dayTag, d, NaN, NaN}]; %#ok<AGROW>
            continue
        end

        EEG = ensure_chanloc_labels(EEG);

        origLabels      = strings(EEG.nbchan,1);
        canonLabels     = strings(EEG.nbchan,1);
        numericAlias    = zeros(EEG.nbchan,1);
        channelSide     = strings(EEG.nbchan,1);
        painLaterality  = strings(EEG.nbchan,1);
        sideNum         = zeros(EEG.nbchan,1);
        ipsiContraNum   = zeros(EEG.nbchan,1);

        for ch = 1:EEG.nbchan
            lab0 = "";
            if ch <= numel(EEG.chanlocs)
                lab0 = string(EEG.chanlocs(ch).labels);
            end

            origLabels(ch) = lab0;
            [canon, numAlias] = canonicalize_monopolar_label(lab0, ch);

            canonLabels(ch) = canon;
            numericAlias(ch) = numAlias;
            EEG.chanlocs(ch).labels = char(canon);

            channelSide(ch) = side_from_label(canon);
            painLaterality(ch) = laterality_relative_to_pain(channelSide(ch), PainSide);
            sideNum(ch) = double(channelSide(ch) == "R");
            ipsiContraNum(ch) = double(painLaterality(ch) == "contra");
        end

        regionByChan = {};
        if ~isempty(regionDef)
            try
                regionByChan = label_regions_from_def(canonLabels, regionDef);
            catch
                regionByChan = {};
            end
        end

        Mono(end+1).EEG = EEG; %#ok<SAGROW>
        Mono(end).DayStr = char(dayTag);

        Meta.DayIndex(end+1) = d;
        Meta.srate(end+1) = EEG.srate;
        Meta.nbchan(end+1) = EEG.nbchan;
        Meta.LabelsByDay{d} = cellstr(canonLabels(:));
        Meta.OrigLabelsByDay{d} = cellstr(strtrim(origLabels(:)));
        Meta.NumericLabelsByDay{d} = numericAlias(:);
        Meta.ChannelSideByDay{d} = cellstr(channelSide(:));
        Meta.LateralityRelPainByDay{d} = cellstr(painLaterality(:));
        Meta.SideLR_NumByDay{d} = sideNum(:);
        Meta.IpsiContra_NumByDay{d} = ipsiContraNum(:);
        if ~isempty(regionByChan)
            Meta.RegionByDay{d} = regionByChan(:);
        else
            Meta.RegionByDay{d} = {};
        end

        Ratings = [Ratings; {dayTag, d, NaN, NaN}]; %#ok<AGROW>
    end

    Ratings = attach_daily_ratings(Ratings, rootPath, subjID, subjectRatings);

    if ~isempty(Mono)
        savePath = fullfile(setFolder, [subjID '_RawEEG.mat']);
        EEG_list = {Mono.EEG}; %#ok<NASGU>
        save(savePath, 'Mono', 'Meta', 'PainSide', 'Ratings', 'EEG_list', '-v7.3');

        if opts.writeSummaries
            summaryPath = fullfile(setFolder, [subjID '_RawEEG_Summary.txt']);
            write_subject_summary(summaryPath, subjID, dayTags, missingDays, failedLoads, Meta, Ratings);
        end

        if opts.writeMergedRatings
            ratingsCsv = fullfile(setFolder, [subjID '_Ratings_Merged.csv']);
            writetable(Ratings, ratingsCsv);
        end

        if fid > 0
            validate_subject_load(Meta, Ratings, fid);
            fprintf(fid, 'Saved RawEEG for %s. Days loaded: %d\n', subjID, numel(Mono));
        end

        warns = string(unique([missingDays, failedLoads]));
        if isempty(warns)
            warns = "";
        else
            warns = strjoin(warns, '; ');
        end

        status(end+1) = struct( ...
            'Subject', subjID, ...
            'DaysLoaded', numel(Mono), ...
            'DaysMissing', numel(missingDays), ...
            'DaysFailed', numel(failedLoads), ...
            'Warnings', warns, ...
            'SavePath', savePath ); %#ok<AGROW>
    else
        warning('No valid EEG loaded for %s', subjID);
        if fid > 0
            fprintf(fid, 'No valid EEG loaded.\n');
        end
        status(end+1) = struct( ...
            'Subject', subjID, ...
            'DaysLoaded', 0, ...
            'DaysMissing', numel(missingDays), ...
            'DaysFailed', numel(failedLoads), ...
            'Warnings', "No valid EEG data", ...
            'SavePath', "" ); %#ok<AGROW>
    end
end

if fid > 0
    fclose(fid);
end

fprintf('EEG loading and label normalization complete.\n');
end

% =========================================================================
% helpers
% =========================================================================

function [fullFilePath, filename] = find_set_file(setFolder, subjID, dayTag)
fullFilePath = "";
filename = "";

dayTag = char(string(dayTag));

candidates = { ...
    sprintf('%s_%s_Sel_2 resamp.set', subjID, dayTag), ...
    sprintf('%s_%s_Sel_2 resampled.set', subjID, dayTag), ...
    sprintf('%s_%s_Sel_2_f resampled.set', subjID, dayTag)};

for i = 1:numel(candidates)
    testPath = fullfile(setFolder, candidates{i});
    if exist(testPath, 'file') == 2
        fullFilePath = string(testPath);
        filename = candidates{i};
        return
    end
end

hits = dir(fullfile(setFolder, sprintf('%s_%s*.set', subjID, dayTag)));
hits = hits(~[hits.isdir]);
if ~isempty(hits)
    fullFilePath = string(fullfile(hits(1).folder, hits(1).name));
    filename = hits(1).name;
end
end

function EEG = ensure_chanloc_labels(EEG)
if ~isfield(EEG,'chanlocs') || isempty(EEG.chanlocs)
    EEG.chanlocs = repmat(struct('labels',''), EEG.nbchan, 1);
elseif numel(EEG.chanlocs) < EEG.nbchan
    EEG.chanlocs(end+1:EEG.nbchan,1) = struct('labels','');
elseif ~isfield(EEG.chanlocs,'labels')
    for k = 1:EEG.nbchan
        EEG.chanlocs(k).labels = '';
    end
end
end

function [canonLab, numLab] = canonicalize_monopolar_label(inputLabel, chIndex)
t = upper(string(inputLabel));
t = regexprep(t, '[\s\-\_\(\)]', '');
t = regexprep(t, '[:\.]+$', '');
t = regexprep(t, '(REF|REFC|AVG)$', '');

if strlength(t) == 0
    [canonLab, numLab] = fallback_label_from_index(chIndex);
    return
end

if all(isstrprop(t,'digit'))
    n = str2double(t);
    if n >= 1 && n <= 9
        canonLab = "L" + n; numLab = n; return
    elseif n >= 10 && n <= 18
        canonLab = "R" + (n - 9); numLab = n; return
    end
end

t = regexprep(t, '^LI', 'L');
t = regexprep(t, '^RI', 'R');

mL = regexp(t, '^L(\d+)$', 'tokens', 'once');
if ~isempty(mL)
    k = str2double(mL{1});
    if k >= 1 && k <= 9
        canonLab = "L" + k; numLab = k; return
    end
end

mR = regexp(t, '^R(\d+)$', 'tokens', 'once');
if ~isempty(mR)
    k = str2double(mR{1});
    if isnan(k) || k < 1 || k > 18
        [canonLab, numLab] = fallback_label_from_index(chIndex);
        return
    end
    if k >= 1 && k <= 9
        canonLab = "R" + k; numLab = 9 + k; return
    elseif k >= 10 && k <= 18
        canonLab = "R" + (k - 9); numLab = k; return
    end
end

[canonLab, numLab] = fallback_label_from_index(chIndex);
end

function [lab, num] = fallback_label_from_index(idx)
if idx <= 9
    lab = "L" + idx;
    num = idx;
else
    lab = "R" + (idx - 9);
    num = idx;
end
end

function s = side_from_label(canonLab)
if startsWith(canonLab, "L")
    s = "L";
elseif startsWith(canonLab, "R")
    s = "R";
else
    s = "";
end
end

function lc = laterality_relative_to_pain(sideChar, painSide)
sideChar = upper(string(sideChar));
painSide = upper(string(painSide));
if sideChar == "" || painSide == ""
    lc = "";
elseif sideChar == painSide
    lc = "ipsi";
else
    lc = "contra";
end
end

function Ratings = attach_daily_ratings(Ratings, rootPath, subjID, subjectRatings)
n = height(Ratings);

if isstruct(subjectRatings) && isfield(subjectRatings, subjID)
    R = subjectRatings.(subjID);
    if istable(R)
        dateCol = first_present(R, {'Date','Day','DayStr'});
        painCol = first_present(R, {'Pain_Intensity','Pain','Intensity'});
        unplCol = first_present(R, {'Unpleasantness','Unpl','Unpleasant'});
        assert(~isempty(dateCol) && ~isempty(painCol) && ~isempty(unplCol), ...
            'Ratings table must contain date, pain, and unpleasantness columns.');
        Tdates = normalize_dates(R.(dateCol));
        Ratings = inner_join_by_date(Ratings, Tdates, R.(painCol), R.(unplCol));
    elseif isfield(R,'Pain_Intensity') && isfield(R,'Unpleasantness')
        PI = R.Pain_Intensity(:);
        UN = R.Unpleasantness(:);
        if numel(PI) < n, PI(end+1:n,1) = NaN; end
        if numel(UN) < n, UN(end+1:n,1) = NaN; end
        Ratings.Pain_Intensity = PI(1:n);
        Ratings.Unpleasantness = UN(1:n);
    end
end

csvPath = fullfile(rootPath, subjID, sprintf('%s_DailyRatings.csv', subjID));
if exist(csvPath,'file') == 2
    try
        R = readtable(csvPath);
        dateCol = first_present(R, {'Date','Day','DayStr'});
        painCol = first_present(R, {'Pain_Intensity','Pain','Intensity'});
        unplCol = first_present(R, {'Unpleasantness','Unpl','Unpleasant'});
        assert(~isempty(dateCol) && ~isempty(painCol) && ~isempty(unplCol), ...
            'DailyRatings CSV must contain date, pain, and unpleasantness columns.');
        Tdates = normalize_dates(R.(dateCol));
        Ratings = inner_join_by_date(Ratings, Tdates, R.(painCol), R.(unplCol));
    catch ME
        warning('Failed to read %s: %s', csvPath, ME.message);
    end
end

try
    PI = read_root_rating_csv(fullfile(rootPath,'Pain_intensity.csv'), subjID);
    if ~isempty(PI)
        TPI = normalize_dates(PI.Date);
        Ratings = inner_join_by_date(Ratings, TPI, PI.Value, nan(height(TPI),1));
    end
    UN = read_root_rating_csv(fullfile(rootPath,'Unpleasantness.csv'), subjID);
    if ~isempty(UN)
        TUN = normalize_dates(UN.Date);
        Ratings = inner_join_by_date(Ratings, TUN, nan(height(TUN),1), UN.Value);
    end
catch ME
    warning('Root-level date merge failed for %s: %s', subjID, ME.message);
end

try
    PID = read_root_rating_by_day(fullfile(rootPath,'Pain_intensity.csv'), subjID, 'Pain_intensity');
    if ~isempty(PID)
        for r = 1:height(PID)
            di = PID.DayIndex(r);
            if di >= 1 && di <= n && ~isnan(PID.Value(r))
                Ratings.Pain_Intensity(di) = PID.Value(r);
            end
        end
    end

    UND = read_root_rating_by_day(fullfile(rootPath,'Unpleasantness.csv'), subjID, 'Unpleasantness');
    if ~isempty(UND)
        for r = 1:height(UND)
            di = UND.DayIndex(r);
            if di >= 1 && di <= n && ~isnan(UND.Value(r))
                Ratings.Unpleasantness(di) = UND.Value(r);
            end
        end
    end
catch ME
    warning('Root-level day-index merge failed for %s: %s', subjID, ME.message);
end
end

function out = read_root_rating_by_day(path, subjID, whichVar)
out = [];
if exist(path,'file') ~= 2, return; end
T = readtable(path);
v = lower(string(T.Properties.VariableNames));
iSubj = find(v=="subject" | v=="subj" | v=="id", 1);
iDay  = find(v=="day" | v=="dayindex" | v=="idx", 1);
iVal  = find(v==lower(string(whichVar)) | v=="value" | v=="rating" | v=="score" | ...
             v=="pain_intensity" | v=="painintensity" | v=="unpleasantness", 1);
if isempty(iSubj) || isempty(iDay) || isempty(iVal), return; end
mask = upper(string(T{:,iSubj})) == upper(string(subjID));
if ~any(mask), return; end
out = table(double(T{mask,iDay}), double(T{mask,iVal}), ...
    'VariableNames', {'DayIndex','Value'});
end

function out = read_root_rating_csv(path, subjID)
out = [];
if exist(path,'file') ~= 2, return; end
T = readtable(path);
v = lower(string(T.Properties.VariableNames));

iDate = find(v=="date" | v=="day" | v=="daystr" | v=="recordeddate", 1);
if isempty(iDate)
    for k = 1:numel(v)
        if contains(v(k), "date")
            iDate = k;
            break
        end
    end
end
if isempty(iDate), return; end

iSubj = find(strcmpi(T.Properties.VariableNames, subjID), 1);
if isempty(iSubj)
    cleaned = regexprep(string(T.Properties.VariableNames), '\W', '');
    target = regexprep(subjID, '\W', '');
    iSubj = find(strcmpi(cleaned, target), 1);
end

if ~isempty(iSubj)
    out = table(T{:,iDate}, double(T{:,iSubj}), 'VariableNames', {'Date','Value'});
    return
end

iS = find(v=="subject" | v=="subj" | v=="id", 1);
iV = find(v=="value" | v=="rating" | v=="score" | v=="pain_intensity" | v=="unpleasantness", 1);
if ~isempty(iS) && ~isempty(iV)
    mask = upper(string(T{:,iS})) == upper(string(subjID));
    if ~any(mask), return; end
    out = table(T{mask,iDate}, double(T{mask,iV}), 'VariableNames', {'Date','Value'});
end
end

function name = first_present(T, candidates)
name = '';
for k = 1:numel(candidates)
    if ismember(candidates{k}, T.Properties.VariableNames)
        name = candidates{k};
        return
    end
end
end

function Tdates = normalize_dates(DateCol)
s = string(DateCol);
DayStr = strings(numel(s),1);

for i = 1:numel(s)
    si = strtrim(s(i));
    dt = NaT;

    if contains(si,'/') || contains(si,'-')
        patterns = {'MM/dd/yyyy','M/d/yyyy','MM-d-yyyy','M-d-yyyy', ...
                    'yyyy-MM-dd','yyyy/M/d','yyyy/MM/dd','dd/MM/yyyy','d/M/yyyy'};
        for p = 1:numel(patterns)
            try
                dt = datetime(si,'InputFormat',patterns{p});
                if ~isnat(dt), break; end
            catch
            end
        end
    else
        d = regexprep(si,'\D','');
        if strlength(d) == 8
            try
                dt = datetime(d,'InputFormat','MMddyyyy');
            catch
            end
        end
    end

    if isnat(dt)
        try
            dt = datetime(si);
        catch
        end
    end

    if isnat(dt)
        DayStr(i) = regexprep(si,'\D','');
    else
        DayStr(i) = string(datestr(dt,'mmddyyyy'));
    end
end

Tdates = table(DayStr, 'VariableNames', {'DayStr'});
end

function Ratings = inner_join_by_date(Ratings, Tdates, PI, UNP)
L = height(Ratings);
key = table(string(Ratings.DayStr), 'VariableNames', {'DayStr'});
Tdates.DayStr = string(Tdates.DayStr);

J = outerjoin(key, [Tdates table(PI(:), UNP(:), ...
    'VariableNames', {'Pain_Intensity','Unpleasantness'})], ...
    'Keys', 'DayStr', 'MergeKeys', true, 'Type', 'left');

if height(J) ~= L
    return
end

if ~ismember('Pain_Intensity', Ratings.Properties.VariableNames)
    Ratings.Pain_Intensity = NaN(L,1);
end
if ~ismember('Unpleasantness', Ratings.Properties.VariableNames)
    Ratings.Unpleasantness = NaN(L,1);
end

if ismember('Pain_Intensity', J.Properties.VariableNames)
    idx = ~isnan(J.Pain_Intensity);
    Ratings.Pain_Intensity(idx) = J.Pain_Intensity(idx);
end
if ismember('Unpleasantness', J.Properties.VariableNames)
    idx = ~isnan(J.Unpleasantness);
    Ratings.Unpleasantness(idx) = J.Unpleasantness(idx);
end
end

function regionByChan = label_regions_from_def(canonLabels, regionDef)
regionByChan = repmat({''}, numel(canonLabels), 1);
lists = {'LeftAnterior','LeftPosterior','RightAnterior','RightPosterior'};
names = {'Anterior','Posterior','Anterior','Posterior'};

for i = 1:numel(canonLabels)
    lab = char(canonLabels(i));
    assigned = false;
    for j = 1:numel(lists)
        if isfield(regionDef, lists{j})
            L = regionDef.(lists{j});
            labShort = string(regexprep(lab,'^LI','L'));
            labShort = string(regexprep(labShort,'^RI','R'));
            if ismember(labShort, string(L)) || ismember(string(lab), string(L))
                regionByChan{i} = names{j};
                assigned = true;
                break
            end
        end
    end
    if ~assigned
        regionByChan{i} = '';
    end
end
end

function write_subject_summary(path, subjID, dayTags, missingDays, failedLoads, Meta, Ratings)
try
    fidS = fopen(path, 'w');
    if fidS < 0, return; end
    fprintf(fidS, 'Subject: %s\n', subjID);
    fprintf(fidS, 'Days expected: %d\n', numel(dayTags));
    fprintf(fidS, 'Days loaded: %d\n', numel(Meta.DayIndex));
    fprintf(fidS, 'Days missing: %d\n', numel(missingDays));
    fprintf(fidS, 'Days failed: %d\n', numel(failedLoads));
    if ~isempty(missingDays), fprintf(fidS, 'Missing count detail available in logs.\n'); end
    if ~isempty(failedLoads), fprintf(fidS, 'Failed count detail available in logs.\n'); end
    if ~isempty(Meta.srate), fprintf(fidS, 'Unique sampling rates: %s\n', strjoin(string(unique(Meta.srate)), ', ')); end
    if ~isempty(Meta.nbchan), fprintf(fidS, 'Unique channel counts: %s\n', strjoin(string(unique(Meta.nbchan)), ', ')); end
    nPI = sum(~isnan(Ratings.Pain_Intensity));
    nUN = sum(~isnan(Ratings.Unpleasantness));
    fprintf(fidS, 'Pain ratings available: %d of %d\n', nPI, height(Ratings));
    fprintf(fidS, 'Unpleasantness ratings available: %d of %d\n', nUN, height(Ratings));
    fclose(fidS);
catch
end
end

function validate_subject_load(Meta, Ratings, logFID)
fprintf(logFID, '\n==== Subject audit: %s ====\n', Meta.Subject);
fprintf(logFID, 'Days loaded: %d\n', numel(Meta.DayIndex));
if ~isempty(Meta.srate)
    fprintf(logFID, 'Unique sampling rates: %s\n', strjoin(string(unique(Meta.srate)), ', '));
end
if ~isempty(Meta.nbchan)
    fprintf(logFID, 'Unique channel counts: %s\n', strjoin(string(unique(Meta.nbchan)), ', '));
end
nPI = sum(~isnan(Ratings.Pain_Intensity));
nUN = sum(~isnan(Ratings.Unpleasantness));
fprintf(logFID, 'Ratings coverage: Pain %d of %d, Unpleasantness %d of %d\n', ...
    nPI, height(Ratings), nUN, height(Ratings));
end