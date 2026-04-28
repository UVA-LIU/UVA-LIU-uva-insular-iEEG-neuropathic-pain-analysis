function status = run_bipolar_qc(rootPath, exclusionFile, subjectIDs, opts)
% apply_bipolar_qc_and_merge_ratings
%
% Apply bipolar channel exclusions and embed day-level ratings.
%
% INPUT
%   rootPath       : root data directory
%   exclusionFile  : CSV file with columns:
%                    Subject, Day, Exclude
%   subjectIDs     : cell array of subject IDs
%
% OPTIONAL (opts)
%   opts.writeLogs        (default true)
%   opts.writeSummary     (default true)
%   opts.outputSuffix     (default '_BipolarEEG_QC_Post.mat')
%
% OUTPUT
%   <subj>_BipolarEEG_QC_Post.mat

% ---------- defaults ----------
if nargin < 1 || isempty(rootPath)
    error('rootPath is required.');
end
if nargin < 2 || isempty(exclusionFile)
    error('exclusionFile is required.');
end
if nargin < 3 || isempty(subjectIDs)
    error('subjectIDs must be provided.');
end
if nargin < 4 || isempty(opts)
    opts = struct();
end

if ~isfield(opts,'writeLogs'), opts.writeLogs = true; end
if ~isfield(opts,'writeSummary'), opts.writeSummary = true; end
if ~isfield(opts,'outputSuffix'), opts.outputSuffix = '_BipolarEEG_QC_Post.mat'; end

if exist(exclusionFile,'file') ~= 2
    error('Exclusion CSV not found: %s', exclusionFile);
end

% ---------- read exclusion table ----------
T = readtable(exclusionFile, 'TextType','string');
requiredCols = ["Subject","Day","Exclude"];
for c = requiredCols
    if ~ismember(c, T.Properties.VariableNames)
        error('CSV missing required column: %s', c);
    end
end

T.Subject = upper(strtrim(string(T.Subject)));
T.Day     = strtrim(string(T.Day));
T.Exclude = string(T.Exclude);

% ---------- load ratings (optional) ----------
Tpain = load_rating_table(fullfile(rootPath,'Pain_intensity.csv'));
Tunpl = load_rating_table(fullfile(rootPath,'Unpleasantness.csv'));

summaryRows = {};

for s = 1:numel(subjectIDs)

    subj = string(subjectIDs{s});
    inDir  = fullfile(rootPath, subj, subj + "_setfiles");
    inFile = fullfile(inDir,  subj + "_BipolarEEG.mat");

    if exist(inFile,'file') ~= 2
        warning('Missing file for %s', subj);
        continue;
    end

    S = load(inFile);
    if ~isfield(S,'EEG_bipolar_list') || isempty(S.EEG_bipolar_list)
        warning('No bipolar data for %s', subj);
        continue;
    end

    Elist = S.EEG_bipolar_list;
    Meta  = [];
    if isfield(S,'Meta'), Meta = S.Meta; end

    PainSide = "";
    if isfield(S,'PainSide'), PainSide = string(S.PainSide); end

    Tsub = T(T.Subject == subj, :);
    EEG_bipolar_QC = cell(size(Elist));

    % ---------- logging ----------
    fid = -1;
    if opts.writeLogs
        logPath = fullfile(inDir, subj + "_QC_Log.txt");
        fid = fopen(logPath,'w');
        fprintf(fid, 'QC log for %s\n', subj);
    end

    for d = 1:numel(Elist)

        E = Elist{d};
        if isempty(E) || ~isfield(E,'data')
            continue;
        end

        rawLabs = E.chanlabels;
        if isempty(rawLabs)
            rawLabs = compose("Ch%d",1:size(E.data,1));
        end
        rawLabs = cellstr(rawLabs);

        % ----- match exclusions -----
        maskDay = match_day_rows(Tsub.Day, d, Meta);
        dRows   = Tsub(maskDay,:);

        reqTokens = split_and_clean(dRows.Exclude);
        [canonReq, canonReqN] = canonicalize_pairs(reqTokens);

        [canonLabs, canonLabsN] = canonicalize_pairs(rawLabs);

        loc = match_pairs(canonReq, canonReqN, canonLabs, canonLabsN);
        exclIdx = unique(loc(loc>0));

        keepIdx = setdiff(1:size(E.data,1), exclIdx);

        % ----- build output -----
        EEGq = struct();
        EEGq.data       = E.data(keepIdx,:);
        EEGq.srate      = getfield_safe(E,'srate',NaN);
        EEGq.chanlabels = rawLabs(keepIdx);

        if isfield(E,'region'), EEGq.region = to_cellstr(E.region(keepIdx)); end
        if isfield(E,'side'),   EEGq.side   = to_cellstr(E.side(keepIdx)); end
        if isfield(E,'axis'),   EEGq.axis   = to_cellstr(E.axis(keepIdx)); end
        if isfield(E,'laterality'), EEGq.laterality = to_cellstr(E.laterality(keepIdx)); end

        % ----- ratings -----
        dayTag = get_day_tag(Meta, S, d);
        PI = get_rating_for_day(Tpain, subj, d, dayTag);
        UN = get_rating_for_day(Tunpl, subj, d, dayTag);

        EEGq.info = struct( ...
            'Subject', char(subj), ...
            'DayIndex', d, ...
            'Pain_Intensity', PI, ...
            'Unpleasantness', UN, ...
            'N_before', size(E.data,1), ...
            'N_kept', numel(keepIdx));

        EEG_bipolar_QC{d} = EEGq;

        if fid > 0
            fprintf(fid,'Day %d: kept %d, removed %d\n', d, numel(keepIdx), numel(exclIdx));
        end

        summaryRows(end+1,:) = {char(subj), d, numel(keepIdx), numel(exclIdx), PI, UN}; %#ok<AGROW>
    end

    if fid > 0
        fclose(fid);
    end

    outFile = fullfile(inDir, subj + opts.outputSuffix);
    save(outFile,'EEG_bipolar_QC','-v7.3');

    fprintf('Saved QC data for %s\n', subj);
end

% ---------- summary ----------
if opts.writeSummary && ~isempty(summaryRows)
    Sum = cell2table(summaryRows, 'VariableNames', ...
        {'Subject','DayIndex','N_Kept','N_Removed','Pain','Unpleasantness'});
    writetable(Sum, fullfile(rootPath,'Bipolar_QC_Summary.csv'));
end

status = [];
end