function T = build_connectivity_longtable(rootPath, subjectIDs, opts)
% build_connectivity_epoch_longtable
%
% Compute epoch-level connectivity (PLV + coherence) and output a
% long-format table for statistical modeling.
%
% INPUT
%   rootPath    : root data directory
%   subjectIDs  : cell array of subject IDs
%
% OPTIONAL (opts)
%   opts.epochLenSec  (default = 5)
%   opts.maxEpochs    (default = 60)
%   opts.welch        struct (win, noverlap, nfft, fmax)
%   opts.autosave     (default = true)
%
% OUTPUT
%   Combined table T
%
% SAVES
%   CONN_Epoch_Level_LongTable.csv
%   CONN_Epoch_Level_LongTable.mat

% ---------------- defaults ----------------
if nargin < 1 || isempty(rootPath)
    error('rootPath is required.');
end
if nargin < 2 || isempty(subjectIDs)
    error('subjectIDs must be provided.');
end
if nargin < 3 || isempty(opts)
    opts = struct();
end

if ~isfield(opts,'epochLenSec'), opts.epochLenSec = 5; end
if ~isfield(opts,'maxEpochs'),   opts.maxEpochs = 60; end
if ~isfield(opts,'welch')
    opts.welch = struct('win',2048,'noverlap',1024,'nfft',8192,'fmax',55);
end
if ~isfield(opts,'autosave'), opts.autosave = true; end

% ---------------- bands ----------------
bands = struct( ...
    'delta',[1 4], ...
    'theta',[4 8], ...
    'alpha',[8 13], ...
    'beta',[13 30], ...
    'gamma',[30 opts.welch.fmax]);

bnames = fieldnames(bands)';

RowBuffer = {};
RowCount  = 0;

% ================================================================
for s = 1:numel(subjectIDs)

    subj = string(subjectIDs{s});
    inDir  = fullfile(rootPath, subj, subj + "_setfiles");
    qcFile = fullfile(inDir, subj + "_BipolarEEG_QC_Post.mat");

    if exist(qcFile,'file') ~= 2
        warning('Missing QC file for %s', subj);
        continue;
    end

    S = load(qcFile);
    if ~isfield(S,'EEG_bipolar_QC') || isempty(S.EEG_bipolar_QC)
        warning('No bipolar QC data for %s', subj);
        continue;
    end

    EEGq = S.EEG_bipolar_QC;

    % determine pain side
    PainSide = "L";
    if isfield(S,'PainSide') && contains(upper(string(S.PainSide)),"R")
        PainSide = "R";
    end

    % ============================================================
    for d = 1:numel(EEGq)

        E = EEGq{d};
        if isempty(E) || ~isfield(E,'data')
            continue;
        end

        X = double(E.data);
        fs = getfield_safe(E,'srate',NaN);

        [nChan,nSamp] = size(X);
        if nSamp < fs
            continue;
        end

        % epoching
        Nepoch = round(opts.epochLenSec * fs);
        nMax   = min(opts.maxEpochs, floor(nSamp / Nepoch));
        if nMax < 1
            continue;
        end

        ind = reshape(1:(Nepoch*nMax), Nepoch, nMax);

        % metadata
        labs  = string(get_labels_str(E,nChan));
        side  = string(get_meta_str(E,'side',nChan));
        axis  = string(get_meta_str(E,'axis',nChan));
        lat   = string(get_meta_str(E,'laterality',nChan));
        region= string(get_meta_str(E,'region',nChan));

        % region groups
        [LA,LP,RA,RP] = region_sets_str(region,nChan);

        if PainSide=="L"
            IpsiA=LA; IpsiP=LP; ContraA=RA; ContraP=RP;
        else
            IpsiA=RA; IpsiP=RP; ContraA=LA; ContraP=LP;
        end

        PairSets = build_pairsets(IpsiA,IpsiP,ContraA,ContraP);

        % Welch grid
        [fUse,useK] = build_freq_grid(Nepoch,opts.welch,fs);

        for ps = 1:numel(PairSets)

            P = PairSets{ps};
            idxPairs = generate_pairs(P);

            for pidx = 1:size(idxPairs,1)

                i1 = idxPairs(pidx,1);
                i2 = idxPairs(pidx,2);

                for e = 1:nMax

                    seg = ind(:,e);

                    x = X(i1,seg) - mean(X(i1,seg));
                    y = X(i2,seg) - mean(X(i2,seg));

                    [Cxy,~] = mscohere(x,y,[],[],[],fs);
                    cUse = Cxy(useK);

                    for bn = bnames
                        band = bn{1};
                        mask = (fUse >= bands.(band)(1) & fUse < bands.(band)(2));

                        % COH
                        row = struct();
                        row.Subject = subj;
                        row.DayIndex = d;
                        row.PairID = labs(i1) + "_" + labs(i2);
                        row.PairType = P.name;
                        row.Band = band;
                        row.Metric = "coh";
                        row.Value = mean(cUse(mask),'omitnan');

                        RowCount = RowCount + 1;
                        RowBuffer{RowCount} = row;

                        % PLV
                        row.Metric = "plv";
                        row.Value  = compute_plv(x,y,fs,bands.(band));

                        RowCount = RowCount + 1;
                        RowBuffer{RowCount} = row;
                    end
                end
            end
        end
    end
end

% ================================================================
T = struct2table([RowBuffer{:}]);

% typing
T.Subject = categorical(T.Subject);
T.PairType = categorical(T.PairType);
T.Band = categorical(T.Band);
T.Metric = categorical(T.Metric);

% save
if opts.autosave
    outDir = fullfile(rootPath,'MixedEffects');
    if ~exist(outDir,'dir'), mkdir(outDir); end

    writetable(T, fullfile(outDir,'CONN_Epoch_Level_LongTable.csv'));
    save(fullfile(outDir,'CONN_Epoch_Level_LongTable.mat'),'T','-v7.3');

    fprintf('Saved connectivity table (%d rows)\n', height(T));
end
end