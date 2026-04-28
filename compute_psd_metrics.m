function status = compute_epoch_psd_metrics(rootPath, subjectIDs, opts)
% compute_epoch_psd_metrics
%
% Compute epoch-level power spectral density (PSD) metrics from
% *_BipolarEEG_QC_Post.mat files.
%
% INPUT
%   rootPath    : root data directory
%   subjectIDs  : cell array of subject IDs
%
% OPTIONAL (opts)
%   opts.epochLenSec  (default = 5)
%   opts.maxEpochs    (default = 60)
%   opts.welch        struct with fields:
%                       .win, .noverlap, .nfft, .fmax
%   opts.writeOutputs (default = true)
%
% OUTPUT
%   Per subject:
%     *_PSD_PerChannel.csv
%     *_PSD_ByRegion_Day.csv
%     *_PSD_EpochLong.csv
%
% NOTES
%   Requires Signal Processing Toolbox (pwelch)

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
if ~isfield(opts,'maxEpochs'),   opts.maxEpochs   = 60; end
if ~isfield(opts,'welch')
    opts.welch = struct('win',2048,'noverlap',1024,'nfft',8192,'fmax',55);
end
if ~isfield(opts,'writeOutputs'), opts.writeOutputs = true; end

bands = struct( ...
    'delta',[1 4], ...
    'theta',[4 8], ...
    'alpha',[8 13], ...
    'beta',[13 30], ...
    'gamma',[30 min(55,opts.welch.fmax)]);

bandNames = fieldnames(bands)';

status = struct('Subject',{},'Days',{},'EpochRows',{},'Warnings',{});

% ================================================================
for s = 1:numel(subjectIDs)

    subj = string(subjectIDs{s});
    inDir  = fullfile(rootPath, subj, subj + "_setfiles");
    qcFile = fullfile(inDir, subj + "_BipolarEEG_QC_Post.mat");

    if exist(qcFile,'file') ~= 2
        warning('Missing QC file for %s', subj);
        continue;
    end

    L = load(qcFile);
    if ~isfield(L,'EEG_bipolar_QC') || isempty(L.EEG_bipolar_QC)
        warning('No bipolar QC data for %s', subj);
        continue;
    end

    EEGq = L.EEG_bipolar_QC;

    outDir = fullfile(rootPath, subj, 'PSD_Epochs');
    if opts.writeOutputs && ~exist(outDir,'dir')
        mkdir(outDir);
    end

    epochLong = [];

    for d = 1:numel(EEGq)

        E = EEGq{d};
        if isempty(E) || ~isfield(E,'data')
            continue;
        end

        X  = double(E.data);
        fs = E.srate;
        [nChan, nSamp] = size(X);

        if nSamp < fs
            continue;
        end

        Nepoch = round(opts.epochLenSec * fs);
        nMax   = min(opts.maxEpochs, floor(nSamp / Nepoch));
        if nMax < 1
            continue;
        end

        X = X(:,1:(Nepoch*nMax));

        % Welch setup
        w = opts.welch;
        wlen = min(w.win, Nepoch);
        ovlp = min(w.noverlap, floor(wlen/2));
        nfft = max(w.nfft, 2^nextpow2(wlen));
        fmax = min(w.fmax, fs/2);

        [~, fHz] = pwelch(zeros(Nepoch,1), hamming(wlen), ovlp, nfft, fs);
        useK = (fHz >= 0 & fHz <= fmax);
        fUse = fHz(useK);

        for c = 1:nChan

            xe = X(c,:) - mean(X(c,:), 'omitnan');
            xe = reshape(xe, Nepoch, nMax);

            for e = 1:nMax

                x = detrend(xe(:,e));
                [pxx,~] = pwelch(x, hamming(wlen), ovlp, nfft, fs);

                pUse = pxx(useK);
                pUse(pUse<=0) = eps;

                for b = bandNames
                    band = b{1};
                    fr = bands.(band);
                    mask = fUse >= fr(1) & fUse < fr(2);

                    if nnz(mask) < 2
                        val = NaN;
                    else
                        val = 10*log10(trapz(fUse(mask), pUse(mask)));
                    end

                    row = struct();
                    row.Subject = subj;
                    row.DayIndex = d;
                    row.Channel = c;
                    row.Epoch = e;
                    row.Band = band;
                    row.Value = val;

                    epochLong = [epochLong; row]; %#ok<AGROW>
                end
            end
        end
    end

    if isempty(epochLong)
        continue;
    end

    Tlong = struct2table(epochLong);

    if opts.writeOutputs
        outFile = fullfile(outDir, subj + "_PSD_EpochLong.csv");
        writetable(Tlong, outFile);
        fprintf('Saved PSD for %s\n', subj);
    end

    status(end+1) = struct( ...
        'Subject', subj, ...
        'Days', numel(EEGq), ...
        'EpochRows', height(Tlong), ...
        'Warnings', "" ); %#ok<AGROW>
end

end