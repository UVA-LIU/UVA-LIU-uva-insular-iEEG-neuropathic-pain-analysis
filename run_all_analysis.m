function run_all_analysis(rootPath, subjectIDs, subjectInfo, regionDef, exclusionFile, opts)
% run_all_analysis
%
% Master script for the insular iEEG chronic neuropathic pain analysis.
%
% This script runs the analysis pipeline in the order used for the manuscript.
%
% Required inputs:
%   rootPath       : root data directory
%   subjectIDs     : cell array of subject IDs
%   subjectInfo    : struct containing subject-specific recording day labels
%   regionDef      : struct defining insular contact regions
%   exclusionFile  : CSV file defining bipolar channel exclusions
%
% Optional input:
%   opts           : struct with optional subfields:
%                    opts.load
%                    opts.bipolar
%                    opts.qc
%                    opts.psd
%                    opts.conn
%
% Example:
%   run_all_analysis(rootPath, subjectIDs, subjectInfo, regionDef, exclusionFile, opts)

if nargin < 1 || isempty(rootPath)
    error('rootPath is required.');
end
if nargin < 2 || isempty(subjectIDs)
    error('subjectIDs is required.');
end
if nargin < 3 || isempty(subjectInfo)
    error('subjectInfo is required.');
end
if nargin < 4
    regionDef = [];
end
if nargin < 5 || isempty(exclusionFile)
    error('exclusionFile is required.');
end
if nargin < 6 || isempty(opts)
    opts = struct();
end

if ~isfield(opts,'load'),    opts.load = struct(); end
if ~isfield(opts,'bipolar'), opts.bipolar = struct(); end
if ~isfield(opts,'qc'),      opts.qc = struct(); end
if ~isfield(opts,'psd'),     opts.psd = struct(); end
if ~isfield(opts,'conn'),    opts.conn = struct(); end

fprintf('\n=== Running full insular iEEG analysis pipeline ===\n');

%% Step 1: Load raw EEGLAB set files
fprintf('\n[1/7] Loading raw EEG from EEGLAB set files...\n');
load_raw_eeg_from_setfiles(rootPath, subjectInfo, [], regionDef, opts.load);

%% Step 2: Build bipolar montage
fprintf('\n[2/7] Building bipolar montage...\n');
build_bipolar_montage(rootPath, subjectIDs, regionDef);

%% Step 3: Apply bipolar QC and merge ratings
fprintf('\n[3/7] Applying bipolar QC and merging ratings...\n');
run_bipolar_qc(rootPath, exclusionFile, subjectIDs, opts.qc);

%% Step 4: Compute PSD metrics
fprintf('\n[4/7] Computing PSD metrics...\n');
compute_psd_metrics(rootPath, subjectIDs, opts.psd);

%% Step 5: Compute connectivity long table
fprintf('\n[5/7] Computing connectivity metrics...\n');
build_connectivity_longtable(rootPath, subjectIDs, opts.conn);

%% Step 6: Run PSD mixed-effects models
fprintf('\n[6/7] Running PSD mixed-effects models...\n');
run_psd_lme_models([], rootPath);

%% Step 7: Run connectivity mixed-effects models
fprintf('\n[7/7] Running connectivity mixed-effects models...\n');
run_connectivity_lme_models([], rootPath);

fprintf('\n=== Full analysis pipeline complete ===\n');

end