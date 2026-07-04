% data-raw/generate_matlab_fixtures.m
%
% Runs the original MATLAB pipeline on the synthetic Philips PMU log file
% bundled in inst/extdata/ and writes CSV fixtures for bit-perfect comparison
% against the R port in mrpheus.
%
% Run from the mrpheus package root:
%   cd /path/to/mrpheus
%   matlab -batch "run data-raw/generate_matlab_fixtures.m"
%
% Or interactively in MATLAB:
%   cd /path/to/mrpheus
%   run data-raw/generate_matlab_fixtures.m
%
% Requires ReadPhilipsScanPhysLog.m and pan_tompkin.m from neonatal_qrs.
% Set NEONATAL_QRS_PATH below to point at your local clone.
%
% Fixtures written to tests/testthat/fixtures/

% ---- Configuration --------------------------------------------------------
NEONATAL_QRS_PATH = fullfile('..', 'neonatal_qrs');  % adjust if needed
LOG_FILE          = fullfile('inst', 'extdata', 'example_physlog.log');
FIXTURE_DIR       = fullfile('tests', 'testthat', 'fixtures');
SFREQ             = 496;
ECG_COL           = 1;   % v1raw

% ---- Validate paths -------------------------------------------------------
if ~exist(LOG_FILE, 'file')
    error('Log file not found: %s\nRun from the mrpheus package root.', LOG_FILE);
end
if ~exist(NEONATAL_QRS_PATH, 'dir')
    error(['neonatal_qrs not found at: %s\n' ...
           'Edit NEONATAL_QRS_PATH in this script.'], NEONATAL_QRS_PATH);
end

addpath(NEONATAL_QRS_PATH);
mkdir(FIXTURE_DIR);

fprintf('Reading: %s\n', LOG_FILE);

% ---- Read physlog ---------------------------------------------------------
data = ReadPhilipsScanPhysLog(LOG_FILE);

n_samples  = size(data.C, 1);
n_channels = size(data.C, 2);
fprintf('Samples : %d\n', n_samples);
fprintf('Channels: %d\n', n_channels);

% Signal matrix — first and last 100 rows
writematrix(data.C(1:100,        :), fullfile(FIXTURE_DIR, 'fixture_C_first100.csv'));
writematrix(data.C(end-99:end,   :), fullfile(FIXTURE_DIR, 'fixture_C_last100.csv'));

% Marker indices
writematrix(data.I.ScannerStart(:), fullfile(FIXTURE_DIR, 'fixture_scanner_start.csv'));
writematrix(data.I.ScannerStop(:),  fullfile(FIXTURE_DIR, 'fixture_scanner_stop.csv'));
writematrix(data.I.VcgOnset(:),     fullfile(FIXTURE_DIR, 'fixture_vcg_onset.csv'));

fprintf('Written: C matrices and marker fixtures\n');

% ---- QRS detection --------------------------------------------------------
ecg = double(data.C(:, ECG_COL));
fprintf('Running pan_tompkin on ECG channel %d (v1raw)...\n', ECG_COL);

[qrs_amp_raw, qrs_i_raw, delay] = pan_tompkin(ecg, SFREQ, 0);

fprintf('R-peaks detected: %d\n', length(qrs_i_raw));
fprintf('Filter delay    : %d samples\n', delay);
fprintf('Mean HR         : %.1f bpm\n', 60 * SFREQ / mean(diff(qrs_i_raw)));

writematrix(qrs_i_raw(:),   fullfile(FIXTURE_DIR, 'fixture_qrs_i_raw.csv'));
writematrix(qrs_amp_raw(:), fullfile(FIXTURE_DIR, 'fixture_qrs_amp_raw.csv'));

fprintf('\nAll fixtures written to: %s\n', FIXTURE_DIR);
