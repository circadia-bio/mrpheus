% generate_fixtures.m
%
% Generates CSV test fixtures by running the original MATLAB pipeline on
% the synthetic Philips PMU log file.
%
% Instructions
% ------------
% 1. Drag these files into MATLAB Online:
%      generate_fixtures.m
%      ReadPhilipsScanPhysLog.m
%      pan_tompkin.m
%      example_physlog.log   (from mrpheus/inst/extdata/)
%
% 2. Run this script (press Run or type: run generate_fixtures.m)
%
% 3. Download the generated CSV files and place them in:
%      mrpheus/tests/testthat/fixtures/

LOG_FILE = 'example_physlog.log';
SFREQ    = 496;
ECG_COL  = 1;   % v1raw

if ~exist(LOG_FILE, 'file')
    error(['Log file not found: %s\n' ...
           'Make sure example_physlog.log is in the same folder.'], LOG_FILE);
end

fprintf('Reading: %s\n', LOG_FILE);
data = ReadPhilipsScanPhysLog(LOG_FILE);

n = size(data.C, 1);
fprintf('Samples : %d\n', n);
fprintf('Channels: %d\n', size(data.C, 2));

% Signal matrix
writematrix(data.C(1:100,      :), 'fixture_C_first100.csv');
writematrix(data.C(end-99:end, :), 'fixture_C_last100.csv');
fprintf('Written: fixture_C_first100.csv, fixture_C_last100.csv\n');

% Markers
writematrix(double(data.I.ScannerStart(:)), 'fixture_scanner_start.csv');
writematrix(double(data.I.ScannerStop(:)),  'fixture_scanner_stop.csv');
writematrix(double(data.I.VcgOnset(:)),     'fixture_vcg_onset.csv');
fprintf('Written: fixture_scanner_start.csv, fixture_scanner_stop.csv, fixture_vcg_onset.csv\n');

% QRS detection
ecg = double(data.C(:, ECG_COL));
fprintf('Running pan_tompkin...\n');
[qrs_amp_raw, qrs_i_raw, delay] = pan_tompkin(ecg, SFREQ, 0);

fprintf('R-peaks : %d\n', length(qrs_i_raw));
fprintf('Delay   : %d samples\n', delay);
fprintf('Mean HR : %.1f bpm\n', 60 * SFREQ / mean(diff(qrs_i_raw)));

writematrix(double(qrs_i_raw(:)),   'fixture_qrs_i_raw.csv');
writematrix(double(qrs_amp_raw(:)), 'fixture_qrs_amp_raw.csv');
fprintf('Written: fixture_qrs_i_raw.csv, fixture_qrs_amp_raw.csv\n');

fprintf('\nDone. Download all fixture_*.csv files and place in:\n');
fprintf('  mrpheus/tests/testthat/fixtures/\n');
