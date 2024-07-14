
% This MATLAB code performs feature extraction from EEG data,
% organizes the results, and saves them to specified directories.
%
% 1. See the feature_params.m file for explenation and parameter seting
%


clc, clear, close all
eeglab
feature_params % Load label parameter file

path = split(pwd(),'feature_extraction');
path = path{1};

DATA_PATH = "D:\Ben Gurion University of Negev Dropbox\Miriam Guen\Miriam Guendelman\Projects\dynamotypes\data\scm_data"; % TODO: Fill the SCM data path
% data_path = fullfile(path, 'example_data', 'source_data'); % Use this to run on example data

%SAVE_PATH = "D:\Ben Gurion University of Negev Dropbox\Miriam Guen\Miriam Guendelman\Projects\dynamotypes\data\scm_data";  % TODO: Fill the feature save path
SAVE_PATH =  fullfile(path, 'example_data', 'features', seizure_control); % TODO: Save path for the output

if ~exist(SAVE_PATH, 'dir')
    mkdir(SAVE_PATH)
end

% Saving the feature extraction selected parameters in the result folder
save(fullfile(SAVE_PATH, 'run_params.mat'))

FILTER = ['_filter_' num2str(highpass_val) '_' num2str(lowpass_val) '_'];
RUN_NAME = [simulate FILTER 'time_from_bif_' num2str(time_from_bif)];

file_list = dir(fullfile(DATA_PATH, '*.mat'));


file_number = length(file_list);

start_bif = table();
end_bif = table();

parfor i = 1:file_number
    name = file_list(i).name;
    n = split(name, ".mat");
    save_path = fullfile(SAVE_PATH, n{1}); % Without the .mat suffix
    file_name = fullfile(DATA_PATH, name);

    % Extract the features
    [scores, onset_score, offset_score, start_table, end_table] = extract_comp(file_name, seizure_control);

    % Arrange and save the file data
    scores.file_name = repmat(string(file_list(i).name), size(scores, 1), 1);
    scores.comp_num = (1:n_comp)';

    onset_score.file_name = repmat(string(file_list(i).name), size(onset_score, 1), 1);
    onset_score.comp_num = (1:n_comp)';

    offset_score.file_name = repmat(string(file_list(i).name), size(offset_score, 1), 1);
    offset_score.comp_num = (1:n_comp)';

    writetable(scores, fullfile(SAVE_PATH, [RUN_NAME '_' name '_general_score_variance_both.csv']))

    writetable(horzcat(start_table, onset_score), fullfile(SAVE_PATH, [RUN_NAME '_' name '_start_bif_variance_both.csv']))

    writetable(horzcat(end_table, offset_score), fullfile(SAVE_PATH, [RUN_NAME '_' name '_end_bif_variance_both.csv']))
end