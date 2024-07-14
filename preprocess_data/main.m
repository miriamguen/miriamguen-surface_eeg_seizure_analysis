% This code requiers the eeglab packege developed with eeglab2022.0 version and tested on matlb R2022a 
% The code performes the following: 
% 1. Automated selection of ictal segments based on seizure annotations,
%    using a time margin of 30% of the seizure lenght on each edge. 
%    This code was developed to mach the EPILEPSIAE binary EEG files, and a
%    patient_list.csv, and seizure_table.csv, and patient_file_list.csv
%    for each patient, example files are provided in the example data folder. 
%
% 2. Preprocessing ictal data, filtering and burst corrections.
% 3. Apply PCA-ICA decomposition on the time seiries
% 4. Randomly select a non ictal segment, of the same length as the ictal
%    filethat is at least 5 minutes away from the seizure onset and offset and
%    run the same analysis on it
% 5. Save the EEG lab data structure with the seizure information



clc, clear all, close all
eeglab

path = split(pwd(),'preprocess_data');
path = path{1};

data_path = fullfile(path, 'example_data', 'sorce_data'); %TODO: your data path
save_path =  fullfile(path, 'example_data', 'scm_data'); %TODO: save path for the output
mkdir(save_path);


% read the patient list 
patient_list = readtable(fullfile(data_path, 'patient_list.csv'));
history_all = 'LOG OF CONTROL AND SEIZURE DIFF';

for i=1:length(patient_list.code) 
    fprintf([patient_list.code{i} '\n'])
    temp = split(patient_list.code{i},'_'); %get the patient_id

    % read seizure list 
    temp_path = dir(fullfile(data_path, 'seizure_tables', [temp{end} '_*_seizures.csv']));
    seizure_list = readtable(fullfile(temp_path.folder, temp_path.name),'PreserveVariableNames',true);
    
    % read file list 
    temp_path = dir(fullfile(data_path, 'pat_file_tables', ['pat_' temp{end} '_*_file_list.csv']));
    file_list = readtable(fullfile(temp_path.folder, temp_path.name),'PreserveVariableNames',true);

    % try to find the seizures in the available data detect seizure 
 
    detect_seizure_and_save(seizure_list, file_list, save_path, data_path);

end

