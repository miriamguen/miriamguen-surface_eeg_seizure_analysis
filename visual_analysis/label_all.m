%% manually label scm
clear all, close all
eeglab


reviewer_name = 'rev_1';


path = split(pwd(),'visual_analysis_app');
path = path{1};

data_path = fullfile(path, 'example_data', 'scm_data'); %TODO: your data path
label_path =  fullfile(path, 'example_data', 'labels');

mkdir(label_path)

file_list = dir(fullfile(data_path, '*.mat'));
file_names = extractfield(file_list, 'name');

skip_table = table;

% if this is not the first labeling session read the saved label file
label_file_path = fullfile(label_path, sprintf('%s_labels.xlsx', reviewer_name));
if exist(label_file_path, 'file') == 2
    opts = detectImportOptions(label_file_path);
    opts = setvartype(opts, {  'freq', 'scm', 'clear_start', 'clear_end'}, 'double');
    opts = setvartype(opts, {'comp_num', 'bif_start', 'bif_end','notes', 'file_name'},'string');
    label_table = readtable(label_file_path, opts);
    labeld_files = unique(label_table.('file_name'));
    file_names = setdiff(file_names, labeld_files);

else
    column_names = {'comp_num', 'bif_start', 'bif_end', 'freq', 'scm', 'clear_start', 'clear_end', 'notes', 'file_name'};
    variable_types = {'string', 'string', 'string', 'double', 'double', 'double', 'double', 'string', 'string'};
    label_table = table('Size', [0, length(column_names)], 'VariableTypes', variable_types, 'VariableNames', column_names);
end


f = msgbox(['you have ' num2str(size(file_names,2)) ' files left to review' ]);

for i = 1:length(file_names)
    disp(file_names{i})
    disp(i)
    % manualy label the componants that seem to includ patterns most similar
    file_name = fullfile(data_path ,file_names{i});
    load(file_name);
    

    title = ['origin: ' seizure_info.origin{1} ' wanted duration in sec: '...
        num2str(seizure_info.sagmant_duration_t)...
        'seizure starts at: ' num2str(seizure_info.time_margin) ' sec'...
        'with the pattern of:  ' seizure_info.pattern{1,1}];
    
    eegplot(clean_seizure_data.data,'srate',...
        clean_seizure_data.srate,'title', ['Clean data ' title],...
        'winlength',50,'position', [0,500,1200,500],...
        'eloc_file',clean_seizure_data.chanlocs );
    
    
    eegplot(clean_seizure_data.icaweights*clean_seizure_data.data,'srate',...
        clean_seizure_data.srate, 'title',['ICs ' title],'winlength',...
        50, 'position', [0,100,1500,500]);
    
    pop_viewprops(clean_seizure_data, 0, [1:seizure_info.n_comp],...
        {'freqrange', seizure_info.bnd}, {}, 1, 'ICLabel');
    
    
    file_starts = file_info.start_ts;
    file_edns   = file_info.start_ts + seconds(file_info.duration_in_sec);
    sagment_onset = seizure_info.onset - seconds(seizure_info.time_margin);
    sagment_end   = seizure_info.offset + seconds(seizure_info.time_margin);
    
    clear dataTable
    tutorialApp; %after pressing send table the app sends the scm table as "dataTable"
    % opts.Interpreter = 'tex';
    if file_starts > sagment_onset
        mh = msgbox({['sagment starts:  ' datestr(sagment_onset)];...
            ['  file starts:   ' datestr(file_starts)]},...
            'Warninng');
        th = findall(mh, 'Type', 'Text');                   %get handle to text within msgbox
        th.FontSize = 12;
        mh.Resize = 'on';
    end
    if file_edns < sagment_end
        mh = msgbox({['sagment ends:  ' datestr(sagment_end)];
            ['   file ends:  ' datestr(file_edns)]},...
            'Warninng');
        th = findall(mh, 'Type', 'Text');                   %get handle to text within msgbox
        th.FontSize = 12;
        mh.Resize = 'on';
    end
    pause()
    close all

    name = file_names{i};
    
    dataTable.file_name = repmat(convertCharsToStrings(name),height(dataTable),1);
    dataTable.('bif_start') = string(dataTable.('bif_start'));
    dataTable.('bif_end') = string(dataTable.('bif_end'));
    dataTable.("notes") = string(dataTable.("notes"));
    dataTable.("comp_num") = string(dataTable.("comp_num"));
    dataTable.("scm") = double(dataTable.("scm"));
    dataTable.("clear_start") = double(dataTable.("clear_start"));
    dataTable.("clear_end") = double(dataTable.("clear_end"));

    label_table = vertcat(label_table, dataTable);
    writetable(label_table, label_file_path);
    
end
