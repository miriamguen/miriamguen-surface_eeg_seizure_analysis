function history = detect_seizure_and_save(fseizure, record_list, save_path, data_path)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input: fseizure - a table documating the patients' seizures               %
%        record_list - a table containing all of the patients' records      %
% The function saves as an mat file the seizure and control seizure         %
% information, including the signals, onset\offst, the channels... and more %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


history = '';
safty_space_sec = 5*60;         % in sec
time_margin_proportion = 0.3;   % as a fraction of the seizure length
record_indexs = fseizure.file_seizure_ind(:);

%set the Clean_data Parameters

n_comp = 8;                      % number of pca componants to use
bnd = [1 40];                     % the frequancy band to use
dsmp = 256;
reref = 1;                       % create avarage ref


for i_seiz = 1:size(fseizure,1) % going through all the seizures
    % onset & offset time of current seizure
    seiz_onset = fseizure.onset(i_seiz);
    seiz_offset = fseizure.offset(i_seiz);
    time_margin = seconds(seiz_offset- seiz_onset)*time_margin_proportion;

    % add 1 to get matlab index from python index
    row_idx = fseizure.file_seizure_ind(i_seiz) + 1;
    record_onset = record_list.start_ts(row_idx);
    f_rate = record_list.sample_freq(row_idx);

    % how long into the record in sec the seizure occured
    % in case the seizure begins vary close to the begining of the
    % recording this will set the first index to the 1st sec
    seconds_until_onset = max(seconds(seiz_onset - record_onset)-time_margin,1);
    % get the record duration
    rec_dur = record_list.duration_in_sec(row_idx);
    % how meny seconds until sample offset, maximal is the rec duration
    seconds_to_offset = min(seconds(seiz_offset-record_onset)+time_margin,rec_dur);
    % calc seizure duration + margins
    seiz_dur_sec =  seconds_to_offset - seconds_until_onset;

    % calc the current sample
    sample_onset = round(seconds_until_onset * f_rate);
    sample_offset = round(sample_onset + (seiz_dur_sec * f_rate));
    % not to exceed recording range
    sample_offset = min(sample_offset, rec_dur * f_rate);

    %Read the EEG file
    header_path = split(record_list.file_path(row_idx), '/');
    header_path = fullfile(data_path, 'data_files', header_path{end});
    [signals,Fs,chan_names,locs] = rawdat2load(header_path);

    % sanity check make sure we have the correct sample rate
    assert(Fs == f_rate,'somthing with the sample rate is wrong!!!');
    %(rec_dur * Fs) = size(signals,2);
    %gets the signals of the current seizure and get the numer of samples in it
    seizure_signals = signals(:,sample_onset:sample_offset);
    nsample_seiz = size(seizure_signals,2);

    % create control_signals
    % this variable indicates whether there are several seizures in the recording
    other_seizures = find(record_indexs == fseizure.file_seizure_ind(i_seiz));

    % Select all the clear intervals including safty space
    onset_all = seconds(fseizure.onset(other_seizures)- record_onset);
    onset_all = (onset_all -safty_space_sec)*Fs;
    onset_all = [onset_all',rec_dur*Fs];
    offset_all = seconds(fseizure.offset(other_seizures)- record_onset);
    offset_all = (offset_all+safty_space_sec)*Fs;
    offset_all = [1 offset_all'];
    % select the maximal interval
    [max_interval, ind] = max(onset_all - offset_all);
    if max_interval < nsample_seiz
        history = sprintf('Inter ictal time not sufficent for control in the recording');
        save(['history.', num2str(i_seiz),'.mat'], "history");
        continue
    end


    % get all the possible start indexes
    sample_range = offset_all(ind):(onset_all(ind)- nsample_seiz);
    N = numel(sample_range);
    randInd = randi(N, 1);
    idx_start = sample_range(randInd);  %randsample(randInd ,1);
    %create the control signals
    ind_end = min((idx_start+nsample_seiz-1),record_list.num_samples(row_idx));
    control_signals = signals(:,idx_start:ind_end);

    %% clean the signal
    icawin = size(seizure_signals,2);% the window_size of the ica


    % the clean seizure data with labeld ics
    info = record_list.file_path(row_idx,1);
    info = info{1,1};
    try
        [clean_seizure_data,~] = clean_data(seizure_signals,Fs,...
            n_comp,bnd,dsmp,icawin,reref, chan_names,locs);
    catch
        history = sprintf('file:\n %s Patient seizure %s\n  Failed preprocessing seizure data \n',info ,i_seiz);
        save(['history.', num2str(i_seiz),'.mat'], "history");
        continue
    end
    % the clean control data with labeled ics
    try
        [clean_control_data,~] = clean_data(control_signals,Fs,...
            n_comp,bnd,dsmp,icawin,reref, chan_names,locs);
        assert(size(clean_seizure_data.data,1) == size(clean_control_data.data,1))
    catch
        history = sprintf('file:\n %s Patient seizure %s\n  Failed preprocessing control data \n',info ,i_seiz);
        save(['history.','_no_control_', num2str(i_seiz),'.mat'], "history");
        continue
    end

    %The relevant seizure info to save
    seizure_info = fseizure(i_seiz,:); %save the relevant seizure info
    %The wanted dureation
    seizure_info.sagmant_duration_t = seconds(seiz_offset-seiz_onset)+2*time_margin;
    %The duration in practice
    seizure_info.sagmant_duration_p = seiz_dur_sec;
    %The time margin used
    seizure_info.time_margin = time_margin;
    %the band used
    seizure_info.bnd = bnd;
    %the number of componants used
    seizure_info.n_comp = n_comp;

    %get the relevant file info
    file_info = record_list(row_idx,:);
    file_info.seconds_until_onset = seconds_until_onset;
    file_info.seconds_to_offset = seconds_to_offset;

    % the template is:  surf30_pat_1_seiz_1
    set = record_list.file_path(row_idx);
    set = split(set{1}, '/raw_data/');
    set = split(set{2}, '/');
    hospital = set{1,1};
    patient = set{2,1};
    id = erase(set{5,1}, '.head');

    % check this out - this is a better string formating method
    dataName = sprintf("%s_%s_%d_id_%s.mat",hospital,patient,i_seiz,id );

    % Save the dataset
    try
        save(fullfile(save_path, dataName),'clean_seizure_data',...
            'clean_control_data','seizure_info', 'file_info');
    catch
        history = sprintf('problem with designated path - please save manually');
        save(['history.',patient , num2str(i_seiz),'.mat'], "history");
        pause
    end


end
end
