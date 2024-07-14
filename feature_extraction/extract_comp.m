% This function extracts features from the components for labeling.
% It also saves the features of each component for further analysis.

function [scores, onset_score, offset_score, start_table, end_table] = extract_comp(file_name, seizure_control)

% This function extracts the seizure start time and features from the signal
feature_params % load label parameter file
data_to_load = ['clean_', seizure_control, '_data'];

clean_data = load(file_name, data_to_load) ;
clean_data = clean_data.(data_to_load);
load(file_name,  'seizure_info');
load(file_name,  'file_info');

%% Get relevant information
% The number of components
n_components = size(clean_data.icaweights,1); % n components
% Get the sample rate
srate = clean_data.srate;
% Normalize component signal
data = clean_data.icaweights*clean_data.data;
% Standardize seizure data
zdata = zscore(data, 1, 2);
% Get the component brain score
scores = score_comp(clean_data);

%% Initialize feature table
var_names = {'bif_time', 'spectral_slope', 'spectral_intercept', 'spectral_rmse', 'spectral_rsquare', 'freq', ...
    'spectral_slope_low', 'spectral_slope_high', 'spectral_intercept_low', 'spectral_intercept_high', 'spectral_dfe', 'spectral_sse', ...
    'width_autocorr', 'height_autocorr','mean_autocorr', 'std_autocorr', 'entropy_autocorr', 'approx_entropy', 'spectral_entropy', ...
    'stds', 'mean', 'median', 'skewness', 'isi_mean', 'isi_std', 'peak_mean', 'peak_std', ...
    'rms_linear_slope', 'rms_linear_intercept', 'rms_linear_adjrsquare', 'rms_linear_rmse', ...
    'rms_linear_slope_low', 'rms_linear_slope_high', 'rms_linear_intercept_low', 'rms_linear_intercept_high', 'rms_linear_dfe', 'rms_linear_sse', ...
    'isi_linear_slope', 'isi_linear_intercept', 'isi_linear_adjrsquare', 'isi_linear_rmse', ...
    'isi_linear_slope_low','isi_linear_slope_high',  'isi_linear_intercept_low', 'isi_linear_intercept_high', 'isi_linear_dfe', 'isi_linear_sse', ...
    'peak_linear_slope', 'peak_linear_intercept', 'peak_linear_adjrsquare','peak_linear_rmse',  ...
    'peak_linear_slope_low', 'peak_linear_slope_high', 'peak_linear_intercept_low','peak_linear_intercept_high', 'peak_linear_dfe', 'peak_linear_sse',  ...
    'rms_suph_stretch', 'rms_suph_intercept', 'rms_suph_adjrsquare', 'rms_suph_rmse', ...
    'rms_suph_stretch_low', 'rms_suph_stretch_high', 'rms_suph_intercept_low', 'rms_suph_intercept_high', 'rms_suph_dfe', 'rms_suph_sse', ...
    'isi_snic_stretch', 'isi_snic_intercept', 'isi_snic_adjrsquare', 'isi_snic_rmse', ...
    'isi_snic_stretch_low', 'isi_snic_stretch_high', 'isi_snic_intercept_low', 'isi_snic_intercept_high', 'isi_snic_dfe', 'isi_snic_sse', ...
    'isi_sh_stretch', 'isi_sh_intercept', 'isi_sh_adjrsquare', 'isi_sh_rmse', ...
    'isi_sh_stretch_low', 'isi_sh_stretch_high', 'isi_sh_intercept_low', 'isi_sh_intercept_high', 'isi_sh_dfe', 'isi_sh_sse',...
    'peak_suph_stretch', 'peak_suph_intercept', 'peak_suph_adjrsquare','peak_suph_rmse',...
    'peak_suph_stretch_low', 'peak_suph_stretch_high', 'peak_suph_intercept_low', 'peak_suph_intercept_high', 'peak_suph_dfe', 'peak_suph_sse', };

score_vars = {'mara', 'brain', 'muscle', 'eye',	'heart', 'line', 'chan',	'other', 'score'};

%var_type = repmat("double", 1, length(var_names));

%% Get the relevant start and end times for analysis
time_from_bif = max([time_from_bif_isi, time_from_bif_peak, time_from_bif_rms]);

seizure_length = seizure_info.sagmant_duration_t - 2 * seizure_info.time_margin;

seizure_start = seizure_info.time_margin + ...
    seconds(seizure_info.onset - seizure_info.analysis_start_time);

if seizure_start < 0
    disp(seizure_start)
end

seizure_end = seizure_start + seizure_length;

%% Analyze start bifurcation
% The range to look for the onset bifurcation
start_time = [max(seizure_start -bif_margin,1),...
    min(seizure_start + bif_margin + time_from_bif, seizure_info.sagmant_duration_p)] ;

start_table = make_empty_table(n_components, var_names);%, var_type);

pre_ictal_var = {'pre_slow_wave_power', 'pre_high_freq_power', 'pre_ictal_entropy', ...
    'pre_n_deflections', 'pre_deflections_variability', 'pre_median_rate', ...
    'pre_last_direction', 'pre_samples_from_last_deflection'};

% pre_var_type = repmat("double", 1, length(pre_ictal_var));

pre_ictal_table = make_empty_table(n_components, pre_ictal_var);%, pre_var_type);
start_table = [start_table, pre_ictal_table];

if start_time(2) - start_time(1) > bif_margin
    
    onset_score = get_local_score(clean_data, start_time(1), start_time(2));
    
    for i =1:n_components

        data = zdata(i,ceil(start_time(1)*srate):ceil(start_time(2)*srate));

        [bif_start_time, start_table] =  analyze_bif(data, srate, 1, start_table, i);
        if  bif_start_time > 0
            bif_start_time = start_time(1) + bif_start_time;
            pre_data = zdata(i,max(ceil((bif_start_time - 4)*srate),1):ceil(bif_start_time*srate));
        else
            pre_data = zdata(i,max(ceil((start_time(1) - 4)*srate),1):ceil(start_time(1)*srate));
        end

        start_table(i, :).bif_time = bif_start_time;


        try
            start_table = pre_ictal_features(pre_data, srate, start_table, i);
        catch
            fprintf('failed preictal feature extraction')
        end
    end
else
    start_table(:, :).bif_time =  repmat(seizure_start-bif_margin,n_comp,1);
    onset_score = make_empty_table(n_components, score_vars );
end

% for the start quantify the pre ictal slow occilations

%pre_ictal_data = zdata((bif_sample - 4 * srate):bif_sample);
%[pxx, f] = pwelch(pre_ictal_data, srate, 10,1:0.01:40,srate);


%% Analyze end bifurcation
% The range to look for the offset bifurcation
end_time = [max(seizure_end - bif_margin - time_from_bif,1),...
    min(seizure_end + bif_margin, seizure_info.sagmant_duration_p)];



end_table = make_empty_table(n_components, var_names);%, var_type);

if end_time(2) - end_time(1) > bif_margin

    offset_score = get_local_score(clean_data, end_time(1), end_time(2));

    for i =1:n_components
        % end_table(i, :).file_name = file_name;
        data = zdata(i,ceil(end_time(1)*srate):ceil(end_time(2)*srate));

        [bif_end_time, end_table] =  analyze_bif(data, srate, 0, end_table, i);
        if  bif_end_time > 0
            bif_end_time = end_time(1) + bif_end_time;
        end

        end_table(i, :).bif_time =  bif_end_time ;
    end
else
    end_table(:, :).bif_time =  repmat(seizure_end + bif_margin, n_comp, 1);
    offset_score = make_empty_table(n_components, score_vars);
end




end


function feature_table = make_empty_table(n_componants, var_names)%, var_type)
feature_table = array2table( NaN(n_componants,  length(var_names)) ,'VariableNames', var_names);%,'VariableTypes',  var_type);
%table('Size', [n_componants , length(var_names)], ...
%  'VariableTypes',  var_type,...
%  'VariableNames',var_names);
end


function  start_table = pre_ictal_features(data, srate, start_table, i)
% This function adds the following pre-ictal features to the onset
% table {slow_wave_power, high_freq_power, pre_ictal_entropy, pre_ictal_phase, ...
%        n_deflections, deflections_variability, median_rate, last_direction, last_gap}

n = length(data); % Length of the data
f = (0:n-1)*(srate/n); % Frequency range
Y = fft(data); % Compute FFT
% Compute Power Spectral Density (PSD)
P = abs(Y).^2/n; % Power spectral density

% Extract the desired frequency range (1-40 Hz)
freq_range = f >= 1 & f <= 40; % Logical index for 1-40 Hz
spectral_power = P(freq_range); % Extract spectral power in 1-40 Hz
f = f(freq_range);

% Sum the spectral power in the 1-40 Hz range
total_power = sum(spectral_power);

relative_power = spectral_power / total_power;

% for the bistability point test if bifurcations are related to a state with predominant slow waves
slow_freq_range = f >= 1 & f <= 4; % Logical index for 1-4 Hz
start_table(i, :).pre_slow_wave_power = sum(relative_power(slow_freq_range));

high_freq_range = f >= 30 & f <= 40; % Logical index for 30-40 Hz
start_table(i, :).pre_high_freq_power = sum(relative_power(high_freq_range));

start_table(i, :).pre_ictal_entropy = approximateEntropy(data);

% get slow wave phase direction (up-down)
% Design a low-pass filter below 4 Hz
filter_order = 4; % Filter order
cutoff_freq = 4; % Cutoff frequency in Hz
[b, a] = butter(filter_order, cutoff_freq/(srate/2), 'low'); % Butterworth filter design

% Apply the low-pass filter
filtered_data = filtfilt(b, a, data);
% Compute the Hilbert transform to get the analytical signal
analytic_signal = hilbert(filtered_data);

% Extract the deflections in the slow ocssilation (peaks)
phase_information = diff(sign(angle(analytic_signal)));
deflections = find(phase_information);
start_table(i, :).pre_n_deflections = length(deflections);
start_table(i, :).pre_deflections_variability = std(diff(deflections));
start_table(i, :).pre_median_rate = median(diff(deflections));
start_table(i, :).pre_last_direction = sign(phase_information(deflections(end)));
start_table(i, :).pre_samples_from_last_deflection = length(data) - deflections(end);

end

