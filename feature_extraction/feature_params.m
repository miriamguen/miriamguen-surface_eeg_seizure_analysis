%% Run Parameters
n_comp = 8;

seizure_control = 'control';          % Control or seizure 
simulate = 'raw';                     % The run alias, e.g., 'raw', 'raw_filtered'
time_from_bif = 'example';            % The time to analyze (historically), empirically modified per feature (see below)

%% Feature Extraction Parameters
lowpass_val = 40;                     % The lowpass filter value for the data
highpass_val = 1;                     % The highpass filter value for the data
bif_margin = 5;                       % The margin around start and end time to look in seconds
min_peak_distance = 1 / lowpass_val;  % The minimum distance allowed between subsequent identified peaks
time_from_bif_isi = 5;                % Time in seconds for the ISI, lower or eqal to bif_margin
time_from_bif_peak = 2;               % Time for peak hight analysis from bifurcation in seconds, lower or eqal to bif_margin
time_from_bif_rms = 2;                % The RMS envelope time parameter in seconds, lower or eqal to bif_margin