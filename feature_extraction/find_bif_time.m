% The function find_bif_time identifies the bifurcation point (a point of change) 
% in the provided data and extracts analysis data around that point. 
% It uses change-point detection to find significant changes in variance 
% within the data and filters the data using lowpass and highpass filters.


function [bif_sample, analysis_data] = find_bif_time( ...
    data, srate, bif_margin)

    % data: the data of the bifurcation possible range
    % srate: the sample rate 
    % bif_margin: the margin around seizure onset label to look for the bifurcation
    % time_from_bif: the window size for features around the changepoint.

feature_params

time_from_bif = max([time_from_bif_isi, time_from_bif_peak, time_from_bif_rms]);

data = lowpass(data,lowpass_val, srate);
data = highpass(data,highpass_val, srate);

bif_data = data(1:(length(data) - srate * time_from_bif));

TF = ischange(bif_data,'variance', 'MaxNumChanges',1);

bif_sample = find(TF);
w_size =  time_from_bif * srate;

% If no changepoint was found, mark time with -1 and use the original onset for analysis
if isempty(bif_sample)
    bif_sample = -1;
    analysis_data = data(srate * bif_margin: min(srate * bif_margin + w_size, length(data)));
else
    analysis_data = data(bif_sample:min(bif_sample + w_size, length(data)));
end