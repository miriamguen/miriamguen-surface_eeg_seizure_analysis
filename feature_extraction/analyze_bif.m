%% This file recives a component, evaluates it for the seizure start time and extract featuere from this point

function [bif_time, feature_table] =  analyze_bif(data, srate, start_or_end, feature_table, comp_num)

% data - the bifurcation data for analysis
% srate - the sample rate
% start_or_end - if to extract onset (1) or offset (0) features
% feature_table - the initielixed feature table
% comp_num

feature_params % load label params

%% Find the start point and select the data for anlysis

zdata = zscore(data);

% start_or_end: for start use 1 for end use 0.
if ~start_or_end
    % Flip the signal for ending
    zdata = zdata(length(zdata):-1:1);
    [bif_sample, analysis_data] = find_bif_time(zdata, srate, bif_margin);
    if  bif_sample ~= -1
        bif_sample = length(zdata) - bif_sample;
    end
else
    [bif_sample, analysis_data] = find_bif_time(zdata, srate, bif_margin);
end

% Update the onset time to seconds
if bif_sample ~= -1
    bif_time = bif_sample / srate;
else
    bif_time = -1;
    bif_sample = bif_margin * srate;
end

%% extract featurs

feature_table = get_features(analysis_data, srate, feature_table, comp_num) ;

end



