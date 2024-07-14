function [start_time, mean_autocorr, normalized_std_autocorr,...
          entropy, entropy_autocorr, stds] = ...
          get_bif_time_and_features(zdata, srate, w, start_or_end, bif_margin)

% zdata: the data of the bifurcation possible range
% srate: the sample rate 
% w: the windoe size for features around the changepoint.
% start_or_end: for start use 1 for end use 0. 
zdata = zscore(zdata);

if ~start_or_end
% if end bifurcation revers the signal for analysis
    zdata = zdata(length(zdata):-1:1);
end
start_data = zdata(1:srate*bif_margin*2);

TF = ischange(start_data,'variance', 'MaxNumChanges',1);

if sum(TF)== 0
    TF = ischange(start_data,'linear', 'MaxNumChanges',1);
end

if sum(TF)== 0
    TF(srate*bif_margin) = 1;
end

% use the changepoints as start points to evaluate features then
% select the best

w_size = w*srate;
start = find(TF); % ceil(1:step:length(zdata));
start_time = start/srate;

if ~start_or_end
% if end bifurcation the revers time is acually the start
    start_time = length(zdata)/srate - start_time ;
end

win = start:min(start+w_size, length(zdata));
data = zdata(win);
[autocor,~] = xcorr(data, length(win), 'coeff');
[~,lcsh] = findpeaks(autocor,"MinPeakProminence",0.1);

mean_autocorr = mean(diff(lcsh))/srate;
normalized_std_autocorr = std(diff(lcsh)/(srate))/mean_autocorr;
entropy = wentropy(data,'shannon');
entropy_autocorr = wentropy(autocor,'shannon');
stds =  std(data(1:min(srate, length(data))));


% high entropy autocorrelation > 100
% std autocorrelation relative to mean - low close to 0 <0.1
%% minimal negative entropy close to 0 >-1000
% mean autocorrelation distance ~ ocilation cycle time mostly (0.1-0.75)
% stds range are related to the signals amplitude (local) > 1
% envelope > 3
% the first index to satisfy the conditions is selected

% option 2: iregular occilassion
%% low entropy (close to -4000) + high stds/envelope


