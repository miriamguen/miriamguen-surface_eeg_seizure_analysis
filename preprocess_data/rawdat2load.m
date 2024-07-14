function [signals,Fs,chan_names,locs] = rawdat2load(fname)
% this function reads the header and binary files:
% input: headr file name 
% output: EEG matrix, Sample Rate  


% Read in header info 
fid = fopen(fname);
header = textscan(fid, '%s%[^\n]');
[name,header] = strtok(header{1}(2:end),'=');
[~, ind] = intersect(name,{'conversion_factor','elec_names',...
    'num_channels','num_samples','sample_freq'});

factor      =  str2double(header{ind(1),:}(2:end));
chan_names  =  split(header{ind(2),:}(3:end-1),',');
channels    =  str2double(header{ind(3),:}(2:end));
samples     =  str2double(header{ind(4),:}(2:end));
Fs          =  str2double(header{ind(5),:}(2:end));

fclose(fid); 
clear fid;


% Read in actual EEG data from the binary file 
fname = strrep(fname,'.head','.data');
file = fopen(fname,'rb');              
signals = fread(file, [channels samples], 'int16');    
fclose(file); clear file;      

% scale raw data by conversion factor
if samples ~= 0
    signals = signals.*factor;   
    eeglab nogui
    [chan_names, ord, locs] = elc2ord(chan_names);%(strrep(fname,'.data','.head'));
    signals = signals(ord,:);
end

end

function [chan_names, ord, locs] = elc2ord(chan_names)
%read standart locs to select electrodes with known locations
locs = readlocs('standard-10-5-cap385_no_eog.elp','filetype' ,'besa'); 
stan_chan_names = upper({locs.labels});
chan_names = upper(chan_names);
[~,ord2,ord] = intersect(stan_chan_names,chan_names);
chan_names = chan_names(ord);
locs = locs(ord2);
end
