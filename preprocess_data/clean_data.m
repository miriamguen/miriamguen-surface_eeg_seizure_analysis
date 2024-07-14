function [clean_data,w,sphr] = Clean_data_1(data,srate,n_comp,bnd,dsmp,icawin,reref, chan_names,locs)
%PIPELINE using EEGLAB functions to pre-process data. ICA component
%rejection is supervised
%
%INPUTS 1)DATA: channels x samples (matlab double) | eeglab structure
%       2)SRATE: sample rate in Hz
%       3)n_comp: number of pca componants to use
%       with channels to be removed for eeglab structs (e.g. ecg channels; [] for none)
%       4)BAND (temporal filtration): e.g. [1 40] or [1 0] (no low-pass)
%       5)DSAMP: a number for downsampling (Hz)
%
%       6)ICAWIN allows to compute ICs only on an intial segment of the
%       data. n (samples) or n with s|sec|sec.|second
%       7)REREF: rereferencinf for EEG only (default = 1)
%       8)chanel names
%       9)channel locations for eeglab
%
%OUTPUTS CLEAN_DATA a channel x sample matlab matrix. If an EEG lab format
%output is desired use:
%   [clean_data,~]: returnes clean data in EEGLAB struct with the ic labels
%
%   [clean_data,w,sphr]:
%          clean_data: as above
%          w = the ica whight matrix
%          sphr = the sphrical data




%% check input

%only for EEG
if ~exist('reref','var') || isempty(reref)
    reref = 1;
end


if ~exist('icawin','var') || isempty(icawin)|| icawin == 0
    icawin = size(data,2);
else
    if isnumeric(icawin)
        icawin = min(icawin,size(data,2));
    else
        error('Please define ICA window length as a numeric value n seconds')
    end
end


%% add EEGLab to path if it is not in it
if isempty(which('eegplot'))
    eeglab
end
eeg_options;
if option_single
    pop_editoptions('option_single',0,'option_computeica',0);
end

%% start preproc

%import data to eeglab data format
EEG = pop_importdata('dataformat','matlab','nbchan',size(data,1),...
    'data',data,'srate',srate,'pnts',size(data,2),'xmin',0);
EEG = eeg_checkset(EEG);

%import channel positions
EEG.chanlocs = locs;
EEG = eeg_checkset(EEG);
%downsample data
EEG = pop_resample(EEG, dsmp);

temp_locs = EEG.chanlocs;
%high pass (DONE IN CLEAN ARTIFACTS)
EEG = pop_eegfiltnew(EEG, bnd(1),0, 1690, 0,[], 0);
%EEG = cleanline(EEG, 'LineFrequencies',[25 50 100], 'LineAlpha',0.05,...
%    'Bandwidth',2, 'SignalType', 'Channels', 'VerboseOutput', 'false');

% Apply clean_artifacts() to reject bad channels and correct continuous data using Artifact Subspace Reconstruction (ASR)
[EEG,HP,BUR] = clean_artifacts(EEG, 'BurstCriterion', 10,...
    'BurstRejection','off','BurstCriterionRefMaxBadChns',0.3,...
    'ChannelCriterion',0.5, 'WindowCriterion', 'off');

%lowpass filter:
EEG = pop_eegfiltnew(EEG, 0, bnd(2), 1690, 0,[], 0); % Lowpass FIR filter

EEG = pop_interp(EEG, temp_locs); %interpulate for consistancy

if reref
    %Re-reference the data to average
    EEG.nbchan = EEG.nbchan+1;
    EEG.data(end+1,:,:) = zeros(1, EEG.pnts, size(EEG.data,3));
    EEG.chanlocs(1,EEG.nbchan).labels = 'initialReference';
    EEG = pop_reref(EEG, []);
    EEG = pop_select( EEG,'nochannel',{'initialReference'});
end


if icawin < size(EEG.data,2)
    full_dat = EEG.data;
    EEG.data = EEG.data(:,1:icawin);
    EEG.pnts = icawin;
    xmax = EEG.xmax;
    times = EEG.times;
end

% run ICA with parameters ajusted to seizure information with our approch
EEG = pop_runica(EEG,'extended',round(n_comp/2),'interupt','off', ...
                  'lrate',0.0049,'momentum',0.3,'stop',1e-9, 'pca',n_comp);

if exist('full_dat','var')
    EEG.data = full_dat;
    EEG.pnts = size(full_dat,2);
    EEG.xmax = xmax;
    EEG.times = times;
end

% return the ica wheigets and spher values
if nargout >1
    w = EEG.icaweights;
    sphr = EEG.icasphere;
end

% automatic IC labeling using IClabel

if ~isreal(EEG.icaweights)
    EEG.icarjct = 'fail';
    EEG.corrupt = 'complex';
else
    try
        EEG = pop_iclabel(EEG, 'default');
        if EEG.xmax > 20
            [mara.artcomps, mara.info] = MARA(EEG);
            EEG.etc.ic_classification.MARA = mara;
        end
    catch
        disp('Automatic IC labeling failed');
        EEG.icarjct = 'fail';
    end
   
    disp(size(data,2)/srate)
    
    % use this to visualy evaluate the outputs
    %pop_viewprops(EEG, 0, [1:n_comp],{'freqrange', bnd}, {}, 1, 'ICLabel' );
    %pop_eegplot( EEG, 1, 0, 1);
    %pop_eegplot( BUR, 1, 0, 1);
    %pop_eegplot( EEG, 0, 0, 1);
    
    clean_data = EEG;
    
end

