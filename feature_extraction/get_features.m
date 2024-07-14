% get features from a given data channel

function feature_table = get_features(data, srate, feature_table, comp_num)

feature_params
n_samples = length(data);

% Statistical features + get the most promoinant peaks pointing up
skewness_data = skewness(data);
flip = sign(skewness_data);
if flip == 0
    flip = 1;
end

feature_table(comp_num, :).skewness = skewness_data*flip;
data = data*flip;
feature_table(comp_num, :).mean = mean(data);
feature_table(comp_num, :).median = median(data);


% features  related to critical slowing
try
    [width, height] = fwhm(data, srate);
catch
    width = -1;
    height = -1;
end

feature_table(comp_num, :).width_autocorr = width;
feature_table(comp_num, :).height_autocorr = height;


% Autocorrealtion features + the entropy
try
    [autocor,~] = xcorr(data, n_samples, 'coeff');
    [~,lcsh] = findpeaks(autocor,"MinPeakProminence",0.1);
    mean_autocorr = mean(diff(lcsh))/srate;
    feature_table(comp_num, :).mean_autocorr = mean_autocorr;
    feature_table(comp_num,:).std_autocorr = std(diff(lcsh)/(srate));
    feature_table(comp_num,:).entropy_autocorr = wentropy(autocor, 'shannon');
catch
    fprintf('autocorr fail')
end


% Spectral slope and entropy
if length(data) > 2*srate
    [pxx, f] = pwelch(data, srate, 10,1:0.01:40,srate);
    try
        [slope, slope_low, slope_high, intercept, intercept_low, intercept_high, rmse, adjrsquare, sse, dfe] =  fit_bif(f',log(pxx)','linear');
        [~, ind] = max(log(pxx) - (slope*f + intercept));
        feature_table(comp_num, :).freq = f(ind);

        feature_table(comp_num, :).spectral_slope = slope;
        feature_table(comp_num, :).spectral_slope_high = slope_high;
        feature_table(comp_num, :).spectral_slope_low = slope_low;
        feature_table(comp_num, :).spectral_intercept = intercept;
        feature_table(comp_num, :).spectral_intercept_high = intercept_high;
        feature_table(comp_num, :).spectral_intercept_low = intercept_low;
        feature_table(comp_num, :).spectral_rmse = rmse;
        feature_table(comp_num, :).spectral_rsquare = adjrsquare;
        feature_table(comp_num, :).spectral_sse = sse;
        feature_table(comp_num, :).spectral_dfe = dfe;

    catch
        fprintf('spectral fit failed')
    end

    try
        feature_table(comp_num, :).spectral_entropy = wentropy(pxx,'shannon');
    catch
        fprintf('spectral entropy failed')
    end
end

% The time series aproximate entropy
feature_table(comp_num, :).approx_entropy = approximateEntropy(data);

% Amplitude evaluation using rms envelope
rms_samples = time_from_bif_rms*srate;
feature_table(comp_num, :).stds = std(data(1:min(rms_samples, length(data))));

if rms_samples <= length(data)
    try
        [yupper,ylower] = envelope(data(1:rms_samples), rms_samples,'rms');
        differance = yupper - ylower;

        x = linspace(1,time_from_bif_rms+1 ,rms_samples);
        [a, a_low, a_high, b, b_low, b_high, rmse, adjrsquare, sse, dfe] = fit_bif(x', differance', 'linear');
        feature_table(comp_num, :).rms_linear_slope = a;
        feature_table(comp_num, :).rms_linear_slope_low = a_low;
        feature_table(comp_num, :).rms_linear_slope_high = a_high;
        feature_table(comp_num, :).rms_linear_intercept = b;
        feature_table(comp_num, :).rms_linear_intercept_low = b_low;
        feature_table(comp_num, :).rms_linear_intercept_high = b_high;
        feature_table(comp_num, :).rms_linear_rmse = rmse;
        feature_table(comp_num, :).rms_linear_adjrsquare  = adjrsquare;
        feature_table(comp_num, :).rms_linear_sse = sse;
        feature_table(comp_num, :).rms_linear_dfe  = dfe;

        [a, a_low, a_high, b, b_low, b_high, rmse, adjrsquare, sse, dfe] = fit_bif(x', differance', 'suph');
        feature_table(comp_num, :).rms_suph_stretch = a;
        feature_table(comp_num, :).rms_suph_stretch_low = a_low;
        feature_table(comp_num, :).rms_suph_stretch_high = a_high;
        feature_table(comp_num, :).rms_suph_intercept = b;
        feature_table(comp_num, :).rms_suph_intercept_low = b_low;
        feature_table(comp_num, :).rms_suph_intercept_high = b_high;

        feature_table(comp_num, :).rms_suph_adjrsquare  = adjrsquare;
        feature_table(comp_num, :).rms_suph_rmse  = rmse;
        feature_table(comp_num, :).rms_suph_sse  = sse;
        feature_table(comp_num, :).rms_suph_dfe  = dfe;
    catch
        fprintf('rms fit failed')
    end
end

% Peaks for ISI fit
[pks,locs] = findpeaks(data, srate,  ...
    "MinPeakProminence",1, ...
    "MinPeakDistance",min_peak_distance);

if length(pks) >= 5
    try
        isi = diff(locs);
        feature_table(comp_num, :).isi_mean = mean(isi);
        feature_table(comp_num, :).isi_std = std(isi);

        feature_table(comp_num, :).peak_mean = mean(pks);
        feature_table(comp_num, :).peak_std = std(pks);


        x = linspace(locs(1), locs(end), length(isi))';
        isi = isi';
        [a, a_low, a_high, b, b_low, b_high, rmse, adjrsquare, sse, dfe] = fit_bif(x,isi, 'linear');
        feature_table(comp_num, :).isi_linear_slope = a;
        feature_table(comp_num, :).isi_linear_slope_low = a_low;
        feature_table(comp_num, :).isi_linear_slope_high = a_high;
        feature_table(comp_num, :).isi_linear_intercept = b;
        feature_table(comp_num, :).isi_linear_intercept_low = b_low;
        feature_table(comp_num, :).isi_linear_intercept_high = b_high;
        feature_table(comp_num, :).isi_linear_rmse = rmse;
        feature_table(comp_num, :).isi_linear_adjrsquare = adjrsquare;
        feature_table(comp_num, :).isi_linear_sse = sse;
        feature_table(comp_num, :).isi_linear_dfe = dfe;

        [a, a_low, a_high, b, b_low, b_high, rmse, adjrsquare, sse, dfe] = fit_bif(x,isi, 'snic');
        feature_table(comp_num, :).isi_snic_stretch = a;
        feature_table(comp_num, :).isi_snic_stretch_low = a_low;
        feature_table(comp_num, :).isi_snic_stretch_high = a_high;
        feature_table(comp_num, :).isi_snic_intercept = b;
        feature_table(comp_num, :).isi_snic_intercept_low = b_low;
        feature_table(comp_num, :).isi_snic_intercept_high = b_high;
        feature_table(comp_num, :).isi_snic_adjrsquare  = adjrsquare;
        feature_table(comp_num, :).isi_snic_rmse  = rmse;
        feature_table(comp_num, :).isi_snic_sse  = sse;
        feature_table(comp_num, :).isi_snic_dfe  = dfe;

        [a, a_low, a_high, b, b_low, b_high, rmse, adjrsquare, sse, dfe] = fit_bif(x,isi, 'sh');
        feature_table(comp_num, :).isi_sh_stretch = a;
        feature_table(comp_num, :).isi_sh_stretch_low = a_low;
        feature_table(comp_num, :).isi_sh_stretch_high = a_high;
        feature_table(comp_num, :).isi_sh_intercept = b;
        feature_table(comp_num, :).isi_sh_intercept_low = b_low;
        feature_table(comp_num, :).isi_sh_intercept_high = b_high;
        feature_table(comp_num, :).isi_sh_adjrsquare  = adjrsquare;
        feature_table(comp_num, :).isi_sh_rmse  = rmse;
        feature_table(comp_num, :).isi_sh_sse  = sse;
        feature_table(comp_num, :).isi_sh_dfe  = dfe;
    catch
        fprintf('isi fit failed')
    end
end


% Peaks for amplitude change fit
% use only the locs included for the peak hight analysis time
pks = pks(locs<=time_from_bif_peak);
locs = locs(locs<=time_from_bif_peak);

% try to fit only with over 5 peaks not less

pks = pks';
if length(pks) >= 5
    try
        x = (linspace(locs(1), locs(end), length(pks)))';
        [a, a_low, a_high, b, b_low, b_high, rmse, adjrsquare, sse, dfe] = fit_bif(x, pks, 'linear');
        feature_table(comp_num, :).peak_linear_slope = a;
        feature_table(comp_num, :).peak_linear_slope_low = a_low;
        feature_table(comp_num, :).peak_linear_slope_high = a_high;
        feature_table(comp_num, :).peak_linear_intercept = b;
        feature_table(comp_num, :).peak_linear_intercept_low = b_low;
        feature_table(comp_num, :).peak_linear_intercept_high = b_high;
        feature_table(comp_num, :).peak_linear_adjrsquare  = adjrsquare;
        feature_table(comp_num, :).peak_linear_rmse  = rmse;
        feature_table(comp_num, :).peak_linear_sse  = sse;
        feature_table(comp_num, :).peak_linear_dfe  = dfe;

        [a, a_low, a_high, b, b_low, b_high, rmse, adjrsquare, sse, dfe] = fit_bif(x, pks, 'suph');
        feature_table(comp_num, :).peak_suph_stretch = a;
        feature_table(comp_num, :).peak_suph_stretch_low = a_low;
        feature_table(comp_num, :).peak_suph_stretch_high = a_high;
        feature_table(comp_num, :).peak_suph_intercept = b;
        feature_table(comp_num, :).peak_suph_intercept_low = b_low;
        feature_table(comp_num, :).peak_suph_intercept_high = b_high;
        feature_table(comp_num, :).peak_suph_adjrsquare  = adjrsquare;
        feature_table(comp_num, :).peak_suph_rmse  = rmse;
        feature_table(comp_num, :).peak_suph_sse  = sse;
        feature_table(comp_num, :).peak_suph_dfe  = dfe;
    catch
        fprintf('amplitude peak fit failed')
    end
end
end




function [a, a_low, a_high, b, b_low, b_high, rmse, adjrsquare, sse, dfe] = fit_bif(x, y, name)

% evaluate for morphological ternd for typical bifurcation
% assumes the end signal is flipped
start_point = [0, y(1)];
lower  = [0,  -inf ];

switch name
    case 'snic' %on ISI
        ft = fittype( 'b + a/sqrt(x)');
    case 'sh' % on ISI
        ft = fittype( 'b + a* log(x/min(x))');
    case 'suph' %on peak amplitude /rms
        ft = fittype( 'b + a*sqrt(x)');
    case 'linear'
        ft = fittype('b + a*x');
        lower  = [-inf, -inf];
        start_point = [0, y(1)];
    case 'log_linear'
        ft = fittype('b + a*x');
        lower  = [-inf, -inf];
        y = log(y);
end

% parametric time scale, the scaling law relates to
% peak dynamics and not the absolute time scale so time points are ordered
options = fitoptions(ft);
options.Lower = lower;
options.StartPoint = start_point;
[fitobject,gof] = fit(x, y, ft,options);

param_ci = confint(fitobject, 0.95);

a = param_ci (:,1);
a_high = a(1);
a_low = a(2);
a = fitobject.a;

b = param_ci (:,2);
b_high = b(1);
b_low = b(2);
b = fitobject.b;

adjrsquare = gof.adjrsquare;
rmse = gof.rmse;
sse = gof.sse;
dfe = gof.dfe;

end