% Define parameters
Fs = 1000;            % Sampling frequency (samples per second)
T = 5/Fs;             % Sampling period (seconds per sample)
L = 1000;             % Length of signal (number of samples)
t = (0:L-1)*T;        % Time vector (seconds)
f = 4;                % Initial frequency of the signal (Hz)
A = 1;                % Amplitude of the signal

% Generate the clean signal with constant cycle time
clean_signal = A * sin(2 * pi * f * t);

% Define noise parameters
noise_level = 0.6;    % Noise level (amplitude of noise)

% Generate noise
noise = noise_level * randn(size(t));

% Add noise to the signal with constant cycle time
noisy_signal = clean_signal + noise;

% Generate the signal with gradually increasing cycle time
gradual_cycle_signal = A * sin(2 * pi * (f + (0.5 * t)) .* t);

% Add noise to the signal with gradually increasing cycle time
noisy_gradual_signal = gradual_cycle_signal + noise;

% Compute the autocorrelation functions
[acf_clean_signal, lags] = xcorr(clean_signal, 'coeff');
[acf_noisy_signal, ~] = xcorr(noisy_signal, 'coeff');
[acf_gradual_cycle_signal, ~] = xcorr(gradual_cycle_signal, 'coeff');
[acf_noisy_gradual_signal, ~] = xcorr(noisy_gradual_signal, 'coeff');

% Plot the signals and their autocorrelations
figure('Position', [100, 100, 1200, 800]);

% Clean signal with constant cycle time
subplot(4,2,1);
plot(t, clean_signal);
title('Clean Signal with Constant Cycle Time');
xlabel('Time (s)');
ylabel('Amplitude');

subplot(4,2,2);
plot(lags/Fs, acf_clean_signal);
title('Autocorrelation of Clean Signal');
xlabel('Lags (s)');
ylabel('Autocorrelation');

% Noisy signal with constant cycle time
subplot(4,2,3);
plot(t, noisy_signal);
title('Noisy Signal with Constant Cycle Time');
xlabel('Time (s)');
ylabel('Amplitude');

subplot(4,2,4);
plot(lags/Fs, acf_noisy_signal);
title('Autocorrelation of Noisy Signal');
xlabel('Lags (s)');
ylabel('Autocorrelation');

% Clean signal with gradually increasing cycle time
subplot(4,2,5);
plot(t, gradual_cycle_signal);
title('Clean Signal with Gradually Decreasing Cycle Time');
xlabel('Time (s)');
ylabel('Amplitude');

subplot(4,2,6);
plot(lags/Fs, acf_gradual_cycle_signal);
title('Autocorrelation of Gradual Cycle Signal');
xlabel('Lags (s)');
ylabel('Autocorrelation');

% Noisy signal with gradually increasing cycle time
subplot(4,2,7);
plot(t, noisy_gradual_signal);
title('Noisy Signal with Gradually Decreasing Cycle Time');
xlabel('Time (s)');
ylabel('Amplitude');

subplot(4,2,8);
plot(lags/Fs, acf_noisy_gradual_signal);
title('Autocorrelation of Noisy Gradual Cycle Signal');
xlabel('Lags (s)');
ylabel('Autocorrelation');
 % Adjust the spacing between subplots
% Save the figure
print('signal_simulation_with_autocorrelation', '-dpng', '-r300' );
%saveas(gcf, 'signal_simulation_with_autocorrelation.png');
