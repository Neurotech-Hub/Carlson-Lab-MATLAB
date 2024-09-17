load('/Users/mattgaidica/Library/CloudStorage/Box-Box/EODs for Matt/070824_CN12_EOD_70ms.eod.mat');
doSave = true;

testNoise = true;
if (testNoise)
    % high amplitude artifact
    eod(1,6000:6050) = 60;
    
    % Define noise amplitude (adjust this value based on desired noise level)
    noiseAmplitude = 5;  % Example: 5% of the signal amplitude
    
    % Add random noise to each signal in eod
    eod(2,:) = eod(2,:) + noiseAmplitude * randn(size(eod(2,:)));
end

% Assuming eod is a 10x13672 matrix and srate is the sampling rate
[numSignals, L] = size(eod);
t = (0:L-1)/srate;  % Time vector
f = (0:L-1)*(srate/L);  % Frequency vector

% Define the window (e.g., Hamming window)
window = hamming(L)';

% Smoothing function (moving average with window size N)
smoothN = 5;  % Adjust window size for smoothing
smoothFFT = @(P) movmean(P, smoothN);

% Initialize array to store peak frequencies
peakFrequencies = zeros(1, numSignals);

% Create a figure for both subplots
rows = 2;
cols = 2;
close all;
figure('Position', [0, 0, 1200, 1200]);

% Create subplot for time-domain signals
subplot(rows, cols, 1);  % Left subplot
hold on;
for i = 1:numSignals
    plot(t, eod(i, :));
end
title('Time-Domain Signals');
xlabel('Time (seconds)');
ylabel('Amplitude');
grid on;
hold off;

% Create subplot for frequency-domain (FFT) signals
subplot(rows, cols, 2);  % Right subplot
hold on;
for i = 1:numSignals
    % Apply windowing to the signal
    windowedSignal = eod(i, :) .* window;
    
    % Compute the FFT of the windowed signal
    Y = fft(windowedSignal);
    P2 = abs(Y/L);  % Two-sided spectrum
    P1 = P2(1:L/2+1);  % Single-sided spectrum
    P1(2:end-1) = 2*P1(2:end-1);  % Correct amplitude
    
    % Apply smoothing to the FFT magnitude
    P1_smooth = smoothFFT(P1);
    
    % Find the peak frequency
    [~, peakIndex] = max(P1_smooth);  % Get index of max value in P1_smooth
    peakFrequencies(i) = f(peakIndex);  % Store the corresponding frequency
    
    % Plot the smoothed FFT magnitude for each signal
    plot(f(1:L/2+1), P1_smooth);
end
% Calculate the average peak frequency
averagePeak = mean(peakFrequencies);

title(sprintf('Smoothed Spectrums, Peak %.2f Hz',averagePeak));
xlabel('Frequency (Hz)');
ylabel('|P1(f)|');
xlim([1, 200]);  % Adjust frequency range as needed
grid on;
xline(averagePeak,'r:','linewidth',2);
hold off;


% Example input parameters
f0 = averagePeak;       % Center frequency (Hz)
BW = 40;        % Bandwidth (Hz)
C1_value = 1e-8; % Chosen C1 value (1 uF)
C2_value = 1e-8; % Chosen C2 value (1 uF)

% Calculate R1, C1, R2, C2
[R1, C1, R2, C2] = calculateRC(f0, BW, C1_value, C2_value);

% Define the transfer function numerator and denominator
num = [R2*C2 0];  % s*R2*C2 (numerator of transfer function)
den = [R1*R2*C1*C2 (R1*C1 + R2*C2) 1];  % (sR2C2 + 1)(sR1C1 + 1)

% Sampling rate (from your data)
Ts = 1/srate;  % Sampling period

% Convert the analog filter to digital using bilinear transformation
[b, a] = bilinear(num, den, 1/Ts);  % Digital filter coefficients

% Initialize figure for filtered signals
% close all;
% figure('Position', [0, 0, 1200, 600]);

% Number of signals and their length
[numSignals, L] = size(eod);
t = (0:L-1)/srate;  % Time vector

% Subplot for time-domain signals
subplot(rows, cols, 3);  
hold on;
for i = 1:numSignals
    % Filter each signal using the digital filter
    filtered_signal = filter(b, a, eod(i, :));
    
    % Plot filtered signal
    plot(t, filtered_signal);
end
title('Filtered Signals in Time Domain');
xlabel('Time (seconds)');
ylabel('Amplitude');
grid on;
hold off;

% Subplot for FFT of filtered signals
subplot(rows, cols, 4);  
hold on;
for i = 1:numSignals
    % Compute the FFT of the filtered signal
    filtered_signal = filter(b, a, eod(i, :));
    Y = fft(filtered_signal);
    
    % Compute the magnitude of the FFT
    P2 = abs(Y/L);  
    P1 = P2(1:L/2+1);  
    P1(2:end-1) = 2*P1(2:end-1);
    
    % Smooth the FFT magnitude using a moving average
    P1_smooth = movmean(P1, smoothN);
    
    % Frequency vector
    f = (0:L-1)*(srate/L);  
    
    % Plot smoothed FFT magnitude
    plot(f(1:L/2+1), P1_smooth);
end
title('Spectrum of Filtered Signals');
xlabel('Frequency (Hz)');
ylabel('|P1(f)|');
xlim([1, 200]);  % Adjust frequency range as needed
grid on;
hold off;
set(gcf,'color','w');

if (doSave)
    exportgraphics(gcf,'RCfilter.jpg');
end
