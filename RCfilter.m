% Directory containing the .mat files
directory = '/Users/gaidica/Library/CloudStorage/Box-Box/EODs for Matt';

% Get a list of all .mat files in the directory
files = dir(fullfile(directory, '*.mat'));

% Initialize a cell array to store all eod matrices temporarily
eod_signals = {};

% Loop over each file to load eod signals
for k = 1:length(files)
    % Construct the full file path
    filepath = fullfile(directory, files(k).name);
    
    % Load the .mat file
    data = load(filepath);
    
    % Check if the 'eod' variable exists in the loaded data
    if isfield(data, 'eod')
        eod_signals{end+1} = data.eod;  % Append eod matrix to the cell array
    else
        warning('Variable ''eod'' not found in %s', files(k).name);
    end
end

% Determine the number of channels (should be 10)
num_channels = size(eod_signals{1}, 1);

% Find the maximum length (number of time points) across all eod signals
max_length = max(cellfun(@(x) size(x, 2), eod_signals));

% Initialize a matrix to hold all centered signals, padded with zeros
compiled_eod = zeros(num_channels * length(eod_signals), max_length);

% Index to keep track of where to place each 10-channel block
row_idx = 1;

% Center all signals within the maximum length (zero-pad both sides)
for k = 1:length(eod_signals)
    current_eod = eod_signals{k};
    current_length = size(current_eod, 2);  % Length of the current signal
    
    % Calculate padding needed to center the signal
    padding_left = floor((max_length - current_length) / 2);
    padding_right = max_length - current_length - padding_left;
    
    % Place the centered signal in the appropriate rows and columns
    compiled_eod(row_idx:row_idx+num_channels-1, (1+padding_left):(max_length-padding_right)) = current_eod;
    
    % Update row index for the next set of signals
    row_idx = row_idx + num_channels;
end

% Display the size of the compiled eod data
fprintf('Compiled eod size: %dx%d\n', size(compiled_eod));

% Save the compiled eod data to a new .mat file
save('compiled_eod_centered.mat', 'compiled_eod');

eod = compiled_eod; % overwrite

%%
doSave = true;
testNoise = false;

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
    plot(t, normalize(eod(i, :)));
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
xlim([1, 20000]);  % Adjust frequency range as needed
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

% Define the transfer function numerator and denominator for your specific circuit
num = [C2 * R2 0];  % Numerator: C2 * R2 * s (s term for Laplace transform)
den = [R1 * C1 * C2 * R2, (R1 * C1 + C2 * R2), 1];  % Denominator: (1 + R1C1s)(1 + C2R2s)

% Sampling rate
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
xlim([1, 20000]);  % Adjust frequency range as needed
grid on;
hold off;

if (doSave)
    exportgraphics(gcf,'RCfilter.jpg');
end
