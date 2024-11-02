% Load the data
load("/Users/gaidica/Documents/MATLAB/Carlson Lab/BB1334_CT_day00_041123.eod.mat");

useEOD = 1;
wave = eod(useEOD,:);

% Plot the original waveform for verification
close all;
figure;
plot(wave);
title('Original Waveform');
xlabel('Sample Points');
ylabel('Amplitude');

% Normalize the waveform between -1 and 1
wave = wave / max(abs(wave));

% Ensure the waveform length matches a target length (e.g., 1000 points)
desired_length = 1000;  % Adjust to your device's preferred length
wave = resample(wave, desired_length, length(wave));  % Resample to the desired length

% Convert the waveform to 16-bit signed integers (range -32768 to 32767)
wave_int = int16(wave * 32767);

% Make sure the waveform is a column vector
wave_int = wave_int(:);

% Load the header data (assuming header_data.mat has been saved previously)
load('header_data.mat', 'header');

% Combine the header with the new waveform data
output_data = [header; wave_int];

% Open a binary file to save the waveform in .bsv format
filename = 'waveform_1000.bsv';
fid = fopen(filename, 'w');

% Write the header and waveform data to the .bsv file
fwrite(fid, output_data, 'int16');

% Close the file
fclose(fid);

% Optional: Display message
fprintf('Waveform with header saved to %s with %d points\n', filename, desired_length);
