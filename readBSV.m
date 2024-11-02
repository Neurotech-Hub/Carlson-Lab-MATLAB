% Open the working .bsv file for reading in binary mode
fid = fopen('/Volumes/ZEBRA/Waveform1.bsv', 'r');

% Read the file contents as 16-bit signed integers (int16)
file_data = fread(fid, 'int16');

% Close the file
fclose(fid);

% Extract the first 70 values (likely a header)
header = file_data(1:70);

% Save the header to a separate file (optional, if needed later)
save('header_data.mat', 'header');  % Save header in a .mat file

% Skip the first 70 values and get the waveform data
waveform_data = file_data(71:end);

% Plot the waveform data for analysis
plot(waveform_data);
title('Waveform Data from .BSV File (int16, Header Skipped)');
xlabel('Sample Points');
ylabel('Amplitude (16-bit signed int)');

% Display the range of the data
min_val = min(waveform_data);
max_val = max(waveform_data);
fprintf('Data range: %d to %d\n', min_val, max_val);

% Optionally, rescale the data to the [-1, 1] range for further analysis
wave_normalized = double(waveform_data) / 32768;
plot(wave_normalized);
title('Normalized Waveform Data (Header Skipped)');
