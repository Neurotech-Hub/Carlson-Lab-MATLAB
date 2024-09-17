function [R1, C1, R2, C2] = calculateRC(f0, BW, C1_value, C2_value)

% Given center frequency and bandwidth
f_L = f0 - BW/2;  % Low cutoff frequency
f_H = f0 + BW/2;  % High cutoff frequency

% Given capacitor values
C1 = C1_value;  % Set a chosen value for C1 (in Farads)
C2 = C2_value;  % Set a chosen value for C2 (in Farads)

% Calculate resistor values based on desired frequencies
R1 = 1 / (2 * pi * f_L * C1);  % Calculate R1 for high-pass
R2 = 1 / (2 * pi * f_H * C2);  % Calculate R2 for low-pass

% Display results
fprintf("\n\nRC VALUES\n");
fprintf('R1 = %.2f kOhms\n', R1 / 1000);  % Convert Ohms to kOhms
fprintf('C1 = %.2f µF\n', C1 * 1e6);      % Convert Farads to µF
fprintf('R2 = %.2f kOhms\n', R2 / 1000);  % Convert Ohms to kOhms
fprintf('C2 = %.2f µF\n', C2 * 1e6);      % Convert Farads to µF

end
