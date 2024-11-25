clc
close all

%% 1. Load and Inspect Your Data
data = readtable('Seismic.csv');  % Replace with your actual file name
head(data);  % View the first few rows of the dataset
plot(data.time, data.seismic_signal_Z); % Plot the vertical (Z) channel
title('Original Seismic Signal');

%% 2. Pre-Processing: Apply Kalman Filter for Noise Reduction
% Initialize Kalman filter parameters
A = 1;  % State transition matrix
H = 1;  % Observation matrix
Q = 0.01;  % Process noise covariance
R = 0.1;  % Measurement noise covariance
P = 1;  % Initial state covariance
x_est = 0;  % Initial state estimate

% Apply Kalman filter
kalman_filtered_signal = zeros(size(data.seismic_signal_Z));

for t = 1:length(data.seismic_signal_Z)
    % Prediction step
    x_pred = A * x_est;
    P_pred = A * P * A' + Q;

    % Update step
    K = P_pred * H' / (H * P_pred * H' + R);
    x_est = x_pred + K * (data.seismic_signal_Z(t) - H * x_pred);
    P = (1 - K * H) * P_pred;

    % Store the filtered signal
    kalman_filtered_signal(t) = x_est;
end

figure;
plot(data.time, kalman_filtered_signal);
title('Kalman Filtered Signal');

%% 3. Interpolation (Cubic Spline)
interpolated_signal = interp1(data.time, kalman_filtered_signal, data.time, 'spline');

figure;
plot(data.time, kalman_filtered_signal, 'b-', data.time, interpolated_signal, 'r--');
legend('Kalman Filtered Signal', 'Interpolated Signal');
title('Kalman Filtered and Interpolated Signal');

%% 4. Empirical Mode Decomposition (EMD)
imfs = emd(interpolated_signal);
figure;
[numRows, numIMFs] = size(imfs);
for i = 1:numIMFs
    subplot(numIMFs,1,i);
    plot(data.time,imfs(:, i));
    title(['IMF ', num2str(i)]);
end

%% 5. Time-Frequency Analysis
figure;
stft(interpolated_signal, 100, 'Window', hamming(256), 'OverlapLength', 128, 'FFTLength', 512);
title('Time-Frequency Analysis (STFT)');

%% 6. Seismic Magnitude Calculation and Cumulative Energy
% Parameters for magnitude calculation
beta = 1.5;  % Constant related to energy-magnitude scaling
alpha = 4.8;  % Offset for magnitude-energy relation

% Magnitude calculation using amplitude
seismic_magnitude = beta * log10(abs(interpolated_signal)) + alpha;

figure;
plot(data.time, seismic_magnitude);
title('Seismic Magnitude Over Time');
ylabel('Magnitude');
xlabel('Time (s)');

%% 7. Estimate Seismic Moment (M₀) and Categorize Event
% Constants for seismic moment
mu = 3e10;  % Shear modulus (in Pascals)
A = 1e6;   % Fault area (in square meters)
D = 1;     % Average slip (in meters)

% Seismic moment calculation
M0 = mu * A * D;  % Seismic Moment (in N*m)

% Calculate Moment Magnitude (Mw)
Mw = (2/3) * log10(M0) - 10.7;
disp(['Moment Magnitude (Mw): ', num2str(Mw)]);

% Categorize based on Mw
if Mw < 2
    severity = 'Micro';
elseif Mw >= 2 && Mw < 3
    severity = 'Minor';
elseif Mw >= 3 && Mw < 4
    severity = 'Light';
elseif Mw >= 4 && Mw < 5
    severity = 'Moderate';
elseif Mw >= 5 && Mw < 6
    severity = 'Medium';
elseif Mw >= 6 && Mw < 7
    severity = 'Strong';
elseif Mw >= 7 && Mw < 8
    severity = 'Major';
else
    severity = 'Great';
end
disp(['Severity of Event: ', severity]);

%% 8. Plot the Seismic Energy over Time
figure;
plot(data.time, cumsum(interpolated_signal.^2));  % Cumulative energy over time
xlabel('Time (s)');
ylabel('Cumulative Seismic Energy');
title('Cumulative Seismic Energy Over Time');

%% 9. Dynamic Threshold for Eruption Detection and Display Threshold
% Set a dynamic threshold based on the mean and standard deviation
mean_signal = mean(interpolated_signal);
std_signal = std(interpolated_signal);
dynamic_threshold = mean_signal + 3 * std_signal;  % 3 standard deviations above the mean

disp(['Dynamic Threshold: ', num2str(dynamic_threshold)]);  % Display threshold

% Detect eruption events
eruption_events = interpolated_signal > dynamic_threshold;

figure;
hold on;
plot(data.time, interpolated_signal);
plot(data.time(eruption_events), interpolated_signal(eruption_events), 'r.');
legend('Interpolated Signal', 'Detected Events');
title('Dynamic Threshold-Based Volcanic Eruption Detection');
