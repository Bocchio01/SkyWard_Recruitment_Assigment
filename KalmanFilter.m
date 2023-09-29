clc
clear variables
close all

%% Section 0: Constants and Data Import

% Constants
P0 = 101325;       % Standard atmospheric pressure at sea level (Pa)
L = -0.0065;       % Temperature lapse rate (K/m)
T0 = 288.15;       % Standard temperature at sea level (K)
g = 9.80665;       % Acceleration due to gravity (m/s²)
M = 0.0289644;     % Molar mass of Earth's air (kg/mol)
R = 8.31446;       % Universal gas constant (J/(mol·K))

% Function to calculate altitude
% computeAltitude = @(pressure) (T0 / L) * (1 - (pressure / P0).^((R * L) / (g * M)));
computeAltitude = @(pressure) - (R * T0)/(g * M) * log(pressure / P0);

% Function to calculate acceleration magnitude
computeAccMagnitude = @(x, y, z) sqrt(x .* x + y .* y + z .* z);

% Load data from 'sim.csv'
data = readtable('sim.csv', 'NumHeaderLines', 1);

% Extract altitude, altitude variance, acceleration, and acceleration variance
alt = computeAltitude(data.Var8);
acc = computeAccMagnitude(data.Var2, data.Var3, data.Var4);

% Calculate moving average for altitude and acceleration
alt_avg = alt;
acc_avg = acc;
for i = 1:10
    alt_avg  = alt_avg (2:end) - diff(alt_avg ) / 2;
    acc_avg = acc_avg(2:end) - diff(acc_avg) / 2;
end
alt_var = ((alt(1:end-i) -alt_avg)' * (alt(1:end-i) -alt_avg)) / length(alt_avg);
acc_var = ((acc(1:end-i) - acc_avg)' * (acc(1:end-i) - acc_avg)) / length(acc_avg);

%% Section 1: Compute Kalman Matrices

% Define matrices and parameters for the Kalman filter
H = [1 0 0; 0 0 1]; % Maps state variables to sensor data
R = 2 * [alt_var 0; 0 acc_var]; % Measurement noise covariance
Q = [0 0 0; 0 0 0; 0 0 1]; % Process noise covariance matrix
T = data.Var1(2) / 1e6; % Time step
A = [1 T 1/2 * T^2; 0 1 T; 0 0 1]; % Maps previous state to next state

P = eye(3); % Initial guess for P

% Kalman filter iterations
for i = 1:20
    K = P * H' / (H * P * H' + R); % Kalman gains
    P = (eye(3) - K * H) * P;
    P = A * P * A' + Q;
end

% Initialize time vector
t = data.Var1;
estimate = zeros(3, length(t));
estimate(:, 1) = [alt(1); 0; acc(1)];

% Kalman filter estimation
for i = 2:length(t)
    estimate(:, i) = A * estimate(:, i-1);
    estimate(:, i) = estimate(:, i) + K * ([alt(i); acc(i)] - H * estimate(:, i));
end

%% Section 2: Plots

tiledlayout(2, 3);

% Plot altitude
nexttile
hold on
grid on
plot(t, alt, 'r')
plot(t, estimate(1, :), 'b')
legend('Raw', 'Estimated')
title('Altitude')
xlabel('t [us]');
ylabel('Altitude [m.s.l.]');
hold off

% Plot velocity
nexttile
hold on
grid on
plot(t, estimate(2, :), 'b')
legend('Estimated')
title('Velocity')
xlabel('t [us]');
ylabel('Velocity [m/s]');
hold off

% Plot acceleration
nexttile
hold on
grid on
plot(t, acc, 'r')
plot(t, estimate(3, :), 'b')
title('Acceleration')
xlabel('t [us]');
ylabel('Acceleration [m/s²]');
legend('Raw', 'Estimated')
hold off

% Plot altitude gain
nexttile
hold on
grid on
plot(t(2:end), diff(estimate(1, :)), 'b')
legend('Altitude gain')
title('Altitude Gain')
xlabel('t [us]');
ylabel('Gain [m]');
hold off

% Plot acceleration gain
nexttile
hold on
grid on
plot(t(2:end), diff(estimate(3, :)), 'b')
legend('Acceleration gain')
title('Acceleration Gain')
xlabel('t [us]');
ylabel('Gain [m/s²]');
hold off


