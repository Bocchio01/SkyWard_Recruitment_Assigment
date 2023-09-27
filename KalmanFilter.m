clc
clear variables
close all

%% Section 0


% Constants
P0 = 101325;       % Standard atmospheric pressure at sea level (Pa)
L = -0.0065;       % Temperature lapse rate (K/m)
T0 = 288.15;       % Standard temperature at sea level (K)
g = 9.80665;       % Acceleration due to gravity (m/s²)
M = 0.0289644;     % Molar mass of Earth's air (kg/mol)
R = 8.31446;       % Universal gas constant (J/(mol·K))

% Calculate altitude
altitude = @(pressure) (T0 / L) * (1 - (pressure / P0).^((R * L) / (g * M)));
acceleration = @(x, y, z) sqrt(x .* x + y .* y + z .* z);

data = readtable('sim.csv', 'NumHeaderLines', 1);

alt = altitude(data.Var8);
accel = acceleration(data.Var2, data.Var3, data.Var4);

%% Section 1: compute Kalman matrixes

H = [1 0 0; 0 0 1]; % maps x (state variables) to z (sensor data)
R = [100 0; 0 100]; % measurement noise covariance
Q = [0 0 0; 0 0 0; 0 0 1]; % process noise covariance matrix
T = data.Var1(2) / 1e6; % time step
% T = 0.02;
A = [1 T 1/2 * T^2; 0 1 T; 0 0 1]; % maps previous state to next state

P = eye(3); % initial guess for p
for i = 1:20
    K = P*H'/(H*P*H' + R); % Kalman gains
    P = (eye(3) - K *H)*P;
    P = A*P*A' + Q;
end

% display(K)
% display(H)
% display(P)

t = data.Var1;
estimate = zeros(3,length(t));
estimate(:,1) = [alt(1); 0; accel(1)];
for i = 2:length(t)
    estimate(:,i) = A*estimate(:,i-1);
    estimate(:,i) = estimate(:,i) + K*([alt(i); accel(i)] - H *estimate(:,i));
end

%% Section 2: plots

nexttile
hold on
title('Altitude')
plot(data.Var1, estimate(1,:))
plot(data.Var1, alt, 'r')
hold off

nexttile
title('Velocity')
plot(data.Var1, estimate(2,:))

nexttile
hold on
title('Acceleration')
plot(data.Var1, estimate(3,:))
plot(data.Var1, accel, 'r')
legend('Filtered', 'Raw')












