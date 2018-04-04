clear all; clc; close all;

u0     = 235.9;
w0     = 0;
q0     = 0;
theta0 = 0;
time   = [0, 1000];
initial_state = [u0, w0, q0, theta0];

% delta u = 10
initial_delta = [0, 0, 0, 0];
[t, y] = ode45(@(t, y) table4_4(t, y, initial_state), time, initial_delta);

figure; hold on; grid on;
title('Trim Condition');
plot(t, y(:, 1), 'LineWidth', 2, 'DisplayName', '\Delta u^E')
plot(t, y(:, 2), 'LineWidth', 2, 'DisplayName', '\Delta w^E')
plot(t, y(:, 3), 'LineWidth', 2, 'DisplayName', '\Delta q')
plot(t, y(:, 4), 'LineWidth', 2, 'DisplayName', '\Delta \theta')
legend('show')
print('trim', '-dpng')
% delta u = 10
initial_delta = [10, 0, 0, 0];
[t, y] = ode45(@(t, y) table4_4(t, y, initial_state), time, initial_delta);

figure; hold on; grid on;
title('\Delta u^E = 10 m/s');
plot(t, y(:, 1), 'LineWidth', 2, 'DisplayName', '\Delta u^E')
plot(t, y(:, 2), 'LineWidth', 2, 'DisplayName', '\Delta w^E')
plot(t, y(:, 3), 'LineWidth', 2, 'DisplayName', '\Delta q')
plot(t, y(:, 4), 'LineWidth', 2, 'DisplayName', '\Delta \theta')
legend('show')
print('deltau', '-dpng')

% delta w = 10
initial_delta = [0, 10, 0, 0];
[t, y] = ode45(@(t, y) table4_4(t, y, initial_state), time, initial_delta);

figure; hold on; grid on;
title('\Delta w^E = 10 m/s');
plot(t, y(:, 1), 'LineWidth', 2, 'DisplayName', '\Delta u^E')
plot(t, y(:, 2), 'LineWidth', 2, 'DisplayName', '\Delta w^E')
plot(t, y(:, 3), 'LineWidth', 2, 'DisplayName', '\Delta q')
plot(t, y(:, 4), 'LineWidth', 2, 'DisplayName', '\Delta \theta')
legend('show')
print('deltaw', '-dpng')

% delta q = 0.1 s^-1
initial_delta = [0, 0, 0.1, 0];
[t, y] = ode45(@(t, y) table4_4(t, y, initial_state), time, initial_delta);

figure; hold on; grid on;
title('\Delta q = 0.1 rad/s');
plot(t, y(:, 1), 'LineWidth', 2, 'DisplayName', '\Delta u^E')
plot(t, y(:, 2), 'LineWidth', 2, 'DisplayName', '\Delta w^E')
plot(t, y(:, 3), 'LineWidth', 2, 'DisplayName', '\Delta q')
plot(t, y(:, 4), 'LineWidth', 2, 'DisplayName', '\Delta \theta')
legend('show')
print('deltaq', '-dpng')

% delta theta = 0.1 
initial_delta = [0, 0, 0, 0.1];
[t, y] = ode45(@(t, y) table4_4(t, y, initial_state), time, initial_delta);

figure; hold on; grid on;
title('\Delta \theta = 0.1 rad/s');
plot(t, y(:, 1), 'LineWidth', 2, 'DisplayName', '\Delta u^E')
plot(t, y(:, 2), 'LineWidth', 2, 'DisplayName', '\Delta w^E')
plot(t, y(:, 3), 'LineWidth', 2, 'DisplayName', '\Delta q')
plot(t, y(:, 4), 'LineWidth', 2, 'DisplayName', '\Delta \theta')
legend('show')
print('deltatheta', '-dpng')

