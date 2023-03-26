clear all

%% Homework 2 Problem 1
% Author: Darrell Hamilton
% Date: 03/10/2023
% Description: This script was created as a "calculator" for Homework Problem 1.
% Calculates optimal values of a, e, and i to minimize omega_dot for an 
% Molniya orbit around Earth

% Initialize Constants
J_2 = 0.00108;          % Earth J_2
R = 6370;               % Radius of Earth
mu = 3.986*10^5;        % Gravitational constant for earth

% Calculate Period
T = (1/3)*3600*24;

% Calculate Semimajor Axis
a = ((T/(2*pi))^2*mu)^(1/3);

% Calculate mean motion
n = sqrt(mu/a^3);

% Determine maximum eccentricity due to altitude constraint
e_max = -(600 + R)/a + 1;

% Calculate inclination
i = acos(sqrt(1/5))*180/pi;

% Plot omega_dot to
e=0.0:0.01:0.99;
R_alt = a*(1-e) - R;
idx = find(R_alt > 600);
e = e(idx);

omega_dot = -3/2*n*J_2*(R/a)^2*cos(i)./(1-e.^2).^2;
plot(e,omega_dot)
title("Earth Molniya Orbit")
ylabel("Longitude of Ascending Node [rad]")
xlabel("e")



