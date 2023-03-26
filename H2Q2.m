clear all
%% Homework 2 Problem 2
% Author: Darrell Hamilton
% Date: 03/10/2023
% Description: This script was created as a "calculator" for Homework Problem 1.
% Calculates optimal values of a, e, and i to minimize omega_dot for an 
% Molniya orbit around Mars


% Initialize Constants
J_2 = 0.00196;      % J_2 for Mars
R = 3390;           % Martian radius (km)
mu = 4.282*10^4;    % Gravitational constant for Mars (km^3/s^2)

% Calculate Period
T = 3600*24 + 39*60 + 35;

% Calculate semimajor axis
a = ((T/(2*pi))^2*mu)^(1/3);
% Calculate mean motion
n = sqrt(mu/a^3);

% Calculate and plot omega_dot
e=0.0:0.01:0.99;
R_alt = a*(1-e) - R;
idx = find(R_alt >= 400);

% Determine maximum eccentricity due to altitude constraint
e_max = -(600 + R)/a + 1;

e = e(idx);

i = acos(sqrt(1/5));
omega_dot = -3/2*n*J_2*(R/a)^2*cos(i)./((1-e.^2).^2);
plot(e,omega_dot)
title("Mars Molniya Orbit")
ylabel("Longitude of Ascending Node [rad]")
xlabel("e")