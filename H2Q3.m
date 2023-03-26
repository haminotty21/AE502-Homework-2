clear all
close all

% Initialize constants and initial conditions
mu = 3.986*10^5;
omega_0 = 90/180*pi;        % radians
i_0 = 1.10654;              % radians
w_0 = 5/180*pi;             % radians
a_0 = 26600;                % km
e_0 =  0.74;    
M_0 = 10/180*pi;            % radians
tspan = [0 100*24*3600];    % seconds

% Compute E_0 via Algoritm 3.1
if M_0 <= pi
    E_0 = M_0 + e_0/2;
elseif M_0 > pi
    E = M_0 - e_0/2;
end

%   Solve Keplar's equation to calculate True Anomaly (theta)
E(1) = E_0;
thres = 10^-8;
FE = thres + 1;
inc = 1;
while(abs(FE) > thres)
    FE = (E(inc) - e_0*sin(E(inc))-M_0)/(1-e_0*cos(E(inc)));
    E(inc+1) = E(inc) - FE;
    inc = inc + 1;
end
E = E(end);
theta_0 = 2*atan2(sqrt(1+e_0)*tan(0.5*E),sqrt(1-e_0));

% Use perigee radius to determine angular momentum
r_p = a_0 * (1-e_0);
h_0 = sqrt((1+e_0)*mu*r_p);

% Form initial condition vector to input into ode45
x0 = [h_0, e_0, theta_0, omega_0, i_0, w_0];


options = odeset('AbsTol',1e-10,'RelTol',1e-13);
[t, x] = ode45(@Gauss_Planetary_Equations, tspan, x0,options);

% Reassign variables for readability
h = x(:,1);
e = x(:,2);
f = x(:,3);
LAN = x(:,4);
i = x(:,5);
w = x(:,6);

% Calculate semi-major axis and mean anomaly from given elements
a = h.^2./mu .* (1./(1-e.^2));

% Convert angle to degrees
LAN = LAN*180/pi;
i = i*180/pi;
w = w*180/pi;

figure
subplot(3,2,1)
plot(t/(3600*24),a )
ylabel('[km]')
xlabel('time (days)')
title('Semi-major Axis (a)')

subplot(3,2,2)
plot(t/(3600*24), e)
xlabel('time (days)')
title('Eccentricity (e)')

subplot(3,2,3)
plot(t/(3600*24),i)
xlabel('time (days)')
ylabel('[deg]')
title('Inclination (i)')

subplot(3,2,4)
plot(t/(3600*24),LAN)
xlabel('time (days)')
title('Longitude of ascending node (Omega)')
ylabel('[deg]')

subplot(3,2,5)
plot(t/(3600*24),w)
xlabel('time (days)')
ylabel('[deg]')
title('Argument of Periapsis (w)')