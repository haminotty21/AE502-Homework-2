function x = Gauss_Planetary_Equations(t,x)
%% Gauss_Planetary_Equations
% Author: Darrell Hamilton
% Date: 3/10/2023
% Description: Take in inputs and forms the Planetary Equations for the
%   elements listed in the input vector x. This is to be used with a
%   numerical integrator such as ode45

% Inputs: x(1) : Angular Momentum             (h)
%         x(2) : eccentricity                 (e)
%         x(3) : true anomaly                 (theta)
%         x(4) : Longitude of Ascending Node  (omega)
%         x(5) : Inclination                  (i)
%         x(6) : Argument of periapsis        (w)

% Outputs:  Evaluated planetary equations for
%           x(1) : Angular Momentum
%           x(2) : eccentricity                 
%           x(3) : true anomaly                 
%           x(4) : Longitude of Ascending Node  
%           x(5) : Inclination                
%           x(6) : Argument of periapsis       

%%
mu = 3.986*10^5;    % Gravitational constant for earth
J_2 = 0.00108;      % J2 pertubation constant for earth
R = 6370;           % Earth's radius (km)

% Assign variables to input values to increase readability
h = x(1);
e = x(2);
theta = x(3);
omega = x(4);
i = x(5);
w = x(6);

% Calculate radius magnitude to satellite
r = h^2/mu * 1/(1 + e*cos(theta));


U = w + theta;
% Use Gauss's Equations to determine perturbations
dh = -3/2 * J_2*mu*R^2 /r^3 * sin(i)^2 * sin(2*U);           % Check last term

de = 3/2*J_2*mu*R^2/h/r^3 ...
*(h^2/mu/r*sin(theta)*(3*sin(i)^2*sin(U)^2 - 1) ...
-sin(2*U)*sin(i)^2*((2+e*cos(theta))*cos(theta)+e));

dtheta = h/r^2 + 3/2*J_2*mu*R^2 / (e*h*r^3) * (h^2/(mu*r) * cos(theta)*(3*sin(i)^2*sin(U)^2 - 1) ...
    + (2 + e* cos(theta))*sin(2*U)*sin(i)^2*sin(theta));
domega = -3 * J_2*mu*R^2 / (h*r^3) * sin(U)^2*cos(i);
di = -3/4 * J_2*mu*R^2 / (h*r^3) * sin(2*U)*sin(2*i);
dw = 3/2 * J_2*mu*R^2 / (e*h*r^3) * (h^2/(mu*r) * cos(theta)*(1 - 3*sin(i)^2*sin(U)^2) ... 
    - (2 + e*cos(theta))*sin(2*U)*sin(i)^2*sin(theta) + 2*e*cos(i)^2*sin(U)^2);


% Reassign x for output
x = [dh; de; dtheta; domega; di; dw];

end


