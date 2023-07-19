close all; clear; 

function three_body_problem()
% Define the initial conditions
x0 = [-0.970, 0.243, 0];
y0 = [0.466, 0.927, -0.523];
vx0 = [0.466, 0.927, -0.523];
vy0 = [-0.970, 0.243, 0];

% Define the mass ratios
m1 = 1;
m2 = 1;
m3 = 1;

% Combine the initial conditions into a single vector
y0 = [x0, y0, vx0, vy0];

% Set up the ODE solver options
tspan = [0, 10];
options = odeset('RelTol', 1e-8, 'AbsTol', 1e-8);

% Define the ODE system function
function dydt = odefun(t, y)
    % Unpack the state variables
    x1 = y(1);
    x2 = y(2);
    x3 = y(3);
    y1 = y(4);
    y2 = y(5);
    y3 = y(6);
    vx1 = y(7);
    vx2 = y(8);
    vx3 = y(9);
    vy1 = y(10);
    vy2 = y(11);
    vy3 = y(12);
    
    % Calculate the distances between the bodies
    r12 = sqrt((x1 - x2)^2 + (y1 - y2)^2);
    r13 = sqrt((x1 - x3)^2 + (y1 - y3)^2);
    r23 = sqrt((x2 - x3)^2 + (y2 - y3)^2);
    
    % Calculate the accelerations of the bodies
    ax1 = -m2 * (x1 - x2) / r12^3 - m3 * (x1 - x3) / r13^3;
    ay1 = -m2 * (y1 - y2) / r12^3 - m3 * (y1 - y3) / r13^3;
    ax2 = -m1 * (x2 - x1) / r12^3 - m3 * (x2 - x3) / r23^3;
    ay2 = -m1 * (y2 - y1) / r12^3 - m3 * (y2 - y3) / r23^3;
    ax3 = -m1 * (x3 - x1) / r13^3 - m2 * (x3 - x2) / r23^3;
    ay3 = -m1 * (y3 - y1) / r13^3 - m2 * (y3 - y2) / r23^3;
    
    % Combine the derivatives into a single vector
    dydt = [vx1; vx2; vx3; vy1; vy2; vy3; ax1; ax2; ax3; ay1; ay2; ay3];
end

% Solve the ODE system using ode45
[t, y] = ode45(@odefun, tspan, y0, options);

% Plot the trajectories of the three bodies
figure;
plot(y(:, 1), y(:, 4), '-b',y(:, 2), y(:, 5), '-r', y(:, 3), y(:, 6), '-g');
title('Trajectories of the three bodies');
xlabel('x');
ylabel('y');
legend('Body 1', 'Body 2', 'Body 3');

end