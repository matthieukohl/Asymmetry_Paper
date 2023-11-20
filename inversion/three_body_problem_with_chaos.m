% chaos 3 body problem 

% Three-Body Problem with Chaos
% Using Fourth-Order Runge-Kutta Method

% Define masses of the bodies
m1 = 1;
m2 = 1;
m3 = 1;

% Define initial positions and velocities
r1_0 = [-0.97000436; 0.24308753];
r2_0 = [0; 0] + randn;
r3_0 = [0.97000436; -0.24308753];
v1_0 = [0.4662036850; 0.4323657300];
v2_0 = [-0.9324073700; -0.8647314600];
v3_0 = [0.4662036850; 0.4323657300];

% Combine initial positions and velocities into state vector
y0 = [r1_0; r2_0; r3_0; v1_0; v2_0; v3_0];

% Define time step and total simulation time
dt = 0.01;
tspan = 0:dt:100;

% Define function for calculating derivatives
f = @(t, y) three_body_derivatives(t, y, m1, m2, m3);

% Use fourth-order Runge-Kutta method to simulate system
[t, y] = ode45(f, tspan, y0);

% Plot the trajectories of the three bodies
figure
hold on
plot(y(:,1), y(:,2), 'r')
plot(y(:,3), y(:,4), 'b')
plot(y(:,5), y(:,6), 'g')
xlabel('x')
ylabel('y')
title('Three-Body Problem with Chaos')
legend('Body 1', 'Body 2', 'Body 3')

% Define function for calculating derivatives
function dydt = three_body_derivatives(t, y, m1, m2, m3)
    % Unpack state vector
    r1 = y(1:2);
    r2 = y(3:4);
    r3 = y(5:6);
    v1 = y(7:8);
    v2 = y(9:10);
    v3 = y(11:12);

    % Calculate distances between bodies
    d12 = norm(r1 - r2);
    d13 = norm(r1 - r3);
    d23 = norm(r2 - r3);

    % Calculate derivatives
    dr1dt = v1;
    dr2dt = v2;
    dr3dt = v3;
    dv1dt = ((m2 * (r2 - r1) / d12^3) + (m3 * (r3 - r1) / d13^3));
    dv2dt = ((m1 * (r1 - r2) / d12^3) + (m3 * (r3 - r2) / d23^3));
    dv3dt = ((m1 * (r1 - r3) / d13^3) + (m2 * (r2 - r3) / d23^3));
    dydt = [dr1dt; dr2dt; dr3dt; dv1dt; dv2dt; dv3dt];
end