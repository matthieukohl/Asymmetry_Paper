% chat gpt test 3-body problem

% Three Body Problem Solver

% Define constants
G = 6.6743e-11; % Gravitational constant
m1 = 5.97e24; % Mass of Earth
m2 = 7.34e22; % Mass of Moon
m3 = 0; % Mass of spacecraft

% Define initial conditions
x1 = 0; % Earth's initial position
y1 = 0;
vx1 = 100; % Earth's initial velocity
vy1 = 0;
x2 = 384400e3; % Moon's initial position
y2 = 0;
vx2 = 0; % Moon's initial velocity
vy2 = 1022;
x3 = 0; % Spacecraft's initial position
y3 = 5000e3;
vx3 = 1000; % Spacecraft's initial velocity
vy3 = 0;

% Define time step and duration of simulation
dt = 60; % Time step (seconds)
tmax = 24*60*60*30; % Duration of simulation (seconds)

% Set up arrays to store position and velocity
nsteps = tmax/dt;
x1s = zeros(nsteps,1);
y1s = zeros(nsteps,1);
vx1s = zeros(nsteps,1);
vy1s = zeros(nsteps,1);
x2s = zeros(nsteps,1);
y2s = zeros(nsteps,1);
vx2s = zeros(nsteps,1);
vy2s = zeros(nsteps,1);
x3s = zeros(nsteps,1);
y3s = zeros(nsteps,1);
vx3s = zeros(nsteps,1);
vy3s = zeros(nsteps,1);

% Initialize position and velocity arrays with initial conditions
x1s(1) = x1;
y1s(1) = y1;
vx1s(1) = vx1;
vy1s(1) = vy1;
x2s(1) = x2;
y2s(1) = y2;
vx2s(1) = vx2;
vy2s(1) = vy2;
x3s(1) = x3;
y3s(1) = y3;
vx3s(1) = vx3;
vy3s(1) = vy3;

% Use fourth-order Runge-Kutta method to update position and velocity at each time step
for i = 2:nsteps
[ax1, ay1, ax2, ay2, ax3, ay3] = acceleration(G,m1,m2,m3,x1s(i-1), y1s(i-1), x2s(i-1), y2s(i-1), x3s(i-1), y3s(i-1));
kx1 = vx1s(i-1)*dt;
ky1 = vy1s(i-1)*dt;
kvx1 = ax1*dt;
kvy1 = ay1*dt;
kx2 = vx2s(i-1)*dt;
ky2 = vy2s(i-1)*dt;
kvx2 = ax2*dt;
kvy2 = ay2*dt;
kx3 = vx3s(i-1)*dt;
ky3 = vy3s(i-1)*dt;
kvx3 = ax3*dt;
kvy3 = ay3*dt;
x1s(i) = x1s(i-1) + (kx1 + 0.5*kvx1)*dt;
y1s(i) = y1s(i-1) + (ky1 + 0.5*kvy1)*dt;
vx1s(i) = vx1s(i-1) + kvx1;
vy1s(i) = vy1s(i-1) + kvy1;
x2s(i) = x2s(i-1) + (kx2 + 0.5*kvx2)*dt;
y2s(i) = y2s(i-1) + (ky2 + 0.5*kvy2)*dt;
vx2s(i) = vx2s(i-1) + kvx2;
vy2s(i) = vy2s(i-1) + kvy2;
x3s(i) = x3s(i-1) + (kx3 + 0.5*kvx3)*dt;
y3s(i) = y3s(i-1) + (ky3 + 0.5*kvy3)*dt;
vx3s(i) = vx3s(i-1) + kvx3;
vy3s(i) = vy3s(i-1) + kvy3;
end

% Plot the motion of the three bodies
figure;
plot(x1s,y1s,'b',x2s,y2s,'r',x3s,y3s,'g');
title('Motion of Three Bodies');
xlabel('x (m)');
ylabel('y (m)');
legend('Earth','Moon','Spacecraft');



% Define function to calculate acceleration
function [ax1, ay1, ax2, ay2, ax3, ay3] = acceleration(G,m1,m2,m3,x1, y1, x2, y2, x3, y3)
    r12 = sqrt((x1-x2)^2 + (y1-y2)^2); % Distance between Earth and Moon
    r13 = sqrt((x1-x3)^2 + (y1-y3)^2); % Distance between Earth and spacecraft
    r23 = sqrt((x2-x3)^2 + (y2-y3)^2); % Distance between Moon and spacecraft
    ax1 = -G*m2*(x1-x2)/r12^3 - G*m3*(x1-x3)/r13^3; % x-component of Earth's acceleration
    ay1 = -G*m2*(y1-y2)/r12^3 - G*m3*(y1-y3)/r13^3; % y-component of Earth's acceleration
    ax2 = G*m1*(x1-x2)/r12^3 - G*m3*(x2-x3)/r23^3; % x-component of Moon's acceleration
    ay2 = G*m1*(y1-y2)/r12^3 - G*m3*(y2-y3)/r23^3; % y-component of Moon's acceleration
    ax3 = G*m1*(x1-x3)/r13^3 + G*m2*(x2-x3)/r23^3; % x-component of spacecraft's acceleration
    ay3 = G*m1*(y1-y3)/r13^3 + G*m2*(y2-y3)/r23^3; % y-component of spacecraft's acceleration
end


