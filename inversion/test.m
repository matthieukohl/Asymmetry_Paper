
% Constants
G = 6.67408e-11; % Gravitational constant
M_sun = 1.989e30; % Mass of the sun
M_earth = 5.97e24; % Mass of the Earth
R_earth = 6.371e6; % Radius of the Earth
AU = 1.496e11; % Astronomical unit

% Initial conditions
x_earth = AU; % Earth starts at 1 AU from the Sun
y_earth = 0;
vx_earth = 0; % Earth starts with zero velocity
vy_earth = 29.78e3; % Earth's initial velocity is its orbital velocity

% Time parameters
t_start = 0;
t_end = 100*365*24*3600; % 100 years in seconds
dt = 3600; % Time step of 1 hour

% Simulation loop
t = t_start;
while t < t_end
    % Calculate distance between Earth and Sun
    r_earth_sun = sqrt(x_earth^2 + y_earth^2);
    
    % Calculate gravitational force on Earth
    F_gravity = G * M_sun * M_earth / r_earth_sun^2;
    
    % Calculate gravitational acceleration on Earth
    ax_earth = -F_gravity * x_earth / r_earth_sun / M_earth;
    ay_earth = -F_gravity * y_earth / r_earth_sun / M_earth;
    
    % Update Earth's position and velocity
    x_earth = x_earth + vx_earth * dt;
    y_earth = y_earth + vy_earth * dt;
    vx_earth = vx_earth + ax_earth * dt;
    vy_earth = vy_earth + ay_earth * dt;
    
    % Check for asteroid impact
    if r_earth_sun < R_earth % Earth collided with the Sun
        disp("The Earth collided with the Sun!");
        break;
    elseif r_earth_sun < 1e8 % An asteroid collided with the Earth
        disp("An asteroid collided with the Earth and ejected it from the solar system!");
        break;
    end
    
    % Update time
    t = t + dt;
end

% Plot Earth's orbit
theta = linspace(0, 2*pi);
x_sun = zeros(size(theta));
y_sun = zeros(size(theta));
x_earth_orbit = AU * cos(theta);
y_earth_orbit = AU * sin(theta);
plot(x_sun, y_sun, 'yo', x_earth_orbit, y_earth_orbit, 'b-', x_earth, y_earth, 'ro');
axis equal;
xlabel('x (m)');
ylabel('y (m)');
title('Earth orbit simulation');
legend('Sun', 'Earth orbit', 'Earth');