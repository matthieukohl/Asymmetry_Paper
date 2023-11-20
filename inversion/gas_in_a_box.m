% Gas in a box simulation using Monte Carlo

% Set parameters
N = 1000; % Number of particles
V = 1; % Volume of the box
L = V^(1/3); % Length of the box
T = 1; % Temperature
dt = 0.01; % Time step
tmax = 10; % Maximum simulation time
kB = 1; % Boltzmann constant
h = 1; % Planck's constant
m = 1; % Mass of a particle
lambda = sqrt(h^2/(2*pi*m*kB*T)); % Thermal de Broglie wavelength

% Initialize particle positions and velocities
x = L*rand(N, 3);
v = randn(N, 3)*sqrt(T);

% Initialize arrays for pressure, temperature, density, and entropy
P = zeros(1, tmax/dt);
T_array = zeros(1, tmax/dt);
rho = zeros(1, tmax/dt);
S = zeros(1, tmax/dt);

% Start simulation
for i = 1:tmax/dt
    
    % Move particles using Monte Carlo method
    for j = 1:N
        x_new = x(j,:) + randn(1,3)*sqrt(dt);
        if all(x_new > 0) && all(x_new < L)
            x(j,:) = x_new;
        end
    end
    
    % Calculate pressure, temperature, density, and entropy
    V = L^3;
    V_free = V; % - sum(x < 0, 1) - sum(x > L, 1); %there is a problem in this line
    P(i) = N*kB*T/V_free;
    T_array(i) = mean(sum(v.^2, 2))/(3*N*kB);
    rho(i) = N/V;
    S(i) = -N*kB*(rho(i)*log(rho(i)*lambda^3) - 5/2);
    
end

% Create plot of pressure, temperature, density, and entropy versus time
figure
subplot(2,2,1)
plot(0:dt:tmax-dt, P)
xlabel('Time')
ylabel('Pressure')
subplot(2,2,2)
plot(0:dt:tmax-dt, T_array)
xlabel('Time')
ylabel('Temperature')
subplot(2,2,3)
plot(0:dt:tmax-dt, rho)
xlabel('Time')
ylabel('Density')
subplot(2,2,4)
plot(0:dt:tmax-dt, S)
xlabel('Time')
ylabel('Entropy')

% Add labels outside of the figure
annotation('textbox', [0.1, 0.93, 0.05, 0.05], 'String', '(a)', 'EdgeColor', 'none')
annotation('textbox', [0.49, 0.93, 0.05, 0.05], 'String', '(b)', 'EdgeColor', 'none')
annotation('textbox', [0.1, 0.43, 0.05, 0.05], 'String', '(c)', 'EdgeColor', 'none')
annotation('textbox', [0.49, 0.43, 0.05, 0.05], 'String', '(d)', 'EdgeColor', 'none')


