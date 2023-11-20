% brownian motion settling with gravity

% Parameters
numParticles = 100;         % Number of particles to simulate
timeSteps = 1000;           % Number of time steps to simulate
diffusionCoefficient = 1;   % Diffusion coefficient
particleRadius = 1;         % Particle radius
gravity = 9.8;              % Acceleration due to gravity

% Initialize particle positions
particlePositions = rand(numParticles, 2);

% Initialize particle velocities
particleVelocities = zeros(numParticles, 2);

% Simulation loop
for i = 1:timeSteps
    % Calculate the Brownian motion-based diffusion component
    diffusionComponent = sqrt(2 * diffusionCoefficient) * randn(numParticles, 2);
    
    % Calculate the gravitational component
    gravityComponent = repmat([0, -gravity], numParticles, 1);
    
    % Update particle positions and velocities
    particlePositions = particlePositions + particleVelocities + diffusionComponent + gravityComponent;
    particleVelocities = particleVelocities + diffusionComponent + gravityComponent;
    
    % Enforce boundary conditions (particles cannot leave the simulation box)
    particlePositions(particlePositions < particleRadius) = particleRadius;
    particlePositions(particlePositions > 1 - particleRadius) = 1 - particleRadius;
    
    % Plot particle positions
    scatter(particlePositions(:, 1), particlePositions(:, 2));
    xlim([0, 1]);
    ylim([0, 1]);
    drawnow;
end
