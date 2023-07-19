clear, clc, close all

% This example shows you how to use two black 
% noise signals to generate a landscape surface

% generate two vectors of black noise
N = 1000;
z1 = arbssnoise(N, -5);
z2 = arbssnoise(N, -5);

% terraforming
Z = z1*z2';

% form the grid
[X, Y] = meshgrid(1:N, 1:N);

% plot the surface
figure(1)
surfc(X, Y, Z)
shading interp
colormap jet
axis tight
set(gca, 'FontName', 'Times New Roman', 'FontSize', 12)
xlabel('Longitude, km')
ylabel('Latitude, km')
zlabel('Altitude, km')
alpha(0.8)