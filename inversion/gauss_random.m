clear;

%% Plot Lambda(r) curves for Intermediate Toy-Model

r=logspace(-4,0,200);
  
u = load(strcat('output/toy_model/gauss_random_f0.mat'));

lambda = u.asymmetry;

figure(1)
semilogx(r,lambda); hold on;
xlabel('Updraught Stability Reduction Factor r');
ylabel('Asymmetry Parameter {\lambda}');
title('Asymmetry Parameter for Gauss-Random f0')
set(gca)

u = load(strcat('output/toy_model/gauss_random_fy.mat'));

lambda = u.asymmetry;

figure(2)
semilogx(r,lambda); hold on;
xlabel('Updraught Stability Reduction Factor r');
ylabel('Asymmetry Parameter {\lambda}');
title('Asymmetry Parameter for Gauss-Random fy')
set(gca)


u = load(strcat('output/toy_model/gauss_random_Omega3.mat'));

lambda = u.asymmetry;

figure(3)
semilogx(r,lambda); hold on;
xlabel('Updraught Stability Reduction Factor r');
ylabel('Asymmetry Parameter {\lambda}');
title('Asymmetry Parameter for Gauss-Random Omega3')
set(gca)

u = load(strcat('output/toy_model/gauss_random_Ra3.mat'));

lambda = u.asymmetry;

figure(4)
semilogx(r,lambda); hold on;
xlabel('Updraught Stability Reduction Factor r');
ylabel('Asymmetry Parameter {\lambda}');
title('Asymmetry Parameter for Gauss-Random Ra3')
set(gca)

u = load(strcat('output/toy_model/sin_f0.mat'));

lambda = u.asymmetry;

figure(5)
semilogx(r,lambda); hold on;
xlabel('Updraught Stability Reduction Factor r');
ylabel('Asymmetry Parameter {\lambda}');
title('Asymmetry Parameter for Sine-Forcing f0')
set(gca)

u = load(strcat('output/toy_model/lambda_k.mat'));
v = load(strcat('output/toy_model/lambda_k_continued.mat'));

k = 1:500;

lambda = cat(1,u.asymmetry, v.asymmetry);

figure(6)
plot(k,lambda); hold on;
xlabel('Wavenumber k');
ylabel('Asymmetry Parameter {\lambda}');
title('Asymmetry Parameter for sin(kx)')
set(gca)




