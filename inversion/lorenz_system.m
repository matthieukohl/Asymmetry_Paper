% Define the Lorenz equations
sigma = 10;
beta = 8/3;
rho = 28;
f = @(t, Y) [sigma*(Y(2)-Y(1)); Y(1)*(rho-Y(3))-Y(2); Y(1)*Y(2)-beta*Y(3)];

% Define the Jacobian matrix analytically
J = @(Y) [-sigma, sigma, 0; rho-Y(3), -1, -Y(1); Y(2), Y(1), -beta];

% Set initial conditions and time range
Y0 = [1, 1, 1];
tspan = [0 50];

% Solve the Lorenz equations using ode45
[t,Y] = ode45(f, tspan, Y0);

% Compute the tangent matrix at each point in the solution
n = size(Y,1);
A = zeros(3,3,n);
for i = 1:n
    A(:,:,i) = J(Y(i,:));
end

% Compute the largest singular values of the tangent matrix at each point
S = zeros(n,3);
for i = 1:n
    S(i,:) = svd(A(:,:,i));
end

% Compute the pointwise growth rate of uncertainty
growth_rate = sqrt(sum(S.^2,2));

% Sort the points by growth rate and plot the 10 largest
[~, idx] = sort(growth_rate, 'descend');
plot3(Y(idx(1:10),1), Y(idx(1:10),2), Y(idx(1:10),3), 'r.', 'MarkerSize', 20);
hold on;
plot3(Y(:,1), Y(:,2), Y(:,3), 'k');
xlabel('x');
ylabel('y');
zlabel('z');
title('Lorenz system with the 10 largest points of uncertainty growth');


% Define the Lorenz system parameters
sigma = 10;
beta = 8/3;
rho = 28;

% Define the initial conditions
x0 = [1; 1; 1];

% Define the time span for the simulation
tspan = [0, 100];

% Define the perturbation vector
eps = 1e-6;
dx0 = eps*[1; 0; 0];

% Solve the system of differential equations using ode45
[t, x] = ode45(@(t,x) lorenz(t,x,sigma,beta,rho), tspan, x0);
[~, x_eps] = ode45(@(t,x) lorenz(t,x,sigma,beta,rho), tspan, x0+dx0);

% Compute the Jacobian matrix
J = zeros(3,3,length(t));
for i = 1:length(t)
    J(:,:,i) = jacobian(t(i),x(i,:)',sigma,beta,rho);
end

% Compute the norm of the perturbation vector over time
norm_dx = zeros(length(t),1);
for i = 1:length(t)
    norm_dx(i) = norm(dx0);
end

% Compute the maximum Lyapunov exponent
lambda_max = max(sum(log(sqrt(sum((J.*J),2))),1)./norm_dx);

% Find the point where the uncertainty growth is largest
[~, idx] = max(sum(log(sqrt(sum((J.*J),2))),1)./norm_dx);
x_max = x(idx,:);

% Plot the trajectory of the Lorenz attractor and the point of maximum uncertainty growth
plot3(x(:,1), x(:,2), x(:,3), 'LineWidth', 1.5);
hold on;
plot3(x_max(1), x_max(2), x_max(3), 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r');
xlabel('x');
ylabel('y');
zlabel('z');
title(sprintf('Lorenz Attractor: max Lyapunov exponent = %.4f', lambda_max));
grid on;

function dxdt = lorenz(t, x, sigma, beta, rho)
% Function for the Lorenz system of differential equations
dxdt = zeros(3,1);
dxdt(1) = sigma * (x(2) - x(1));
dxdt(2) = x(1) * (rho - x(3)) - x(2);
dxdt(3) = x(1) * x(2) - beta * x(3);
end

function J = jacobian(t, x, sigma, beta, rho)
% Function for the Jacobian matrix of the Lorenz system
J = [-sigma, sigma, 0; rho-x(3), -1, -x(1); x(2), x(1), -beta];
end

% % Generate random numbers
% data = randn(100);
% 
% % Define smoothing kernel
% sigma = 2;
% n = 6*sigma + 1;
% kernel = fspecial('gaussian', [n n], sigma);
% 
% % Apply smoothing filter
% data_smooth = conv2(data, kernel, 'same');
% 
% % Plot the original and smoothed data
% figure;
% subplot(1,2,1);
% imagesc(data);
% title('Original Data');
% axis equal tight;
% colormap jet;
% colorbar;
% 
% subplot(1,2,2);
% imagesc(data_smooth);
% title('Smoothed Data');
% axis equal tight;
% colormap jet;
% colorbar;
% 
% % Jeopardy game
% 
% % Define categories and questions
% categories = {'Famous Scientists', 'US Presidents', 'Geography', 'Sports', 'Movies'};
% questions = {...
%     'This scientist is known for his theory of relativity.', ...
%     'This president served during the Civil War and was assassinated.', ...
%     'This country is home to the Great Wall.', ...
%     'This athlete won six Olympic gold medals in sprinting.', ...
%     'This movie won the Best Picture Oscar in 2020.'};
% answers = {...
%     'Who is Albert Einstein?', ...
%     'Who is Abraham Lincoln?', ...
%     'What is China?', ...
%     'Who is Usain Bolt?', ...
%     'What is "Parasite"?'};
% 
% % Randomly select a category and question
% num_categories = length(categories);
% num_questions = length(questions);
% category_index = randi(num_categories);
% question_index = randi(num_questions);
% category = categories{category_index};
% question = questions{question_index};
% answer = answers{question_index};
% 
% % Display the clue
% clue = [category ': ' question];
% fprintf('Clue: %s\n', clue);
% 
% % Get user input
% response = input('Your answer: ', 's');
% 
% % Check the answer
% if strcmpi(response, answer)
%     fprintf('Correct!\n');
% else
%     fprintf('Sorry, the correct answer is: %s\n', answer);
% end
% 
% 
