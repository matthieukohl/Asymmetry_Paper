%% Compare toy-model results to r-simulations and abs-simulations

% r-simulations wavenumbers

r = [0.01,0.02,0.05,0.1,0.2,0.3,0.4,0.6,0.8,1.0];

n_W_r = [13.5621,13.1434,12.7994,11.7083,9.9645,8.7613,8.0796,7.4070,6.9377,6.2743];

n_RHS_r = [29.8946,24.7720,23.5225,19.5714,15.2045,12.4606,11.3807,9.9027,9.0329,8.1522];

lambda_omega_r = [0.7034,0.7035,0.7073,0.7184,0.7012,0.6598,0.6274,0.5752,0.5305,0.4997];


lambda_omegaQG_r= [0.7439,0.7361,0.7202,0.7093,0.6874,0.6536,0.6258,0.5797,0.5342,0.5020];

%lambda_omegaQG_r = [0.7489,0.7333,0.7111,0.6939,0.6716,0.6391,0.6125,0.5740,0.5261,0.5020]; 

% abs-simulations wavenumbers

list = {'abs0.4_T85','abs0.6_T85','abs0.7_T85','abs0.8_T85','abs0.9_T85',...
    'abs1.0_T85','abs1.2_T85','abs1.4_T85','abs1.6_T85','abs1.8_T85',...
    'abs2.0_T85','abs2.5_T85','abs3.0_T85','abs4.0_T85','abs6.0_T85'};

surf_temp = [270.0311,277.5996,280.5327,283.1391,285.4284,287.6549,290.8506...
293.6622, 296.1045, 298.1490, 300.0797, 303.7376, 306.5808 ...
310.7297, 316.1383];

n_W_abs = [11.9953,12.5895,12.6759,13.1430,12.9818,13.1746,12.9951,13.4855,...
13.2733,13.1440,13.2256,12.7880,12.2098,10.7010,9.7016];

n_RHS_abs = [21.6518,22.1836,21.8106,22.3484,21.9664,22.4217,22.0013,22.4879,...
22.6347,22.6507,22.4226,22.3218,21.1196,19.5400,19.2665];

r_abs = [0.8113,0.7093,0.6707,0.6214,0.6000,0.5561,0.5299,0.4609,0.4425,0.3959,...
    0.3944,0.3204,0.3116,0.2757,0.2462];

r_abs_500 = [0.8296,0.7396,0.6776,0.6410,0.5950,0.5339,0.4676,0.3686,0.3044,...
    0.2386,0.2001,0.1244,0.0836,0.0399,0.0086];

lambda_omega_abs= [0.6176,0.6363,0.6371,0.6421,0.6598,0.6496,0.6517,0.6560,...
    0.6587,0.6611,0.6642,0.6761,0.6708,0.6560,0.6527];

lambda_omegaQG_abs = [0.5463,0.5679,0.5780,0.5846,0.6068,0.6027,0.6135,0.6206,...
    0.6273,0.6398,0.6429,0.6673,0.6765,0.6801,0.7013];

lambda_omega_abs_500 = [0.5845,0.6235,0.6429,0.6510,0.6825,0.6713,0.6837,0.6894,...
    0.6986,0.7012,0.7041,0.7137, 0.7107,0.6873,0.6699];

lambda_omegaQG_abs_500 = [0.5338,0.5549,0.5745,0.5810,0.6136,0.6131,0.6343,...
    0.6447,0.6615,0.6783,0.6847,0.7132,0.7352,0.7337,0.7560];

lambda_omega_abs_500_fulldays = [0.5959,0.6234,0.6383,0.6519,0.6676,0.6775,...
    0.6829,0.6918,0.6987,0.6996,0.7021,0.7030,0.7078,0.6929,0.6678];

lambda_omegaQG_abs_500_fulldays = [0.5408,0.5571,0.5695,0.5823,0.5983,0.6150,0.6313,...
    0.6464,0.6634,0.6709,0.6837,0.7053,0.7307,0.7380,0.7599];

% Pre-iniatilize 

lambda_r = zeros(length(r),1);
lambda_abs = zeros(length(list),1);

% Calculate lambda_theory based on the toy-model

% r-simulation

for ii=1:length(r)
    lambda_r(ii) = toy_model1(r(ii),n_W_r(ii));
end

lambda_r_rand = [0.7535, 0.7535,0.7518,0.7184,0.6639,0.6259,...  
0.6038, 0.5383, 0.5217,0.5003];

% abs-simulation

for jj=1:length(list)
   lambda_abs(jj) = toy_model1(r_abs_500(jj),n_W_abs(jj));
end

% lambda_abs_rand = [0.4710,0.5560,0.5690,0.5626,0.4330,0.6503,0.5888,...
% 0.6297,0.6351,0.6267,0.6415,0.7451,0.6629,0.7939, 0.8061];

lambda_abs_rand = [0.5425,0.5781,0.5917,0.5857,0.5826,0.6167,0.6294,0.6607,...
    0.6583,0.6867,0.6988,0.7271,0.7771,0.7961,0.8386];

lambda_abs_rand_500 = [0.5714,0.5751,0.6055,0.6108,0.6182,0.6425,0.6701,...
0.6942,0.7129,0.7197,0.7542,0.7945,0.8256,0.8471,0.8760];

%lambda_abs_rand = [0.5550,0.5795,0.5796,0.6178,0.6137,0.6344,0.6271,...
%0.6533,0.6675,0.6894,0.6931,0.7413,0.7668,0.7911,0.8273];

% Plot Asymmetries r-simulation

figure(1)
semilogx(r,lambda_omega_r,'Color','blue','linewidth',1.6); hold on;
semilogx(r,lambda_omegaQG_r,'Color','blue','LineStyle','--','linewidth',1.6); hold on
semilogx(r,lambda_r,'Color','red','LineStyle','--','linewidth',1.6); hold on;
semilogx(r,lambda_r_rand,'red','LineStyle','-.','linewidth',1.6); hold on;
xlabel('Reduction factor r','interpreter','latex')
ylabel('Asymmetry parameter $\lambda$','interpreter','latex')
legend('GCM','QG','1-D Toy Model','White Noise')
title('Toy-model prediction')
set(gca,'fontsize', 14);
set(gca, 'TickDir', 'out','Box', 'off','Layer', 'top')
set(gca,'linewidth',1.5)
ylim([0.5,1])
legend boxoff
%saveas(gcf,'Plots_AOFD/Toy_Model_Predictions1','epsc')

figure(2)
semilogx(r,1000/(cosd(46)*6371)*n_W_r,'linewidth',1.6)
xlabel('Reduction factor r','interpreter','latex')
ylabel('Nondimensional k','interpreter','latex')
title('Dominant Wavenumber ')
set(gca,'fontsize', 14);
set(gca, 'TickDir', 'out','Box', 'off','Layer', 'top')
set(gca,'linewidth',1.5)
ylim([1.4 3.5])
%saveas(gcf,'Plots_AOFD/Dominant_Wavenumber1','epsc')

% Plot Asymmetries abs-simulation

figure(3)
plot(surf_temp,lambda_omega_abs_500_fulldays,'Color','blue'); hold on;
plot(surf_temp,lambda_omegaQG_abs_500_fulldays,'Color','blue','LineStyle','--'); hold on
plot(surf_temp,lambda_abs,'Color','red','LineStyle','--'); hold on;
%plot(surf_temp,lambda_abs_rand,'red','LineStyle','-.'); hold on;
xlabel('T_{surf}(K)')
ylabel('Asymmetry parameter \lambda')
legend('omega','QG','Toy Model','Location','Northwest')
title('abs-simulation')
%saveas(gcf,'Plots_Generals/final/Toy_Model_sin_abs','epsc')


% figure(4)
% plot(surf_temp,n_W_abs)
% xlabel('r')
% ylabel('Wavenumber')
% title('Wavenumber abs simulation')