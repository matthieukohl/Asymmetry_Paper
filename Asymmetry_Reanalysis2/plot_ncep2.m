% plot lambda from ncep2 using paul's script and my script

close all; clear

% Northern-Hemisphere

load('data_ncep2/asymmetry_paul_xt.mat');

[~,lat30] = min(abs(lat-30));
[~,lat70] = min(abs(lat-70));

[~,p500] = min(abs(level-500));

lambda_paul_ncep2 = squeeze(mean(lambda_paul(lat70:lat30,p500,:),1));

clear('level','lat','p500','lat30','lat70','lambda_paul');

load('data_ncep2/asymmetry_xt_erai_nh_paul.mat');

[~,lat30] = min(abs(lat-30));
[~,lat70] = min(abs(lat-70));

[~,p500] = min(abs(level-500));

lambda_paul_erai = squeeze(mean(lambda_paul(lat70:lat30,p500,:),1));

clear('level','lat','p500','lat30','lat70','lambda_paul');

load('data_ncep2/asymmetry_xt_erai_low_res_paul.mat');

[~,lat30] = min(abs(lat-30));
[~,lat70] = min(abs(lat-70));

[~,p500] = min(abs(level-500));

lambda_paul_erai_low_res = squeeze(mean(lambda_paul(lat70:lat30,p500,:),1));

clear('level','lat','p500','lat30','lat70','lambda_paul');

load('data_ncep2/asymmetry_xt_erai_low_res_92_01_paul.mat');

[~,lat30] = min(abs(lat-30));
[~,lat70] = min(abs(lat-70));

[~,p500] = min(abs(level-500));

lambda_paul_erai_low_res_92_01 = squeeze(mean(lambda_paul(lat70:lat30,p500,:),1));

clear('level','lat','p500','lat30','lat70','lambda_paul');

load('data_ncep2/asymmetry_xt_erai_92_95_nh_paul.mat');

[~,lat30] = min(abs(lat-30));
[~,lat70] = min(abs(lat-70));

[~,p500] = min(abs(level-500));

lambda_paul_erai_92_95 = squeeze(mean(lambda_paul(lat70:lat30,p500,:),1));

clear('level','lat','p500','lat30','lat70','lambda_paul');


figure(1)
plot(lambda_paul_ncep2,'Linewidth',1.5); hold on
plot(lambda_paul_erai,'Linewidth',1.5,'Color','r'); hold on;
plot(lambda_paul_erai_low_res,'linewidth',1.5,'Color','r','linestyle',':'); hold on;
plot(lambda_paul_erai_92_95,'linewidth',1.5,'Color','g'); hold on;
plot(lambda_paul_erai_low_res_92_01,'linewidth',1.5,'Color','g','Linestyle',':'); hold on;
ylabel('\lambda')
set(gca,'xtick',1:12,...
 'xticklabel',{'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'})
set(gca,'Fontsize',10)
legend('NCEP2(92-01)','ERAI(10-14)','ERAI (10-14)-Low Res','ERAI(92-95)','ERAI(92-01)-Low Res','Location','Northwest'); legend boxoff
title('Lambda 30-70 NH')

% Southern-Hemisphere

load('data_ncep2/asymmetry_paul_xt.mat');

[~,lat30] = min(abs(lat+30));
[~,lat70] = min(abs(lat+70));

[~,p500] = min(abs(level-500));

lambda_paul_ncep2 = squeeze(mean(lambda_paul(lat30:lat70,p500,:),1));

clear('level','lat','p500','lat30','lat70','lambda_paul');

load('data_ncep2/asymmetry_xt_erai_sh_paul.mat');

[~,lat30] = min(abs(lat+30));
[~,lat70] = min(abs(lat+70));

[~,p500] = min(abs(level-500));

lambda_paul_erai = squeeze(mean(lambda_paul(lat30:lat70,p500,:),1));

clear('level','lat','p500','lat30','lat70','lambda_paul');

load('data_ncep2/asymmetry_xt_erai_low_res_paul.mat');

[~,lat30] = min(abs(lat+30));
[~,lat70] = min(abs(lat+70));

[~,p500] = min(abs(level-500));

lambda_paul_erai_low_res = squeeze(mean(lambda_paul(lat30:lat70,p500,:),1));

clear('level','lat','p500','lat30','lat70','lambda_paul');

load('data_ncep2/asymmetry_xt_erai_low_res_92_01_paul.mat');

[~,lat30] = min(abs(lat+30));
[~,lat70] = min(abs(lat+70));

[~,p500] = min(abs(level-500));

lambda_paul_erai_low_res_92_01 = squeeze(mean(lambda_paul(lat30:lat70,p500,:),1));

clear('level','lat','p500','lat30','lat70','lambda_paul');

load('data_ncep2/asymmetry_xt_erai_92_95_sh_paul.mat');

[~,lat30] = min(abs(lat+30));
[~,lat70] = min(abs(lat+70));

[~,p500] = min(abs(level-500));

lambda_paul_erai_92_95 = squeeze(mean(lambda_paul(lat30:lat70,p500,:),1));

clear('level','lat','p500','lat30','lat70','lambda_paul');

figure(2)
plot(lambda_paul_ncep2,'Linewidth',1.5); hold on
plot(lambda_paul_erai,'Linewidth',1.5,'Color','r'); hold on;
plot(lambda_paul_erai_low_res,'linewidth',1.5,'Color','r','linestyle',':'); hold on;
plot(lambda_paul_erai_92_95,'linewidth',1.5,'Color','g'); hold on;
plot(lambda_paul_erai_low_res_92_01,'linewidth',1.5,'Color','g','Linestyle',':'); hold on;
ylabel('\lambda')
set(gca,'xtick',1:12,...
 'xticklabel',{'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'})
set(gca,'Fontsize',10)
legend('NCEP2(92-01)','ERAI(10-14)','ERAI (10-14)-Low Res','ERAI(92-95)','ERAI(92-01)-Low Res','Location','Northeast'); legend boxoff
title('Lambda 30-70 SH')


% figure(2)
% plot(lambda_paul_ncep2,'Linewidth',1.5); hold on
% plot(lambda_paul_erai,'Linewidth',1.5); hold on;
% plot(lambda_paul_erai_low_res,'linewidth',1.5);
% ylabel('\lambda')
% set(gca,'xtick',1:12,...
%  'xticklabel',{'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'})
% set(gca,'Fontsize',12)
% legend('NCEP2(92-01)','ERAI(10-14)','ERAI (10-14)-Low Res','Location','North'); legend boxoff
% title('Lambda 40-60 SH')

% Plot lambda vs lat


load('data_ncep2/asymmetry_paul_xt.mat');

[~,p500] = min(abs(level-500));

lambda_ncep2_djf = squeeze(mean(lambda_paul(:,p500,[12,1,2]),3));
lambda_ncep2_jja = squeeze(mean(lambda_paul(:,p500,[6,7,8]),3));
clear('lat','level','lambda_paul')

%load('data_ncep2/asymmetry_xt_erai_low_res_paul.mat');

load('data_ncep2/asymmetry_xt_erai_low_res_92_01_paul.mat');

[~,p500] = min(abs(level-500));

lambda_erai_low_res_djf = squeeze(mean(lambda_paul(:,p500,[12,1,2]),3));
lambda_erai_low_res_jja = squeeze(mean(lambda_paul(:,p500,[6,7,8]),3));

figure(3)
plot(lat,lambda_ncep2_djf,'Color','b','linewidth',1.5); hold on;
plot(lat,lambda_erai_low_res_djf,'Color','b','linewidth',1.5,'linestyle',':'); hold on;
plot(lat,lambda_ncep2_jja,'linewidth',1.5,'Color','r'); hold on;
plot(lat,lambda_erai_low_res_jja,'linewidth',1.5,'Color','r','linestyle',':'); hold on;
xlabel('lat'); xlim([-90 90]); 
legend('djf-ncep2','djf-erai','jja-ncep2','jja-erai')
legend boxoff
title('500hPa - 92/01(low-resolution)')



