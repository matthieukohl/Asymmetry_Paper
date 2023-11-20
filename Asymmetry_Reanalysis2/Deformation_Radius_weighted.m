function [L_mid,L_500,L_av] = Deformation_Radius_weighted(sigma,f,level,lat,top,weight)

% Calculate the Deformation Radius to nondimensionalize Wavenumber

sigma = squeeze(sum(sum(weight.*sigma,1),2)./(sum(sum(weight,1),2)));

% Midpoint Method

% Static stability

surf = 1;
[~,I] = min(abs(level-(level(surf)-(level(surf)-level(top))/2)));

S_mid = sigma(I);

[~,p500] = min(abs(level-50000));

S_500 = sigma(p500);

[~,plow] = min(abs(level-92500));

S_av = mean(sigma(plow:top));

% delta_p

[~,l2] = min(abs(level-(level(surf)-(level(surf)-level(top))/4)));

[~,l1] = min(abs(level-(level(surf)-3*(level(surf)-level(top))/4)));

delta_p_mid = abs(level(l1)-level(l2));

delta_p_av = abs(level(plow)-level(top))/2;

% deformation radius

L_mid = sqrt(S_mid/2)*delta_p_mid/mean(f);

L_500 = sqrt(S_500/2)*delta_p_mid/mean(f);

L_av = sqrt(S_av/2)*delta_p_av/mean(f);


end