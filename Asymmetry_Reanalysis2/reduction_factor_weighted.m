function [r_mid,r_500,r_av] = reduction_factor_weighted(r,level,top,weight)

% Calculate the Deformation Radius to nondimensionalize Wavenumber

r = squeeze(sum(sum(weight.*r,1),2)./sum(sum(weight,1),2));

% Midpoint Method

% Static stability

surf = 1;
[~,I] = min(abs(level-(level(surf)-(level(surf)-level(top))/2)));

r_mid = r(I);

% 500

[~,p500] = min(abs(level-50000));

r_500 = r(p500);

% av

[~,plow] = min(abs(level-92500));

r_av = mean(r(plow:top));

end