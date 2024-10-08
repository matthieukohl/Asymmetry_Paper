function [S_mid,S_500,S_av] = Static_Stability(sigma,level,top)

% Calculate the Deformation Radius to nondimensionalize Wavenumber

sigma = squeeze(mean(mean(sigma,1),2));

% Midpoint Method

% Static stability

surf = 1;
[~,I] = min(abs(level-(level(surf)-(level(surf)-level(top))/2)));

S_mid = sigma(I);

[~,p500] = min(abs(level-50000));

S_500 = sigma(p500);

[~,plow] = min(abs(level-92500));

S_av = mean(sigma(plow:top));
end