function [r_mid,r_500,r_av] = reduction_factor_conditioned_on_ascent(r,omega,level,top)

% condition r only on ascent at 500hPa

[~,p500] = min(abs(level-50000));
rr = zeros(size(r,3),1);

for kk = 1:size(r,3)
r_xy = r(:,:,kk);
rr(kk) = squeeze(mean(mean(r_xy(omega(:,:,p500)<0))));
end

clear('r')
r = rr;
      
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