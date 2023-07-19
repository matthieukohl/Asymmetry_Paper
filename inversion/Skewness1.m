function output = Skewness1(field)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

field_average = mean(mean(field(:,:,:,:),2),3);
field_prime = field - field_average;


m = squeeze(-mean(mean(field_prime(:,:,:,:).^3,2),3)./(mean(mean(field_prime(:,:,:,:).^2,2),3)).^(3/2));
output = nanmean(m(:,:),1);


end

