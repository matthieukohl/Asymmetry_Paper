function output = Skewness(field)

% Calculate the vertically-averaged Skewness of Field (latitude,longitude,level,time) 
% using lat/long averaging operators

% field_average = mean(mean(field(:,:,:,:),1),2);
% 
% field_average = repmat(field_average,size(field,1),size(field,2));
% 
% field_prime = field - field_average;
% 
% m = squeeze(mean(mean(field_prime(:,:,:,:).^3,1),2)./(mean(mean(field_prime(:,:,:,:).^2,1),2)).^(3/2));
% output = -nanmean(m(:,:),1);
% %output = -m(12,:);

field_average = nanmean(nanmean(field(:,:,:,:),1),2);

field_average = repmat(field_average,size(field,1),size(field,2));

field_prime = field - field_average;

m = squeeze(nanmean(nanmean(field_prime(:,:,:,:).^3,1),2)./(nanmean(nanmean(field_prime(:,:,:,:).^2,1),2)).^(3/2));
output = -nanmean(m(:,:),1);
%output = -m(12,:);

end



 