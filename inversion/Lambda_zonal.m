function output = Lambda_zonal(field)
%% Calculate vertically-averaged Asymmetry Factor lambda Following O'Gorman 2011 for Field(lat,long,level,time) using lat/long averages

field_average = mean(field(:,:,:,:),2);
%field_average = repmat(field_average,size(field,1),size(field,2));
field_prime = field-field_average;

field_up = min(field(:,:,:,:),0);
field_up_average = mean(field_up(:,:,:,:),2);
%field_up_average = repmat(field_up_average,size(field,1),size(field,2));
field_up_prime = field_up-field_up_average;

m = squeeze(mean(field_prime(:,:,:,:).*field_up_prime(:,:,:,:),2)./(mean(field_prime(:,:,:,:).^2,2)));
output = squeeze(nanmean(nanmean(m(:,:,:),1),2));

%output = m(12,:);

end
