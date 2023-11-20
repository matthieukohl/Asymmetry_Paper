function output = Lambda1(field)

field_average = mean(mean(field(:,:,:,:),2),3);
field_prime = field-field_average;

field_up_prime = min(field_prime(:,:,:,:),0);

m = squeeze(mean(mean(field_prime(:,:,:,:).*field_up_prime(:,:,:,:),2),3)./(mean(mean(field_prime(:,:,:,:).^2,2),3)));
output = nanmean(m(:,:),1);


end

