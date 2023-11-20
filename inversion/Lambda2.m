function output = Lambda2(field)
%% Calculate vertically-averaged Asymmetry Factor lambda Following O'Gorman 2011 for Field(lat,long,level,time) using lat/long averages


%% Calculate vertically-averaged Asymmetry Factor lambda Following O'Gorman 2011 for Field(lat,long,level,time) using lat/long averages

field_average = mean(mean(field(:,:,:,:),1),2);
field_average = repmat(field_average,size(field,1),size(field,2));
field_prime = field-field_average;

field_up = min(field(:,:,:,:),0);
field_up_average = mean(mean(field_up(:,:,:,:),1),2);
field_up_average = repmat(field_up_average,size(field,1),size(field,2));
field_up_prime = field_up-field_up_average;

output = mean(mean(field_prime(:,:,:,:).*field_up_prime(:,:,:,:),1),2)./(mean(mean(field_prime(:,:,:,:).^2,1),2));

end

