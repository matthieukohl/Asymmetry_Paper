function output = Lambda_zonal(field)
% Calculate Lambda using zonal-mean, averaged over latitude and height

field_average = mean(field,2);
field_average = repmat(field_average,1,size(field,2));
field_prime = field-field_average;

field_up = min(field,0);
field_up_average = mean(field_up,2);
field_up_average = repmat(field_up_average,1,size(field,2));
field_up_prime = field_up-field_up_average;

m = squeeze(mean(field_prime.*field_up_prime,2)./(mean(field_prime.^2,2)));

%output = mean(m,3);
output = mean(mean(m));

end