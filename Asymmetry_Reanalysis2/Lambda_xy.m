function output = Lambda_xy(field)
% Calculate Lambda(lat,time) using xy-mean
% field(lon,lat,level,time)

field_average = mean(mean(field,1),2);
field_prime = field-field_average;

field_up = min(field,0);
field_up_average = mean(mean(field_up,1),2);
field_up_prime = field_up-field_up_average;

m = mean(mean(field_prime.*field_up_prime,1),2)./(mean(mean(field_prime.^2,1),2));

output = squeeze(m); % output(level);

end