function output = Lambda_x(field)
% Calculate Lambda(lat,time) using x-mean
% field(lon,lat,level,time)

field_average = mean(field,1);
field_prime = field-field_average;

field_up = min(field,0);
field_up_average = mean(field_up,1);
field_up_prime = field_up-field_up_average;

m = mean(field_prime.*field_up_prime,1)./(mean(field_prime.^2,1));


output = squeeze(m); % output(lat,level);

end