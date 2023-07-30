function output = Lambda_weighted_vs_p(field,lat)
% Calculate vertically-averaged Asymmetry Factor lambda Following O'Gorman 2011 for Field(lat,long,level,time) 
% using weighted lat/long averages

weight = repmat(cosd(lat),1,size(field,2),size(field,3),size(field,4));
norm = sum(sum(weight(:,:,1,1)));

field_up = min(field(:,:,:,:),0);
        
field_avg = sum(sum(field.*weight,1),2)/norm;
field_up_avg = sum(sum(field_up.*weight,1),2)/norm;

field_prime = field-field_avg;
field_up_prime = field_up-field_up_avg;

field_cov = sum(sum(field_prime.*field_up_prime.*weight,1),2)/norm;
field_2 = sum(sum(field_prime.^(2).*weight,1),2)/norm;

pholder = field_cov./field_2;

output = nanmean(squeeze(pholder),2);

end