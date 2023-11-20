function [var_out] = rearrange_era_2d(var_in,latspan,lonspan)

% Rearrange the field var such that var(lat,lon,k) starting from low to
% high latitude and surface to top


[var_out] = ...
 deal(zeros(length(latspan), length(lonspan)));




var_out(:,:) = var_in(:,:)'; 



% flip lat

for ll =1:length(lonspan)
      
var_out(:,ll) = flip(var_out(:,ll));
  
end


    


end