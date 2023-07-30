function [var_out] = rearrange_era(var_in,latspan,lonspan,level)

% Rearrange the field var such that var(lat,lon,k) starting from low to
% high latitude and surface to top


[var_out] = ...
 deal(zeros(length(latspan), length(lonspan), length(level)));



for k = 1 : length(level)
    var_out(:,:,k) = var_in(:,:,k)'; 
end


% flip level

for ii = 1:length(latspan)
    for ll =1:length(lonspan)
        
            var_out(ii,ll,:) = flip(var_out(ii,ll,:));
            
       
    end
end

% flip lat

for ll =1:length(lonspan)
      
 for kk = 1:length(level)
 var_out(:,ll,kk) = flip(var_out(:,ll,kk));
 end
       
end


    


end