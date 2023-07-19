function [var_out] = rearrange_single_level(var_in,latspan,lonspan,timespan)

% Rearrange the field var such that var(lat,lon,k,t) starting from low to
% high latitude and surface to top

sz = size(var_in);
ndim = length(sz);

if ndim ==3

[var_out] = ...
 deal(zeros(length(latspan), length(lonspan), length(timespan)));


for t = 1 : length(timespan)
   
        var_out(:,:,t) = var_in(:,:,t)'; 
  
end



% flip lat

for ll =1:length(lonspan)
        for tt=1:length(timespan)
            
            var_out(:,ll,tt) = flip(var_out(:,ll,tt));
            
        end
end

else
    
[var_out] = ...
 deal(zeros(length(latspan), length(lonspan), length(timespan)));


for t = 1 : length(timespan)
 var_out(:,:,t) = var_in(:,:,t)'; 
end

% flip lat

for ll =1:length(lonspan)
    for tt=1:length(timespan)
        var_out(:,ll,tt) = flip(var_out(:,ll,tt));
    end
end
 
end

end