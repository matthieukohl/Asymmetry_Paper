function [var_out] = rearrange(var_in,latspan,lonspan,level,timespan)

% Rearrange the field var such that var(lat,lon,k,t) starting from low to
% high latitude and surface to top

sz = size(var_in);
ndim = length(sz);

if ndim ==4

[var_out] = ...
 deal(zeros(length(latspan), length(lonspan), length(level), length(timespan)));


for t = 1 : length(timespan)
    for k = 1 : length(level)
        var_out(:,:,k,t) = var_in(:,:,k,t)'; 
    end
end

% flip level

for ii = 1:length(latspan)
    for ll =1:length(lonspan)
        for tt=1:length(timespan)
            var_out(ii,ll,:,tt) = flip(var_out(ii,ll,:,tt));
            
        end
    end
end

% flip lat

for ll =1:length(lonspan)
        for tt=1:length(timespan)
            for kk = 1:length(level)
            var_out(:,ll,kk,tt) = flip(var_out(:,ll,kk,tt));
            end
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

