function [output] = weighted_avg(field,lat)

% weighted horizontal average

weight = repmat(cosd(lat),1,size(field,2));
norm = sum(sum(weight));
[output] = deal(zeros(size(field,3),size(field,4)));


for k=1:size(field,3)
   for t=1:size(field,4)
    
        
     output(k,t) = sum(sum(field(:,:,k,t).*weight(:,:),1),2)./norm';
     
   end
end

output = mean(output(:,:),1);

end

