
function output = Skewness(field,lat)

% Calculate the vertically-averaged Skewness of Field (latitude,longitude,level,time) 
% using weighted lat/long averaging operators

weight = repmat(cosd(lat),1,size(field,2));
norm = sum(sum(weight));
[field_prime]= deal(zeros(size(field,1),size(field,2),size(field,3),size(field,4)));
[field_2,field_3,pholder]=deal(zeros(size(field,1),size(field,2)));

for t=1:size(field,4)
    for k=1:size(field,3)
        
     field_avg = sum(sum(squeeze(field(:,:,k,t)).*weight(:,:)))./norm';
     
     field_prime(:,:,k,t) = field(:,:,k,t)-field_avg;
  
     field_2(k,t) = sum(sum(squeeze(field_prime(:,:,k,t).^2).*weight(:,:)))./norm';
     
     field_3(k,t) = sum(sum(squeeze(field_prime(:,:,k,t).^3).*weight(:,:)))./norm';
     
     pholder(k,t) = -field_3(k,t)./field_2(k,t).^(3/2);
     
    end
end

output = nanmean(pholder,1);

end