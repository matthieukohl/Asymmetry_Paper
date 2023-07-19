function d_d2pz = d_d2p(z, level)

sizez = size(z);
Nk = length(level);
d_d2pz = zeros(sizez);
level = double(level);

if (length(sizez) == 2) && ((sizez(1) == 1) ||(sizez(2) == 1))
    
% upper and lower boundary conditions
dp1 = level(2)-level(1); dp2 = level(3)-level(2);

d_d2pz(1)  = (z(1)*dp2+z(3)*dp1-z(2)*(dp1+dp2))/(dp1^2*dp2);

dpN = level(Nk)-level(Nk-1); dpNm1 = level(Nk-1)-level(Nk-2);

d_d2pz(Nk) = (z(Nk)*dpNm1+z(Nk-2)*dpN - z(Nk - 1)*(dpN+dpNm1)) / (dpN^2*dpNm1);

% inner points

for k = 2 : Nk - 1
    dp1 = level(k + 1) - level(k);
    dp2 = level(k) - level(k - 1);
    d_d2pz(k) = (z(k+1)*dp2+z(k-1)*dp1-z(k)*(dp1+dp2))/(0.5*dp1^2*dp2+0.5*dp1*dp2^2);
end   

else if length(sizez) == 3
   
% upper and lower boundary conditions
dp1 = level(2)-level(1); dp2 = level(3)-level(2);

d_d2pz(: ,:, 1)  = (z(:,:,1)*dp2+z(:,:,3)*dp1-z(:,:,2)*(dp1+dp2))/(dp1^2*dp2);

dpN = level(Nk)-level(Nk-1); dpNm1 = level(Nk-1)-level(Nk-2);

d_d2pz(:,:, Nk) = (z(:, :,Nk)*dpNm1+z(:,:,Nk-2)*dpN - z(:, :, Nk - 1)*(dpN+dpNm1)) / (dpN^2*dpNm1);

% inner points

for k = 2 : Nk - 1
    dp1 = level(k + 1) - level(k);
    dp2 = level(k) - level(k - 1);
    d_d2pz(: ,:, k) = (z(:,:,k+1)*dp2+z(:,:,k-1)*dp1-z(:,:,k)*(dp1+dp2))/(0.5*dp1^2*dp2+0.5*dp1*dp2^2);
end

else

        disp('error with dimension in function d_dp()');
        
    end
end
        
end