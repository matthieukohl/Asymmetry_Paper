clear;
latmax = 90; % degrees
latmin = -90;
lonmax = 360; % degrees
lonmin = 0;

r=0.1;

for i=1:3
    


path_temp = strcat('/disk7/mkohl/interpolation/output1/temp_day00',num2str(25 * i),'h00.nc');
path_u = strcat('/disk7/mkohl/interpolation/output1/u_day00',num2str(25 * i),'h00.nc');
path_v = strcat('/disk7/mkohl/interpolation/output1/v_day00',num2str(25 * i),'h00.nc');
path_omega = strcat('/disk7/mkohl/interpolation/output1/omega_day00',num2str(25 * i),'h00.nc');
path_vor = strcat('/disk7/mkohl/interpolation/output1/vor_day00',num2str(25 * i),'h00.nc');

lat_series_original = ncread(path_temp, 'latitude');
lon_series_original = ncread(path_temp, 'longitude');
[lat_indices, lon_indices] = latlonindices(lat_series_original, lon_series_original, latmin, latmax, lonmin, lonmax);
plevels = double(ncread(path_temp, 'level'));

Omega = 7.2921e-5;
R = 6371000.0;
Ra = 287.04;
cp = 1005.7;
kappa = Ra/cp;
dt = 3600.0 * 6.0;
window_x = 1;
window_y = 1;

lon_series = lon_series_original(lon_indices);
lat_series = lat_series_original(lat_indices);
dphi = (lat_series(2) - lat_series(1)) / 180.0 * 3.1415926;
dlambda = (lon_series(2) - lon_series(1)) / 180.0 * 3.1415926;

event_latspan = [83:111];
event_lonspan = [1:256];
event_timespan = [1:100];
level_indices = [5:30];
level = plevels(level_indices);
p_up = 20000;
p_w = 5000;
p_up_index = find(level==p_up);
Reduction = reduction_factor(r,level);
lat = lat_series_original(event_latspan);
lon = lon_series_original(event_lonspan);
phi = lat / 180.0 * 3.1415926;
lambda = lon / 180.0 * 3.1415926;
[~, Phi] = meshgrid(lambda, phi);

start = [event_lonspan(1), event_latspan(1), level_indices(1), event_timespan(1)];
count = [length(event_lonspan), length(event_latspan), length(level_indices), length(event_timespan)];

T_in= ncread(path_temp,'temp_plevel', start, count);
u_in= ncread(path_u,'temp_plevel', start, count);
v_in= ncread(path_v,'temp_plevel', start, count);
omega_in= ncread(path_omega,'temp_plevel', start, count);
%vor_in= ncread(path_vor,'vor_plevel', start, count);

[T, u, v, omega,vor] = ...
        deal(zeros(length(event_latspan), length(event_lonspan), length(level_indices), length(event_timespan)));
 
    for k = 1 : length(level_indices)
        for t = 1 : length(event_timespan)
            T(:, :, k, t) = T_in(:, :, k, t)';
            u(:, :, k, t) = u_in(:, :, k, t)';
            v(:, :, k, t) = v_in(:, :, k, t)';
            omega(:, :, k, t) = omega_in(:, :, k, t)';
            %vor(:, :, k, t) = vor_in(:, :, k, t)';
        end
    end

f0 = 2 * Omega * sin(lat_series_original(97) / 180 * 3.1415926);
beta = 2 * Omega / R * cos(lat_series_original(97) / 180 * 3.1415926);
[A, B, C, J, J1, J2, Qx, Qy, ...
    du_dlambda, dv_dlambda, du_dphi, dv_dphi, ...
    dT_dlambda, dT_dphi, sigma_accu,sigma_accu_smoothed,sigma_r,sigma_r_smoothed] = deal(zeros(length(event_latspan), length(event_lonspan), length(level), length(event_timespan)));
    
[sigma, sigma_eff,dtheta_dp_ma, T_avg,theta_avg] = deal(zeros(length(level), length(event_timespan)));
[rhs, rhs_smoothed, omega_QG, omega_b, sigma_old,sigma_new,sigma_new_smoothed] = deal(zeros(length(event_latspan), length(event_lonspan), length(level), length(event_timespan)));

for t = 1 : length(event_timespan)

T_avg(:, t) = reshape(mean(mean(T(:, :, :, t), 1), 2), size(level));



% % effective static stability: sigma_eff
% temp_omega = squeeze(omega(:, :, :, t));
% [~, ~, theta_avg] = eff_stat_stab(level, T_avg(:, t), 1); % this routine is just used to calculate potential temperature                


% spatially varying static stability: sigma_accu
% this quantity is mainly used to calculate J term
Level = repmat(reshape(level, 1, 1, length(level)), length(event_latspan), length(event_lonspan), 1);
theta = T(:, :, :, t) .* (1e5 ./ Level) .^ kappa;
sigma_accu(:, :, :, t) = -Ra * T(:, :, :, t) ./ (Level .* theta) .* d_dp(theta, level);

% static stability: sigma
% the following formula is more accurate since it takes into account the moist part in the exponent 
% when calculating the potential temperature
theta_avg (:,t) = reshape(mean(mean(theta, 1), 2), size(level));
sigma(:, t) =  - Ra * T_avg(:, t) ./ (level .* theta_avg(:,t)) .* gradient(theta_avg(:,t), level);


%deal with points where sigma_accu is not smooth
sigma_smooth_window = 3;    
if sigma_smooth_window ~= 0
    for k = 1 : length(level)
        sigma_accu_smoothed(:, :, k, t) = smooth2a(sigma_accu(:, :, k, t), sigma_smooth_window, sigma_smooth_window);
    end
else
    sigma_accu_smoothed(:, :, :, t) = sigma_accu(:, :, :, t);
end


% term A (Q vector)

for k = 1 : length(level)

    % probably the ug and vg field need smoothing as well...
    du_dlambda(:, :, k, t) = d_dlambda(u(:, :, k, t), dlambda);
    dv_dlambda(:, :, k, t) = d_dlambda(v(:, :, k, t), dlambda);
    du_dphi(:, :, k, t) = d_dphi(u(:, :, k, t), dphi);
    dv_dphi(:, :, k, t) = d_dphi(v(:, :, k, t), dphi);

    dT_dlambda(:, :, k, t) = d_dlambda(T(:, :, k, t), dlambda);
    dT_dphi(:, :, k, t) = d_dphi(T(:, :, k, t), dphi);
    Qx(:, :, k, t) = 1 / R^2 * (...
        (1 ./ cos(Phi).^2) .* du_dlambda(:, :, k, t) .* dT_dlambda(:, :, k, t) + ...
        (1 ./ cos(Phi))    .* dv_dlambda(:, :, k, t) .* dT_dphi(:, :, k, t));
    Qy(:, :, k, t) = 1 / R^2 * (...
        (1 ./ cos(Phi))    .* du_dphi(:, :, k, t)    .* dT_dlambda(:, :, k, t) + ...
                              dv_dphi(:, :, k, t)    .* dT_dphi(:, :, k, t));
    div_temp = div(Qx(:, :, k, t), Qy(:, :, k, t), Phi, dphi, dlambda);
    A(:, :, k, t) = 2 * Ra / level(k) * div_temp;
end
clear('temp_x', 'temp_y', 'div_temp');

% term B (planetary vorticity)

temp = d_dp(v(:, :, :, t), level);
for j = 1 : length(lat)
    B(j, :, :, t) = f0 * beta * temp(j, :, :);
end
clear('temp');

% right hand side

rhs(:, :, :, t) = A(:, :, :, t) + B(:, :, :, t); 


% % deal with points where rhs is not smooth
% rhs_smooth_window = 2;    
% if rhs_smooth_window ~= 0
%     for k = 1 : length(level)
%         rhs_smoothed(:, :, k, t) = smooth2a(rhs(:, :, k, t), rhs_smooth_window, rhs_smooth_window);
%     end
% else
%     rhs_smoothed(:, :, :, t) = rhs(:, :, :, t);
% end


%Inversion

%omega_QG(:, :, :, t) = SIP_sigma(R, f0, length(lon), length(phi), length(level), phi, lambda, level, ...
                                   %sigma_accu_smoothed(:, :, :, t), dphi, dlambda, rhs(:, :, :, t),  omega(:,:,:,t));

omega_QG(:, :, :, t) = SIP_diagnostic(R, f0, length(lon), length(phi), length(level), phi, lambda, level, ...
                                                            sigma(:, t), dphi, dlambda, rhs(:,:,:,t),  omega_b(:,:,:,t));

                                              
for l=1:length(event_latspan)
    for m=1:length(event_lonspan)
      for n= 1:length(level)
          
         sigma_old(l, m, n, t) = 1*sigma(n,t); 
        
            
      end
    end
end

for s=1:10
       
 for l=1:length(event_latspan)
    for m=1:length(event_lonspan)
      for n= 1:length(level)
          if omega_QG(l,m,n,t)>=0
           sigma_new(l, m, n, t) = 1*sigma(n,t); %1*sigma_accu(l,m,n,t); % %
          else 
           sigma_new(l, m, n, t) = Reduction(n)*sigma(n,t); % Reduction(n)*sigma_accu(l,m,n,t); 
          end 
      end
    end
 end

 %Take Arithmetic Mean of Old and New Field for smoother transition
 
 sigma_new = 0.5*(sigma_new+sigma_old);
 
 

 
%deal with points where sigma_accu is not smooth
sigma_smooth_window = 2;    
if sigma_smooth_window ~= 0
    for k = 1 : length(level)
        sigma_new_smoothed(:, :, k, t) = smooth2a(sigma_new(:, :, k, t), sigma_smooth_window, sigma_smooth_window);
    end
else
    sigma_new_smoothed(:, :, :, t) = sigma_new(:, :, :, t);
end

for m=1:3
    
if sigma_smooth_window ~= 0
    for k = 1 : length(level)
        sigma_new_smoothed(:, :, k, t) = smooth2a(sigma_new_smoothed(:, :, k, t), sigma_smooth_window, sigma_smooth_window);
    end
else
    sigma_new_smoothed(:, :, :, t) = sigma_new_smoothed(:, :, :, t);
end
end

%  for n=1:length(level)
%      
%      
%     sigma_new_smoothed(:,:,n,t)= imgaussfilt(sigma_new(:,:,n,t),2);
%     
%     for m=1:2
%         
%     sigma_new_smoothed(:,:,n,t)= imgaussfilt(sigma_new_smoothed(:,:,n,t),1);
%     
%     end
%     
%  end
                               
omega_QG(:, :, :, t) = SIP_sigma(R, f0, length(lon), length(phi), length(level), phi, lambda, level, ...
                                   sigma_new_smoothed(:, :, :, t), dphi, dlambda, rhs(:, :, :, t),  omega_b(:,:,:,t));   

                               
sigma_old = sigma_new;
                               
 end 


                                                        
% for s=1:1
%        
%  for l=1:length(event_latspan)
%     for m=1:length(event_lonspan)
%       for n= 1:length(level)
%           if omega_QG(l,m,n,t)>=0
%            sigma_r(l, m, n, t) = 1*sigma(n,t); %1*sigma_accu(l,m,n,t); % %
%           else 
%            sigma_r(l, m, n, t) = Reduction(n)*sigma(n,t); % Reduction(n)*sigma_accu(l,m,n,t); 
%           end 
%       end
%     end
%  end
% 
% 
%  for n=1:length(level)
%      
%      
%     sigma_r_smoothedd(:,:,n,t)= imgaussfilt(sigma_r(:,:,n,t),1);
%     
%     for i=1:1
%         
%     sigma_r_smoothedd(:,:,n,t)= imgaussfilt(sigma_r_smoothedd(:,:,n,t),1);
%     
%     end
%     
%  end
% 
% 
%  
% % windowWidth = 9; % Some odd number.
% % kernel = ones(1, windowWidth) / windowWidth;
% % for t=1:length(event_timespan)
% %  for l=1:length(event_latspan)
% %     for n=1:length(level)
% %         
% %         sigma_r_smoothedd(l,:,n,t) = conv(sigma_r(l,:,n,t), kernel, 'same');
% %         
% %     end 
% %  end
% % end
%  
% % %deal with points where sigma_accu is not smooth
% % sigma_smooth_window = 3;    
% % if sigma_smooth_window ~= 0
% %     for k = 1 : length(level)
% %         sigma_r_smoothed(:, :, k, t) = smooth2a(sigma_r(:, :, k, t), sigma_smooth_window, sigma_smooth_window);
% %     end
% % else
% %     sigma_r_smoothed(:, :, :, t) = sigma_r(:, :, :, t);
% % end
%  
% % %deal with points where sigma_accu is not smooth
% % sigma_smooth_window = 3;    
% % if sigma_smooth_window ~= 0
% %     for k = 1 : length(level)
% %         sigma_r_smoothed(:, :, k, t) = smooth2a(sigma_r(:, :, k, t), sigma_smooth_window, sigma_smooth_window);
% %     end
% % else
% %     sigma_r_smoothed(:, :, :, t) = sigma_r(:, :, :, t);
% % end
% 
% % test repeated smoothing vs. larger window
% 
% % sigma_smooth_window = 4;    
% % if sigma_smooth_window ~= 0
% %     for k = 1 : length(level)
% %         sigma_r_smoothedd(:, :, k, t) = smooth2a(sigma_r(:, :, k, t), sigma_smooth_window, sigma_smooth_window);
% %     end
% % else
% %     sigma_r_smoothedd(:, :, :, t) = sigma_r(:, :, :, t);
% % end
% % 
% % 
% % for i=1:1
% %     
% % if sigma_smooth_window ~= 0
% %     for k = 1 : length(level)
% %         sigma_r_smoothedd(:, :, k, t) = smooth2a(sigma_r_smoothedd(:, :, k, t), sigma_smooth_window, sigma_smooth_window);
% %     end
% % else
% %     sigma_r_smoothedd(:, :, :, t) = sigma_r_smoothedd(:, :, :, t);
% % end
% % 
% % end
% 
% % omega_QG(:, :, :, t) = SIP_sigma(R, f0, length(lon), length(phi), length(level), phi, lambda, level, ...
% %                                    sigma_r_smoothed(:, :, :, t), dphi, dlambda, rhs(:, :, :, t),  omega(:,:,:,t));   
%                                
% omega_QG(:, :, :, t) = SIP_sigma(R, f0, length(lon), length(phi), length(level), phi, lambda, level, ...
%                                    sigma_r_smoothedd(:, :, :, t), dphi, dlambda, rhs(:, :, :, t),  omega(:,:,:,t));   
% 
% 
% figure(1)
% plot(omega_QG(15,:,12,1)); hold on;
% 
% figure(2)
% plot(sigma_r(15,:,12,1)); hold on;
% plot(sigma_r_smoothedd(15,:,12,1));
%                                
% end 
                               

% f_new(:,:,:,t) = omega_QG(:,:,:,t);
% 
% M= abs(f_new(:,:,:,t)-f_old(:,:,:,t));
% 
% V= mean(M(:));
% 
% count = 0;
% 
% while V>0.0001 && count<15
%  count = count +1;                              
%  for l=1:length(event_latspan)
%     for m=1:length(event_lonspan)
%       for n= 1:length(level) 
%           if f_new(l,m,n,t)>=0
%            sigma_r(l, m, n, t) = 1*sigma_accu(l,m,n,t);
%           else 
%            sigma_r(l, m, n, t) =  Reduction(n)*sigma_accu(l,m,n,t); 
%           end 
%       end
%     end
%  end
% 
% % deal with points where sigma_accu is not smooth
% sigma_smooth_window = 5;    
% if sigma_smooth_window ~= 0
%     for k = 1 : length(level)
%         sigma_r_smoothed(:, :, k, t) = smooth2a(sigma_r(:, :, k, t), sigma_smooth_window, sigma_smooth_window);
%     end
% else
%     sigma_r_smoothed(:, :, :, t) = sigma_r(:, :, :, t);
% end
%  
% f_old(:,:,:,t) = f_new(:,:,:,t);
% 
% f_new(:, :, :, t) = SIP_sigma(R, f0, length(lon), length(phi), length(level), phi, lambda, level, ...
%                                    sigma_r_smoothed(:, :, :, t), dphi, dlambda, rhs(:, :, :, t),  omega(:,:,:,t));  
% 
% M = abs(f_new(:,:,:,t)-f_old(:,:,:,t));
% 
% V = mean(M(:));
% 
% figure(7)
% contour(f_new(:,:,12,1)); hold on
% colorbar
% 
% figure(8)
% contour(omega(:,:,12,1)); hold on
% colorbar
%                                
% end 
% 
% omega_QG(:,:,:,t)=f_new(:,:,:,t);

end 

save(strcat('output/rescale_0.1/sigma/relax_mean/smooth_cross/bc_zero/Total_rhs_smooth0_sigma_smooth2.3_iter10_day',num2str(25*i),'.mat'),'omega','omega_QG','rhs','sigma','sigma_new','sigma_new_smoothed');


end
