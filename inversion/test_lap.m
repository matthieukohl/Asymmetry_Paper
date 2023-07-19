% test Lap function and the wavenumber metric

close all; clear;

N = 100;
L = 2*pi;
x = linspace(0,2*pi,N+1);
dx = abs(x(2)-x(1));
dy = abs(x(2)-x(1));
y = linspace(0,2*pi,N+1);
%[X,Y] = meshgrid(x,y);
[X,Y] = meshgrid(x(1:end-1),y(1:end-1));

f = sin(X).*sin(Y);
f(f>0) = sin(2*X(f>0)).*sin(2*Y(f>0));
f = f-mean(mean(f));

cmax = max(abs(f(:)));
nint = 20;
cint = cmax/nint;
contourf(f,-cmax:cint:cmax,'EdgeColor','none');
xlabel('x'); ylabel('y');
colorbar
colormap(redblue(nint));

%f = sin(X);

lap_f_theory = -20*f;

lap_f = Lap(N,N,dx,dy)*f(:);

lap_f = reshape(lap_f,N,N);

cmax = 2;
nint = 20;
cint = cmax/nint;

% figure('Renderer', 'painters', 'Position', [10 10 1200 400])
% 
% subplot(1,2,1)
% 
% contourf(lap_f_theory,-cmax:cint:cmax,'EdgeColor','none'); colorbar;
% colormap(redblue(nint));
% xlabel('x');
% ylabel('y');
% title('\rm Theory')
% caxis([-cmax cmax])
% 
% subplot(1,2,2)
% 
cmax = max(abs(lap_f(:)));
nint = 20;
cint = cmax/nint;

contourf(lap_f,-cmax:cint:cmax,'EdgeColor','none'); colorbar
colormap(redblue(nint));
xlabel('x');
ylabel('y')
title('\rm Numerical')
caxis([-cmax cmax])
% 
% figure
% 
% theory_diff = abs(lap_f-lap_f_theory);
% 
% contourf(theory_diff,'EdgeColor','none'); 
% title('\rm difference')
% colorbar
% xlabel('x');
% ylabel('y');
% 
% max(theory_diff(:))

% plot various wavenumbers

k = sqrt(0.5*mean(mean(abs(lap_f)))/(mean(mean(abs(f)))))
k_up = sqrt(0.5*mean(mean(abs(lap_f(f>0))))/(mean(mean(abs(f(f>0))))))
k_down = sqrt(0.5*mean(mean(abs(lap_f(f<0))))/(mean(mean(abs(f(f<0))))))

% ratio = lap_f./f;
% ratio(abs(ratio)>1000) = NaN;
% 
% k = sqrt(0.5*nanmean(nanmean(abs(ratio))))
% k_up = sqrt(0.5*nanmean(nanmean(abs(ratio(f>0)))))
% k_down = sqrt(0.5*nanmean(nanmean(abs(ratio(f<0)))))
