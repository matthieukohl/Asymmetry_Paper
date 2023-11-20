% test wavenumber: 2/k = 1/ku + 1/kd

close all; clear

kd = 1;
ku = 2*kd;

wd = 1;
wu = wd*ku/kd;

xu = linspace(0,pi/ku,100);
xd = linspace(-pi/kd,0,100);

wuu = wu*sin(ku*xu);
wdd = wd*sin(kd*xd);

xx = [xd,xu(2:end)];
ww = [wdd,wuu(2:end)];

wuu_dd = -ku^2*wuu;
wdd_dd = -kd^2*wdd;
ww_dd = [wdd_dd,wuu_dd(2:end)];

figure
plot(xu,wu*sin(ku*xu)); hold on; 
plot(xd,wd*sin(kd*xd));

k = sqrt(mean(abs(ww_dd))/mean(abs(ww)))
k_up = sqrt(mean(abs(ww_dd(ww>0)))/mean(abs(ww(ww>0))))
k_down = sqrt(mean(abs(ww_dd(ww<0)))/mean(abs(ww(ww<0))))


x = linspace(0,2*pi,100);
f = sin(x)-0.8;
