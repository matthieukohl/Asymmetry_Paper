% test sinusoid
wd = 1;
kd = 0.5;
ku = 2;
wu = ku*wd/kd; % constraint by mass conservation

dx = 1e-3;
Nd = pi/(kd*dx)+1;
Nd = round(Nd);
Nu = pi/(ku*dx)+1;
Nu = round(Nu);

% define the domain
L = pi/kd + pi/ku;
xd = linspace(-pi/kd,0,Nd);
xu = linspace(0,pi/ku,Nu);

x = [xd,xu(2:end)];
w = [wd*sin(kd*xd),wu*sin(ku*xu(2:end))];

% lambda 

%asymmetry(w)

w_p = w;
w_up = [zeros(1,Nd),wu*sin(ku*xu(2:end))];
w_up_p = w_up - 1/L*trapz(x,w_up);

num = 1/L*trapz(x,w_p.*w_up_p);
den = 1/L * trapz(x,w_p.^2);

num/den

% prediction 

Lu = pi/ku;
Ld = pi/kd;
1 - Lu/(Lu+Ld)

kd/(kd^2/ku+kd)