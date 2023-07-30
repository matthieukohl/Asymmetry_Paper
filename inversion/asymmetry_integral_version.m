function  lambda = asymmetry_integral_version(x,w)

% Calculates the Asymmetry Factor

w_av = trapz(x,w);
w_pr = w-w_av;

w_up = max(w,0);
w_up_av = trapz(x,w_up);
w_up_pr = w_up-w_up_av;

lambda = trapz(x,w_pr.*w_up_pr)/(trapz(x,w_pr.^2));

end