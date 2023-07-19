function  lambda = asymmetry(w)

% Calculates the Asymmetry Factor

w_av = mean(w);
w_pr = w-w_av;

w_up = max(w,0);
w_up_av = mean(w_up);
w_up_pr = w_up-w_up_av;

lambda = mean(w_pr.*w_up_pr)/(mean(w_pr.^2));

end
