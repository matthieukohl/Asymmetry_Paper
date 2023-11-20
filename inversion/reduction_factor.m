function R = reduction_factor(r,level)

%% Reduction factor following O'Gorman 2018
p_up = 20000;
p_w = 5000;

R=ones(size(level));

for l=1:size(level)

R(l) = r+0.5*(1-r)*(1-tanh((level(l)-p_up)/p_w));

end 

end

