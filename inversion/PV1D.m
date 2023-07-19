function dqdt = PV1D(t,q,adv,diab,level)
% Integrate the ODE based on Hoskins equation 74a
% d_t q = adv*dq/dp + diab*q

dq_dp = d_dp(q,level);

%dqdt = diab.*q;%

%dqdt = -adv.*dq_dp; %diab.*q;

dqdt = diab.*q - adv.*dq_dp;

end

