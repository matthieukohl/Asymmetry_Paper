% calculat the toy model envelope

r = logspace(-7,0,10);
k = linspace(1,30,10);
ind_envelope = zeros(size(r));
k_envelope = zeros(size(r));
lambda_envelope = zeros(size(r));
lambda_envelope_high_rossby = zeros(size(r));

lambda = zeros(length(r),length(k));

for ii = 1:length(r)
    for jj = 1:length(k)
        
lambda(ii,jj) = toy_model_nondimensional(r(ii),k(jj));
    end
end

for ii = 1:length(r)
    
lambda_diff = diff(lambda(ii,:));
ind = find(lambda_diff<0.001,1,'first');

ind_envelope(ii) = ind;
k_envelope(ii) = k(ind);

end

% make plot of the envelope at low and high rossby

figure

for ii = 1:length(r)
lambda_envelope(ii) = lambda(ii,ind_envelope(ii));
end

semilogx(r,lambda_envelope,'r'); hold on;
xlabel('Reduction factor r')
ylabel('\lambda')
title('\rm Toy-Model \lambda-Envelope')

eps = 10;

for ii = 1:length(r)
lambda_envelope_high_rossby(ii) = toy_model_nondimensional(r(ii),k_envelope(ii)/sqrt(1+eps));
end

semilogx(r,lambda_envelope_high_rossby,'b');
xlabel('r')
ylabel('\lambda')
ylim([0.5 1])