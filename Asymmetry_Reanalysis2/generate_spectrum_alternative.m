function [x] = generate_spectrum_alternative(spec,L)
% generate a spatial field with prescribed spectral density
% generate a prescribed power spectrum

N = length(spec);
dx = L/N;
k = 2*pi/L*[0:N/2,-N/2+1:-1];
P = spec;
P = P';

A = sqrt(P/(2*N*dx)) .* exp( 2*pi*1i * rand(N,1) );

% For even-length cases, the first sample has no
% equivalent positive frequency, so we force it
% to be real too.
A(1) = abs( A(1) );

% DC sample must be real.
n0 = floor( N/2 ) + 1;
A(n0) = P(n0);



% Force the spectrum to be conjugate symetric.
if( mod(N,2) == 0 )
  A(2:n0-1) = flipud( conj(A(n0+1:end)) );
else
  A(1:n0-1) = flipud( conj(A(n0+1:end)) );
end

% Time domain data.
x = fft(A);
x = x-mean(x);

% test
power = (2*dx/N)*abs(fft(x)).^2;

% figure
% loglog(k,power,'b'); hold on;
% loglog(k,spec,'k--');

end

