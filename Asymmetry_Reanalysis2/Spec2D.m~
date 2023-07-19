function [modk, spectrum] = Spec2D(X, Y, U, V, k)
% 2DSpec: calculates the one-dimensional spectrum given two dimensional velocities at a height index k. Where z_k = k'th index of the vector of height, z.

% Calculation of Power Spectra

% Define wavenumbers in x,y.
M=size(X,1); %Largest allowed wavenumber is half the number of grid points, though, no?
N=size(Y,1);

M = M/2; 
N = N/2;

Lx=X(end);
Ly=Y(end);
Lmax = max([Lx Ly]);

k1=((1:1:M)-(M/2))';
k2=((1:1:N)-(N/2))';
%k1 = (Lmax/Lx)* k1;
%k2 = (Lmax/Ly)* k2;
k1 = 2*pi/Lx*k1;
k2 = 2*pi/Ly*k2;

% Fourier Transform the fields:
uk = fftshift(fft2(U(:,:,k)));
vk = fftshift(fft2(V(:,:,k)));

% Calculates E(\vec{k})
E = abs(uk).^2+abs(vk).^2;

% Calculate one-dimensional (isotropic) spectrum
maxk=round(sqrt((k1(end)^2+k2(end)^2)));
spectrum=zeros(maxk,1); 
modk=zeros(maxk,1);

for i=1:1:M
    for j=1:1:N
        wavenumber=round(sqrt((k1(i)^2+k2(j)^2)));
        if wavenumber==0
           wavenumber=wavenumber+1;
           spectrum(wavenumber)=spectrum(wavenumber)+E(i,j);
           modk(wavenumber)=wavenumber; 
        else
           spectrum(wavenumber)=spectrum(wavenumber)+E(i,j);
           modk(wavenumber)=wavenumber;
        end
    end
end

end
