clear, clc, close all

% noise parameters
fs = 44100;         % sampling frequency, Hz
T = 60;             % signal duration, s
N = round(fs*T);    % number of samples
t = (0:N-1)/fs;     % time vector
alpha = -3;         % PSD slope of -30 dB/dec i.e., black noise

% noise generation
x = arbssnoise(N, alpha);

% calculate the noise PSD
winlen = round(5*fs);
noverlap = round(winlen/5);
nfft = winlen;
win = hanning(winlen, 'periodic');
[Pxx, f] = pwelch(x, win, noverlap, nfft, fs, 'onesided');
Pxx = 10*log10(Pxx);

% plot the noise in the time domain
figure(1)
subplot(2, 1, 1)
plot(t, x, 'r')
grid on
xlim([0 max(t)])
ylim([-1.1*max(abs(x)) 1.1*max(abs(x))])
set(gca, 'FontName', 'Times New Roman', 'FontSize', 14)
xlabel('Time, s')
ylabel('Amplitude')
title('The generated noise signal in the time domain')

% plot the noise PSD
subplot(2, 1, 2)
semilogx(f, Pxx, 'r', 'LineWidth', 2)
hold on
semilogx(f, Pxx(f==1)+10*log10(f.^alpha), 'b', 'LineWidth', 2)
grid on
xlim([1 max(f)])
set(gca, 'FontName', 'Times New Roman', 'FontSize', 14)
xlabel('Frequency, Hz')
ylabel('Magnitude, dBV{^2}/Hz')
title('PSD of the generated noise signal')
legend('Empirical PSD', 'Theoretical PSD')