%% System Identification assignment
clear all; close all; clc;

% Load data
r = 1:10000;
[u,y] = assignment_sys_25(r);

% Plotting data
figure(1)
hold on
plot(r,u)
plot(r,y)
hold off
grid on
title('Data')
legend('u','y')
xlabel('Sample [n]')
ylabel('Amplitude [-]')

%% FFT
ffty = fft(y,length(r))

figure(2)
semilogx(db(abs(r)),ffty)
grid on
xlim([0 length(r)/2])
ylim([-5000 5000])
xlabel('Frequency [Hz?? check Ts!]')
ylabel('Amplitude [dB]')
