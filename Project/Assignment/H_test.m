clear all; close all hidden; clc;

% wn = pi;
% k = 1;
% b = 0.005;
% q = tf('q');
% H = (wn^2*k)/(q^2 + 2*b*wn*q + wn^2);
% H = (1)/(q^2 + 1.7*q^1 + 0.95);
% %H = (q^-2 + 1.7*q^-1 + 0.95)/1;
% H = 1 - 0.5*q^-1 - 1*q^-2;
% H
% pzmap(H)
% 
% w_vect = linspace(0,pi,1000);
% Fresp = squeeze(freqresp(H,w_vect));

% f1 = figure();
% set(f1,'position',[900 550 450 400]);
% subplot(211)
% plot(w_vect,(abs(Fresp)));
% xlim([0 pi]);
% xlabel("Frequency rad/s");
% ylabel("Magnitude");
% grid on;
% subplot(212)
% plot(w_vect,angle(Fresp));
% xlim([0 pi]);
% xlabel("Frequency rad/s");
% ylabel("Phase");
% grid on;
% sgtitle("Bode plot H filter");

% f2 = figure();
% set(f2,'position',[1400 550 450 400]);
% bode(H);
% grid on;

N = 2^16;
r = zeros(1,N);

[u,y] = assignment_sys_25(r);

WINDOW = hann(N/256);
[Pxx,W] = pwelch(y,WINDOW);

% f3 = figure();
% set(f3,'position',[900 50 450 400]);
% plot(W,Pxx);
% grid on;
% xlabel('frequency');
% ylabel('Magnitude');
% title('Welch Power Spectral Density Estimate of H_0(q)');
% xlim([0 pi]);


f4 = figure();
semilogx(W,db(Pxx));
grid on;
ylabel('Magnitude [dB]');
title('Power Spectral Density H_0(q)');
xlim([0 pi]);
xlabel('Normalized frequency [0, \pi]')
xticks([0 0.1*pi 0.25*pi 0.5*pi 0.75*pi pi]);
xticklabels({'0','0.1\pi','0.25\pi','0.5\pi','0.75\pi','\pi'});
