clear all; close all hidden; clc;

N = 3000;
n = (0:N-1);
Ts = 1/(N);
t = n*Ts;

%PRBS1
%r = zeros(1,N);
% r(1) = (randi([0 1])-0.5)*4;
% p = 0.5;
% minclock = 2;
% clock = 0;
% for i = 2:N
%     if rand() >= p && mod(clock,minclock) == 0
%         r(i) = -r(i-1);
%         clock = 0;
%     else
%         r(i) = r(i-1);
%     end
%     clock = clock + 1;
% end

%PRBS2
Nc = 2;
R = 1;
white = randn(N/(Nc*R),1);
r = repmat(2*sign(white(ceil([1:(N/R)]/Nc))),[R,1]);

%PRBS3
%B = 1/4
%r = idinput(N,"PRBS",[0 B],[-2 2]);
%r = zeros(1,N);

[u,y] = assignment_sys_25(r);
[u2,y2] = assignment_sys_25(r);


AutoCorrelate(r,'r',1,25);
AutoCorrelate(u,'u',1,25);

disp("Energy r:"+sum(r.^2));
disp("Energy u:"+sum(u.^2));

WINDOW = hann(60);
[Pr,~] = pwelch(r,WINDOW);
[Pu,w] = pwelch(u,WINDOW);

F = tf([0.505 1.01 0.505],[1  0.7478 0.2722],1,'Variable','q^-1');
Fresp = abs(squeeze(freqresp(F,w)));
[~,I] = min(abs(Fresp - 10^(-3/20)));

figure();
plot(w,(Pr),w,(Pu),w,Fresp,w(I),Fresp(I),'.r','MarkerSize',15);
grid on;
ylabel("Power");
xlabel('Normalized frequency [0, \pi]');
xlim([0 pi]);
xticks([0 0.25*pi 0.5*pi 0.75*pi pi]);
xticklabels({'0','0.25\pi','0.5\pi','0.75\pi','\pi'});
legend("Power r","Power u","Filter F","-3dB point");
%%
figure();
plot(t,r,'r',t,u,'b');
grid on;
legend("r","u");

R = PowerSpectrum(r,"R",1,1,false,false);

U = PowerSpectrum(u,"U",1,1,false,false);

Y = PowerSpectrum(y,"Y",1,1,false,false);

G = etfe(iddata(y,u,1),[]);
figure();
bode(G);
grid on;
