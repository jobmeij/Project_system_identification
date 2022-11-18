clear all; close all hidden; clc;

N = 2048;
n = (0:N-1);
Nf = 128;
Ts = 1/(N);
Fs = 1/Ts;
t = n*Ts;
ws = (2*pi)/N;
Fmax = N;
Fsteps = round(Fmax/Nf);

r = zeros(1,N);
type = 2;
switch type
    case 1 %MULTI SIN
        for i = 1:Nf
            A = 0.06;
            F = Fsteps*(i);
            %p = ((2*pi)/(Nf+1)) * (i-1);
            p = rand(1,1)*2*pi;
            r = r + A*cos(2*pi*(F/Fs) * n + p);
        end
    case 2 %PRBS
        r(1) = (randi([0 1])-0.5)*4;
        p = 0.5;
        minclock = 1;
        clock = 0;
        for i = 2:N
            if rand() >= p && mod(clock,minclock) == 0
                r(i) = -r(i-1);
                clock = 0;
            else
                r(i) = r(i-1);
            end
            clock = clock + 1;
        end 
end

[u,y] = assignment_sys_25(r);

Gest = etfe(iddata(y,u,1),[]);
figure();
bode(Gest);
grid on;

WINDOW = hann(round(N/30));
[Pr,~] = pwelch(r,WINDOW);
[Pu,~] = pwelch(u,WINDOW);
[Py,w] = pwelch(y,WINDOW);

[U,~] = fft_hannwindow(u,4);
[Y,w2] = fft_hannwindow(y,4);

G = Y./U;
G2 = squeeze(freqresp(Gest,w2));

figure();
subplot(211);
plot(w2,20*log10(abs(G2)));
grid on;
ylabel("Magnitude [dB]");
xlabel('Normalized frequency [0, \pi]');
xlim([0 pi]);
xticks([0 0.25*pi 0.5*pi 0.75*pi pi]);
xticklabels({'0','0.25\pi','0.5\pi','0.75\pi','\pi'});
legend("G");
subplot(212);
plot(w2,(angle(G2)));
grid on;
ylabel("Phase [rad]");
xlabel('Normalized frequency [0, \pi]');
xlim([0 pi]);
xticks([0 0.25*pi 0.5*pi 0.75*pi pi]);
xticklabels({'0','0.25\pi','0.5\pi','0.75\pi','\pi'});

figure();
plot(w,20*log10(Py),w,20*log10(Pu));
grid on;
ylabel("Magnitude [dB]");
xlabel('Normalized frequency [0, \pi]');
xlim([0 pi]);
xticks([0 0.25*pi 0.5*pi 0.75*pi pi]);
xticklabels({'0','0.25\pi','0.5\pi','0.75\pi','\pi'});
legend("Py","Pu");

figure();
plot(w,20*log10(abs(Py./Pu)));
grid on;
ylabel("Magnitude [dB]");
xlabel('Normalized frequency [0, \pi]');
xlim([0 pi]);
xticks([0 0.25*pi 0.5*pi 0.75*pi pi]);
xticklabels({'0','0.25\pi','0.5\pi','0.75\pi','\pi'});