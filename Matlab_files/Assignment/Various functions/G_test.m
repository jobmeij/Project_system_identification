clear all; close all hidden; clc;

% Manual G approximation

s = tf('s');

T1 = pi*0.05;
G1 = 1/(s+T1);

wn2 = pi*0.22;
b2 = 0.08;
G2 = (wn2^2)/(s^2 + 2*b2*wn2*s + wn2^2);

wn3 = pi*0.5;
b3 = 0.08;
G3 = (wn3^2)/(s^2 + 2*b3*wn3*s + wn3^2);

G = 3*G1*G2*G3;
Gz = c2d(G,1);

w_vect = linspace(0,pi,1000);

% System Input

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
CrossCorrelate(u,y,1,200);
Gest = etfe(iddata(y,u,1),200);

Fresp_G = squeeze(freqresp(G,w_vect));
Fresp_Gest = squeeze(freqresp(Gest,w_vect));

f1 = figure();
set(f1,'position',[900 550 450 400]);
subplot(211)
plot(w_vect,(abs(Fresp_G)));
xlim([0 pi]);
ylim([0 60]);
xlabel("Frequency rad/s");
ylabel("Magnitude");
grid on;
subplot(212)
plot(w_vect,angle(Fresp_G));
xlim([0 pi]);
xlabel("Frequency rad/s");
ylabel("Phase");
grid on;
sgtitle("Bode plot G_{manual}");

f2 = figure();
set(f2,'position',[1400 550 450 400]);
subplot(211)
semilogx(w_vect,20*log10(abs(Fresp_G)));
xlim([0 pi]);
xlabel("Frequency rad/s");
ylabel("Magnitude");
grid on;
subplot(212)
semilogx(w_vect,angle(Fresp_G));
xlim([0 pi]);
xlabel("Frequency rad/s");
ylabel("Phase");
grid on;
sgtitle("Bode plot G_{manual}");

f3 = figure();
set(f3,'position',[900 50 450 400]);
subplot(211)
plot(w_vect,(abs(Fresp_Gest)));
xlim([0 pi]);
ylim([0 60]);
xlabel("Frequency rad/s");
ylabel("Magnitude");
grid on;
subplot(212)
plot(w_vect,angle(Fresp_Gest));
xlim([0 pi]);
xlabel("Frequency rad/s");
ylabel("Phase");
grid on;
sgtitle("Bode plot G_{est}");

f4 = figure();
set(f4,'position',[1400 50 450 400]);
subplot(211)
semilogx(w_vect,20*log10(abs(Fresp_Gest)));
xlim([0 pi]);
xlabel("Frequency rad/s");
ylabel("Magnitude");
grid on;
subplot(212)
semilogx(w_vect,angle(Fresp_Gest));
xlim([0 pi]);
xlabel("Frequency rad/s");
ylabel("Phase");
grid on;
sgtitle("Bode plot G_{est}");
