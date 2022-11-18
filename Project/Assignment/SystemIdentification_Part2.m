clear all; close all hidden; clc;

N = 1024;
n = (0:N-1);
Nf = 128;
Ts = 1/(N);
Fs = 1/Ts;
t = n*Ts;
ws = (2*pi)/N;
Fmax_V = [512*0.7];

legendtext = strings(length(Fmax_V),1);
Gresp = [];
for f = 1:length(Fmax_V)

    Fmax = Fmax_V(f);
    Fsteps = round(Fmax/Nf);
    r = zeros(1,N);
    type = 1;
    switch type
        case 1 %MULTI SIN
            for i = 1:Nf
                A = 0.06;
                F = Fsteps*(i);
                %p = ((2*pi)/(Nf+1)) * (i-1);
                p = rand(1,1)*2*pi;
                r = r + A*cos(2*pi*(F/Fs) * n + p);
            end
        case 2 %SINGLE SIN
            A = 2;
            F = 32;
            r = A*cos(2*pi*(F/Fs) * n);
        case 3 %WHITE NOISE
            r = (randi([0 1],N,1)-0.5)*2;
        case 4 %PRBS
            r(1) = (randi([0 1])-0.5)*4;
            p = 0.5;
            minclock = 2;
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
        case 5 %SIN SWEEP
            for i = 1:Nf
                length = N/Nf;
                time = t((i-1)*length+1:(i)*length);
                A = 1;
                f = fsteps*i;
                r((i-1)*length+1:(i)*length) = A*sin(2*pi*f * time );
            end
        case 6 %DC
            r = ones(1,N);
        case 7 %Ramp
            r = n/100;
        case 8 %Zeros

    end

    %r = r + 0.1;

    [u,y] = assignment_sys_25(r);

    % figure(1);
    % p = plot(t,r,'r',t,u,'b');
    % grid on;
    % legend("r","u");

    R = PowerSpectrum(r,"R",1,1,false,false);
    disp("Number of peaks in R: " + length(R.Pw(R.Pw>mean(R.Pw))));

    U = PowerSpectrum(u,"U",1,1,false,false);
    disp("Number of peaks in U: " + length(U.Pw(U.Pw>mean(U.Pw))));

    %F = tf([0.505 1.01 0.505],[1  0.7478 0.2722],1,'Variable','q^-1');
    %bode(F);

    %uF = lsim(F,r);
    %PowerSpectrum(u);

    %Y = PowerSpectrum(y,"Y",1,1,false,false);

    G = etfe(iddata(y,u,1),[]);
    [Gresp_ W] = (freqresp(G));
    Gresp_ = squeeze(Gresp_);
    Gresp = [Gresp Gresp_];
    
    legendtext(f) = string("G_{etfe}");
end

figure(5);

subplot(211)
semilogx(W,db(Gresp));
grid on
title('Bode diagram')
xlabel('Normalized frequency [0, \pi]')
ylabel('Magnitude [dB]')
xlim([0 pi])
xticks([0 0.1*pi 0.25*pi 0.5*pi 0.75*pi pi]);
xticklabels({'0','0.1\pi','0.25\pi','0.5\pi','0.75\pi','\pi'});
leg = legend(legendtext);
subplot(212)
semilogx(W,angle(Gresp))
grid on
xlabel('Normalized frequency [0, \pi]')
ylabel('Phase [rad]');
xlim([0 pi])
xticks([0 0.1*pi 0.25*pi 0.5*pi 0.75*pi pi]);
xticklabels({'0','0.1\pi','0.25\pi','0.5\pi','0.75\pi','\pi'});
ylim([-pi pi])
yticks([-pi -0.5*pi 0 0.5*pi pi]);
yticklabels({'-\pi','-0.5\pi','0','0.5\pi','\pi'});



% figure(6)
% set(gcf,'Position',[100 100 1200 600])
% subplot (2,2,[1,2])
% % Power spectrum r and u
% hold on
% plot(R.wvect,abs(R.A))
% plot(U.wvect,abs(U.A))
% hold off
% grid on
% title('Spectrum r and u')
% legend('r','u')
% xlabel('Normalized frequency [0, \pi]')
% ylabel('Magnitude [-]')
% xlim([0 pi])
% %ylim([0 0.07])
% xticks([0 0.25*pi 0.5*pi 0.75*pi pi]);
% xticklabels({'0','0.25\pi','0.5\pi','0.75\pi','\pi'});
% 
% subplot (2,2,[3,4])
% % Power spectrum y
% plot(Y.wvect,abs(Y.A))
% grid on
% title('Spectrum y')
% ylabel('Magnitude [-]')
% xlabel('Normalized frequency [0, \pi]');
% xlim([0 pi]);
% xticks([0 0.25*pi 0.5*pi 0.75*pi pi]);
% xticklabels({'0','0.25\pi','0.5\pi','0.75\pi','\pi'});

% Plotting for report assignment

% figure(7)
% set(gcf,'Position',[100 100 1200 600])
% subplot (2,2,1)
% % Power spectrum r
% plot(R.wvect,abs(R.A))
% grid on
% title('Spectrum r')
% xlabel('Normalized frequency [0, \pi]')
% ylabel('Magnitude [-]')
% xlim([0 pi])
% ylim([0 0.07])
% xticks([0 0.25*pi 0.5*pi 0.75*pi pi]);
% xticklabels({'0','0.25\pi','0.5\pi','0.75\pi','\pi'});
% 
% subplot (2,2,2)
% % Power spectrum u
% plot(U.wvect,abs(U.A))
% grid on
% title('Spectrum u')
% xlabel('Normalized frequency [0, \pi]')
% xlim([0 pi])
% ylim([0 0.07])
% xticks([0 0.25*pi 0.5*pi 0.75*pi pi]);
% xticklabels({'0','0.25\pi','0.5\pi','0.75\pi','\pi'});
% 
% subplot (2,2,[3,4])
% % Power spectrum y
% plot(Y.wvect,abs(Y.A))
% grid on
% title('Spectrum y')
% xlabel('Normalized frequency [0, \pi]')
% ylabel('Magnitude [-]')
% xlim([0 pi])
% xticks([0 0.25*pi 0.5*pi 0.75*pi pi]);
% xticklabels({'0','0.25\pi','0.5\pi','0.75\pi','\pi'});
