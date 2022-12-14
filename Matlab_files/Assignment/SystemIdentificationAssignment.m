%% System Identification assignment
clear all; close all; clc;


Fs = 1;                     % Sample frequency of signal (assumed)
Ts = 1/Fs;
nSimulations = 1 %100;
FreqRes = 1/16;
StartFreqOffset = 0;

tic 

Gains = zeros(nSimulations,1);     % Pre-allocate memory
for i = 1:nSimulations

% Load data
t = 1:1e3;                    % Length of data [samples] 
% r = sin((StartFreqOffset+0*i*FreqRes)*t);
% r = chirp(t,1/512,length(t),1,'logarithmic');
% r = randi([0 1],length(t),1);
r = 4*(randi([0 1],length(t),1)-0.5);
% r = t;
[u,y] = assignment_sys_25(r);   % Extracting data from encrypted file

% Plotting data
if (1)
    figure(1)
    hold on
    plot(t,r)
    plot(t,u)
    plot(t,y)
    hold off
    grid on
    title('Data')
    legend('r','u','y')
    xlabel('Sample [n]')
    ylabel('Amplitude [-]')
end

% FFT
fftr = fft(r,length(t));
fftr = (1/length(t))*abs(fftr);
% Determining power spectrum of r
nr = length(r);                         % Number of samples
fr = (0:nr-1)*(Fs/nr)';                 % Frequency range
powerr = abs(fftr).^2/nr;               % Power of the DFT

fftu = fft(u,length(t));
fftu = (1/length(t))*abs(fftu);
% Determining power spectrum of u
nu = length(u);                         % Number of samples
fu = (0:nu-1)*(Fs/nu)';                 % Frequency range
poweru = abs(fftu).^2/nu;               % Power of the DFT

ffty = fft(y,length(t));
ffty = (1/length(t))*abs(ffty);
% Determining power spectrum of y
ny = length(y);                         % Number of samples
fy = (0:ny-1)*(Fs/ny)';                 % Frequency range
powery = abs(ffty).^2/ny;               % Power of the DFT

if (1)
    figure(2)
    subplot 211
    semilogx(t,db(fftr))
    hold on
    semilogx(t,db(fftu))
    semilogx(t,db(ffty))
    hold off
    grid on
    ylabel('Amplitude [dB]')
    legend('r','u','y','Location','Best')
    title('Bodeplot input-output data')
    
    subplot 212
    semilogx(t,angle(fftr))
    hold on
    semilogx(t,angle(fftu))
    semilogx(t,angle(ffty))
    hold off
    grid on
    ylabel('Angle [deg]')
    xlabel('Frequency [Hz?? check Ts!]')
    
    figure(10)
    semilogx(fr,db(powerr))
    hold on
    semilogx(fu,db(poweru))
    semilogx(fy,db(powery))
    hold off
    grid on
    xlabel('Frequency [??]')
    ylabel('Amplitude [-]')
end

% Determine gain and phase from u to y
GainSys = (1/length(t))*(abs(ffty)./abs(fftu));
PhaseSys = (angle(ffty)./angle(fftu))*(180/pi);

% Determine gain of system from peaks of u and y
Gain = (10^(max(db(ffty))/20))/(10^(max(db(fftu))/20));

% Save gain and input frequency of u 
Gains(i) = db(Gain);
Frequency(i) = (StartFreqOffset+i*FreqRes);

end

toc

if (1)
    figure(3)
    plot(Frequency,Gains)
    grid on
    xlabel('Frequency [??]')
    ylabel('Gain [dB]')
    title('Frequency response')
    
    figure(4)
    subplot 211
    semilogx(1:length(GainSys),GainSys)
    grid on
    ylabel('Gain [dB]')
    title('Bode system')
    
    subplot 212
    semilogx(1:length(PhaseSys),PhaseSys)
    grid on
    ylabel('Phase [deg]')
    xlabel('Frequency [??]')
end

%% Fit PDF to data from y (noise)
pdFunc = fitdist(y,'Normal')
yMin = min(y);
yMax = max(y);
pdfRes = 0.1;
yVals = yMin:pdfRes:yMax;
pdfy = pdf(pdFunc,yVals);

figure()
plot(yVals,pdfy,'LineWidth',2)
grid on
xlabel('Value y')
ylabel('Distribution')
title('Distribution of values on y')

%% Tfest function - estimating G0(q) and plotting bode
data = iddata(y,u,Ts);
nz = 6;
np = 6;
G0_est = tfest(data,nz,np)        % NOTE: estimating with 6,6 gives best results

figure(15)
hold all
bode(G0_est)
grid on
title({'Estimation of G_0(q)','Tfest: Input u and output y'})
