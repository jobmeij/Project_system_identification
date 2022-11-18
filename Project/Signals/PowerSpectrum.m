function Output = PowerSpectrum(U,name,Tsin,NwindowWin,power,db)
if exist('Tsin','var')
    Ts = Tsin;
else
    Ts = 1;
end
if exist('NwindowWin','var')
    NwindowW = 2*ceil(NwindowWin/2)-1;
else
    NwindowW = 1;
end
if exist('name','var')
    text = "spectrum "+name;
else
    text = 'spectrum';
end
if exist('db','var')
    plotdb = db;
else
    plotdb = false;
end
if exist('power','var')
    plotpower = power;
else
    plotpower = true;
end

U = reshape(U,[1 length(U)]);
if mod(length(U),2) ~= 0
    U = [U 0];
end

N = length(U);
t = (0:N-1)*Ts;
ws = (2*pi)/(Ts);
l = 0:N/2;
wvect = (l/N)*ws;

Uw = zeros(length(wvect),1);
for i = 1:length(wvect)
    Uwi = 0;
    for d = -floor(NwindowW/2):floor(NwindowW/2)
        fi = max(min(i+d,length(wvect)),1);
        w = wvect(fi);
        Uwi = sum(U.*exp(-sqrt(-1)*w*t)) + Uwi;
    end
    Uw(i) = Uwi/NwindowW;
end
Uw = [Uw(1); Uw(2:end)*2]/N;
Pw = (1/(2*N^2))*(abs(Uw).^2);

f = figure();
set(f,'position',[1200 600 500 400]);
if plotpower
    if plotdb
        plot(wvect,20*log10(abs(Pw)));
        ylabel('Power [dB]');
    else
        plot(wvect,Pw);
        ylabel('Power');
    end
else 
    if plotdb
        plot(wvect,20*log10(abs(Uw)));
        ylabel('Magnitude [dB]');
    else
        plot(wvect,abs(Uw));
        ylabel('Magnitude');
    end
end

xlabel('Frequency');
title(text);
grid on;

Output.A = Uw;
Output.Pw = Pw;
Output.wvect = wvect;
end
