function Output = SystemBode(Sys,ws)

wvect = ws:ws:ws*(length(Sys));

Mag = 20*log10(abs(Sys));
Phase = angle(Sys);

f = figure();
set(f,'position',[1200 600 500 400]);
subplot(211);
semilogx(wvect,Mag);
ylabel('Magnitude [dB]');
grid on;
subplot(212);
semilogx(wvect,Phase);
xlabel('Frequency');
ylabel('Phase [rad]');
grid on;
sgtitle("Bode");

Output.Mag = Mag;
Output.Phase = Phase;
Output.wvect = wvect;
end
