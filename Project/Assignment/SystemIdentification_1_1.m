clear all;
close all hidden;
clc;

F = tf([0.505 1.01 0.505],[1  0.7478 0.2722],1,'Variable','q^-1')

figure();
bode(F);
grid on;

w_vect = linspace(pi/100,pi,3000);
Fresp = squeeze(freqresp(F,w_vect));

figure();
subplot(211)
plot(w_vect,20*log10(abs(Fresp)));
xlabel("Frequency rad/s");
ylabel("Magnitude (dB)");
grid on;
subplot(212)
plot(w_vect,angle(Fresp));
xlabel("Frequency rad/s");
ylabel("Phase");
grid on;
sgtitle("Bode plot Butterworth filter");


