%% Delay test
clear all; close all; clc

N = 200;
t = 0:N;
r = [zeros(10,1); ones(length(t)-10,1)*2];

ut = zeros(size(r));
yt = zeros(size(r));

A = 500;
for i = 1:A
    [u,y] = assignment_sys_25(r);
    ut = ut+u;
    yt = yt+y;
end

u = ut/A;
y = yt/A;

CrossCorrelate(y,u,"y u");
AutoCorrelate(y,"y");
AutoCorrelate(u,"u");

figure();
plot(t,r,'g.-',t,u,'b.-',t,y,'r.-');
grid on;

legend("Reference","Input","Output");