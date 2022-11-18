function CrossCorrelate(U,Y,Tsin,Tin)
N = length(U);
if exist('Tin','var')
    T = Tin;
else
    T = N;
end
if exist('Tsin','var')
    Ts = Tsin;
else
    Ts = 1;
end
tauvect = (0:T-1)*Ts;
Rvect = zeros(length(tauvect),1);
for i = 1:length(tauvect)
    tau = i;
    Ytau = Y(tau:end);
    Rvect(i) = (1/(N-tau))*sum(U(1:end-tau+1).*Ytau);
end

f = figure();
set(f,'position',[700 600 500 400]);
plot(tauvect,Rvect);
xlabel('Tau');
ylabel('Correlation');
title('Cross Correlation');
grid on;
end


