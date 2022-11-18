%% System identification homework assignment - part 6: Theoretical verification
clear all; close all; clc

N = 1000;
n = (0:N-1);
Ts = 1/(N);
t = n*Ts;

%PRBS2
Nc = 2;
R = 1;
white = randn(N/(Nc*R),1);
r = repmat(2*sign(white(ceil([1:(N/R)]/Nc))),[R,1]);

nk = 1;
nb = 5;
nf = 5;
Ntheta = nb+nf;

Iterations = 100;

Model.F = zeros(Iterations, nf);
Model.B = zeros(Iterations, nb);
VariancesT = zeros(Iterations, nf+nb);

for s = 1:2

text = ["random initial guess","median initial guess"];
disp(text(s));
    
B0 = [zeros(1,nk) median(Model.B)];
F0 = [1 median(Model.F)];
M0 = idpoly([],B0,[],[],F0);

for i = 1:Iterations
    [u,y] = assignment_sys_25(r);
    Data = iddata(y,u,Ts);
    switch s
        case 1
            M = oe(Data,[nb nf nk],oeOptions('EnforceStability',true));
            Covariance = getcov(M);
            VariancesT(i,:) = diag(Covariance);
        case 2
            M = oe(Data,M0,oeOptions('EnforceStability',true));
    end
    Model.F(i,:) = M.F(2:end);
    Model.B(i,:) = M.B(nk+1:end);
end

VarT = median(VariancesT);

f = figure();
Ax = axes();
grid on;
%ylim([-10 10]);

v = ones(size(Model.F,1),1);
for i = 1:Ntheta
    if i <= size(Model.B,2)
        p = i;
        Bline = line(Ax,'XData',i*v,'YData',Model.B(:,p),'Marker','.',...
            'color','b','Markersize',10,'LineStyle','none');
    else
        p = i-size(Model.B,2);
        Fline = line(Ax,'XData',i*v,'YData',Model.F(:,p),'Marker','.',...
            'color','r','Markersize',10,'LineStyle','none');
    end
end

MeanB = mean(Model.B);
VarB = var(Model.B);
MeanF = mean(Model.F);
VarF = var(Model.F);
MeanC = [MeanB MeanF];
VarC = [VarB VarF];

hold on;
variancebar = errorbar(1:Ntheta,MeanC,VarC,'LineStyle','none');
legend([Fline Bline variancebar],"Parameters F","Parameters B","Variance");
xticklabels(["B_{"+string(0:size(Model.B,2)-1)+"}" "F_{"+string(1:size(Model.F,2))+"}"]);
xticks(1:Ntheta);
xlabel("Parameter");
ylabel("Value");
title("Parameters "+text(s));
saveas(f,fullfile([pwd '\Figures6'],"Parameters "+text(s)),'fig'); %Auto save figures
saveas(f,fullfile([pwd '\Figures6'],"Parameters "+text(s)),'png'); %Auto save figures

f = figure();
plot(VarC,'b.','MarkerSize',20);
hold on;
plot(VarT,'r.','MarkerSize',20);
grid on;
xticklabels(["B_{"+string(0:size(Model.B,2)-1)+"}" "F_{"+string(1:size(Model.F,2))+"}"]);
xticks(1:Ntheta);
xlabel("Parameter");
ylabel("Variance");
title("Variance "+text(s));
legend("Calculated variance","Theoretical variance");
saveas(f,fullfile([pwd '\Figures6'],"Variance "+text(s)),'fig'); %Auto save figures
saveas(f,fullfile([pwd '\Figures6'],"Variance "+text(s)),'png'); %Auto save figures

end

G = tf([zeros(1,nk) Model.B(1,:)],[1 Model.F(1,:)],1,'Variable','q^-1');
figure();
bode(G);
grid on;

