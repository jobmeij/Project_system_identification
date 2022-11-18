%% System identification homework assignment - part 5: validation
clear all; close all; clc

N = 3000;
n = (0:N-1);
Ts = 1/(N);
t = n*Ts;

%PRBS2
Nc = 2;
R = 1;
white = randn(N/(Nc*R),1);
r = repmat(2*sign(white(ceil([1:(N/R)]/Nc))),[R,1]);

% G = K(B/F)
% H = C/D
nb_V = [2 2 3 3 5 8 9 11];
nc = 1;
nd = 2;
nf_V = [3 4 3 4 5 8 10 12];
nk = 1;

Nrepetitions = 5;
Nsimulations = 1;

for rep = 1:Nrepetitions
    
    Cost_V = zeros(length(nb_V),Nsimulations);
    Cost2_V = zeros(length(nb_V),Nsimulations);

    [u,y] = assignment_sys_25(r);
    data = iddata(y,u,Ts);

    for s = 1:Nsimulations

        [u2,y2] = assignment_sys_25(r);
        data2 = iddata(y2,u2,Ts);

        for i = 1:length(nb_V)
            nb = nb_V(i);
            nf = nb_V(i);

            %OeSys = oe(data,[nb nf nk]);
            BjSys = bj(data,[nb nc nd nf nk]);

            % Determine residuals same Dataset
            [E,R] = resid(data,BjSys);
            epsilon = E.OutputData;

            % Determine residuals different Dataset
            [E2,R2] = resid(data2,BjSys);
            epsilon2 = E2.OutputData;

            %Cost 1/N sum(epsilon^2) from 0 to N-1
            N = length(epsilon);
            Cost_V(i,s) = (1/N)*sum(epsilon.^2);  

            N2 = length(epsilon2);
            Cost2_V(i,s) = (1/N2)*sum(epsilon2.^2);
        end
    end

    Cost = min(Cost_V,[],2);
    Cost2 = min(Cost2_V,[],2);

    figure();
    %C = plot(Cost_V,'b.','MarkerSize',20);
    Cmin = plot(Cost,'b.-','MarkerSize',20);
    hold on;
    %C2 = plot(Cost2_V,'r.','MarkerSize',10);
    C2min = plot(Cost2,'r.-','MarkerSize',10);
    grid on;
    title('Cost over models');
    xlabel('Model');
    xticklabels("M_{"+(1:8)+"}");
    ylabel('Cost');
    legend([Cmin(1) C2min(1)],'Cost','Cost with new data set');
end


