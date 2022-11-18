%% System identification homework assignment - part 5: validation
clear all; close all; clc


Fs = 1;                     % Sample frequency of signal (assumed)
Ts = 1/Fs;
Tend = 3000;
nSimulations = 1; %100;
FreqRes = 1/16;
t = 0:Ts:Tend;
N = length(t);

% Generating PRBS input signal
if true % PRBS input
    %PRBS
    n = (0:N-1);
    r = zeros(1,N);         % Pre-allocate
    
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
end
%

% Running simulation: obtaining input-output data of model
Gains = zeros(nSimulations,1);     % Pre-allocate memory
for i = 1:nSimulations
    % Generate input and output data
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
end
%

% Now data is generated:
% r = input filter
% u = output filter
% y = output system

% Part 5.1: using the PRBS input and model data, generate different order
% models for G(q,theta).

% Orders to test:
% M1: nb = 2, nf = 3
% M2: nb = 2, nf = 4
% M3: nb = 3, nf = 3
% M4: nb = 3, nf = 4
% M5: nb = 5, nf = 5
% M6: nb = 8, nf = 8
% M7: nb = 9, nf = 10
% M8: nb = 11, nf = 12
OeModelVector =     [2, 3;
                     2, 4;
                     3, 3;
                     3, 4;
                     5, 5;
                     8, 8;
                     9, 10;
                     11, 12;];

% Creating iddata object
data = iddata(y,u,Ts);
%                 
                 
% Fit all OE models
for i = 1:length(OeModelVector)
    % Fit model to data
    nb = OeModelVector(i,1);
    nf = OeModelVector(i,2);
    nk = 0;
    OeSys = oe(data,[nb nf nk]);
    %
    
    % Determine residuals
    [E,R] = resid(data,OeSys);
    %
      
    % Plot auto correlation
    LgndStr = ['n_b=',num2str(nb),', n_f=',num2str(nf)];
    
    figure(33)
    % Change figure position and size
    set(gcf,'Position',[0 0 1200 900])
    
    % Plot
    subplot 211
    hold all
    plot(R(:,1,1),'DisplayName',LgndStr)
    grid on
    legend('-DynamicLegend')
    title('Auto correlation')
    
    % Plot cross correlation
    subplot 212
    hold all
    plot(R(:,1,2))
    grid on
    title('Cross correlation')
    xlabel('Samples delay [n]')
    
    % Save cost
    if (i == 1)
        OeMSE = [i, OeSys.Report.Fit.MSE];
    else
        OeMSE = [OeMSE; i,OeSys.Report.Fit.MSE];
    end
end

% Determine V(theta)
figure(34)
hold all
plot(OeMSE(:,1),OeMSE(:,2));
grid on
title('Cost over models')
xlabel('Model number [n]')
ylabel('Cost')




%% Determining model (old code)
% See matlab resid function:
% https://www.mathworks.com/help/ident/ref/resid.html

if (0)
    close all;
    
    % Creating iddata object
    data = iddata(y,u,Ts);
    %
    
    % Create Box-Jenkins model
    nb = 5;
    nc = 0;
    nd = 0;
    nf = 7;
    nk = 0;
    sys = bj(data,[nb nc nd nf nk])
    %
    
    % Compute residuals
    [E,R] = resid(data,sys)
    % Where:
    % E = error
    % R(:,:,1) = autocorrelation
    % R(:,:,2) = cross-correlation
    %
    
    if true
        figure(9)
        plot(E)
        grid on
        
        figure(10)
        subplot 211
        hold on
        plot(R(:,1,1))
        plot(R(:,2,1))
        hold off
        grid on
        subplot 212
        hold on
        plot(R(:,1,2))
        plot(R(:,2,2))
        hold off
        grid on
        
        figure(11)
        I = impulseest(E);
        showConfidence(impulseplot(I,20),3)
    end
    
    % Determine orders to try
    nbval = 1:10;
    nfval = nbval;
    
    % Pre-allocate
    R11avg = zeros(length(nbval),2);
    R12avg = zeros(length(nbval),2);
    R21avg = zeros(length(nbval),2);
    R22avg = zeros(length(nbval),2);
    FitQuality = zeros(length(nbval),2);
    %
    
    %
    for nb = nbval
        for nf = nfval
            
            %
            sys = bj(data,[nb nc nd nf nk]);
            %
            
            FitQuality(nb,nf) = sys.Report.Fit.FitPercent;
            
            % Determine residuals
            [E,R] = resid(data,sys);
            %
            
            % Determine average values of autocorrelation
            R11avg(nb,nf) = sum(abs(R(:,1,1).^2))/(length(R(:,1,1))-1);
            R12avg(nb,nf) = sum(abs(R(:,1,2).^2))/(length(R(:,1,2))-1);
            R21avg(nb,nf) = sum(abs(R(:,2,1).^2))/(length(R(:,2,1))-1);
            R22avg(nb,nf) = sum(abs(R(:,2,2).^2))/(length(R(:,2,2))-1);
            %
            
        end
    end
    
    % Plotting of results
    figure()
    subplot 221
    surf(R11avg)
    title('R(:,1,1) average')
    xlabel('nb')
    ylabel('nf')
    zlabel('Sum')
    %
    subplot 222
    surf(R12avg)
    title('R(:,1,2) average')
    xlabel('nb')
    ylabel('nf')
    zlabel('Sum')
    %
    subplot 223
    surf(R21avg)
    title('R(:,2,1) average')
    xlabel('nb')
    ylabel('nf')
    zlabel('Sum')
    %
    subplot 224
    surf(R22avg)
    title('R(:,2,2) average')
    xlabel('nb')
    ylabel('nf')
    zlabel('Sum')
    %
    
    %
    figure()
    surf(FitQuality)
    xlabel('nb')
    ylabel('nf')
    zlabel('Fit quality [%]')
    title('Quality of model fit over number of poles/zeros')
    title('Determining fit quality of models w. varying nb and nf')
    %
    
end


