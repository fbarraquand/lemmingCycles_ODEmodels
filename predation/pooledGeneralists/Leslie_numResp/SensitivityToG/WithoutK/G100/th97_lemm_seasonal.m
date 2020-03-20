%%% Model of Turchin and Hanski (Am Nat 1997) adapted to lemmings so that we can model the Traill island population with it
%%% Barraquand & Henden 2011
%%% Seasonal version
%%% Main thing yet to do: it needs a longer winter period and a smaller summer period...

%%%% NB: Play on H to change period/amplitude without too unrealistic assumptions...

%%%% Add-on FB 04/06/2014. Output of percentage of mortality to see what is due to the generalist and specialists, with such reproduction schedule
%%%% Without K

clc 
clear all
close all

%%%%%% Parameters
global r r_min r_max K G H C D s Q

r_max = 6; % max value
r_min = 0.5; % max growth rate bad season 
%K = 500; % much more than lemming obs densities at Traill island, negligible
G = 100; % max value (30-40, low but plausible)
%%% Note: maybe output this from the empirical data
H = 2; %3.5; %2.75; % or 5.0, half saturation generalists, maybe between 3 and 1 but this includes the numerical response as well...
C = 1000; % max consumption rate, per year, maybe lower
D = 0.1; % Maybe very low
s = 1.75; % 1.5 here this is the max, average = half
Q = 100; % 30 from Turchin and Hanski 1997, should hold roughly for lemmings because stoats are larger
% But maybe more, this seems to enlarge the range of possible 

y0=[0.1 0.001 0];
tstart = 0.0;
tstop = 30.0;

tspan=[tstart tstop];                        % timespan for the numerical solution
options= odeset('Reltol',1e-3,'NonNegative',[1 2]);
[tout,yout] = ode45(@th_lemm_seasonal_func, tspan , y0,options);       % ode solver 

% Plotting
subplot(211)
plot(tout,yout(:,1))
ylabel('Lemming')
subplot(212)
plot(tout,yout(:,2))
ylabel('Stoat')
xlabel('Time')

figure,
subplot(2,1,1)  
[AX,H1,H2] = plotyy(tout,yout(:,1),tout,yout(:,2),'semilogy');
ylim=[0 tstop];
y2lim=[0 1]; %%% isnt that y2lim?
set(get(AX(1),'Ylabel'),'String','Lemming Density','Color','k')
set(get(AX(2),'Ylabel'),'String','Stoat Density','Color','r')
set(AX,'Xlim',[0 tstop])
set(get(AX(1),'Xlabel'),'String','Time [Years]','Color','k')
set(AX(1),'YColor','k','Ylim',[0 100])
set(AX(2),'YColor','r','Ylim',[0 1])
set(H1,'LineStyle','-','Color','k','LineWidth',2)
set(H2,'LineStyle','-.','Color','r','LineWidth',2)

%%% Output mortalities

% Threshold version of seasonality
Winter=2*ones(length(tout),1);
 tmod = mod(tout,1.0);
Winter((tmod>0.25)&(tmod<0.75))=0;
% Compute mortalities
MReg =  0*yout(:,1); %%% Perhaps that one was not computed well before
MGen =  G*yout(:,1).*(1-Winter)./(yout(:,1).^2+H^2);% G*(1-Winter)*y(1)^2/(y(1)^2+h^2)% NB: G is the max here, not the mean
MSpe = C*yout(:,2)./(yout(:,1)+D);
Mtot=MReg+MSpe+MGen;

% Plotting
subplot(311)
plot(tout,yout(:,1),'LineWidth',2)
ylabel('Lemming')
subplot(312)
plot(tout,yout(:,2),'LineWidth',2)
ylabel('Stoat')
subplot(313)
plot(tout,MReg./Mtot,'b',tout,MGen./Mtot,'k',tout,MSpe./Mtot,'r','LineWidth',2)
axis([tstart tstop 0 1])
ylabel('Variable mortality (fraction)')
xlabel('Time')
title('Seasonal r and seasonal G')

figure,
subplot(2,1,1)  
[AX,H1,H2] = plotyy(tout,yout(:,1),tout,yout(:,2),'semilogy');
ylim=[0 tstop];
y2lim=[0 1]; %%% isnt that y2lim?
set(get(AX(1),'Ylabel'),'String','Lemming Density','Color','k')
set(get(AX(2),'Ylabel'),'String','Stoat Density','Color','r')
set(AX,'Xlim',[0 tstop])
set(get(AX(1),'Xlabel'),'String','Time [Years]','Color','k')
set(AX(1),'YColor','k','Ylim',[0 100])
set(AX(2),'YColor','r','Ylim',[0 1])
set(H1,'LineStyle','-','Color','k','LineWidth',2)
set(H2,'LineStyle','-.','Color','r','LineWidth',2)

DeathReg = yout(:,1).*MReg;
PredSpe = yout(:,1).*MSpe;
PredGen = yout(:,1).*MGen;
subplot(2,1,2)
%semilogy(tout,DeathReg,'-k',tout,PredSpe,'-r',tout,PredGen,'-b','LineWidth',2); % when regulation is on
semilogy(tout,PredSpe,'-r',tout,PredGen,'-b','LineWidth',2); 
%legend('regulation','specialist','generalists')% regulation on
legend('specialist','generalists')
ylabel('Amount killed')
%print(figure(3),'-dpdf','-r300','TimeSeriesLogScale_LemmingLikeTH97')
print(figure(3),'-dpdf','-r300','TimeSeriesLogScale_LemmingLikeTH97_SinusoidalSeasonality')



