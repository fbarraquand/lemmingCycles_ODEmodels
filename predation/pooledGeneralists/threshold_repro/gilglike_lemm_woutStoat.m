%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% FB 10/11/2013 - Trying to get the essence of the Gilg model vs the Turchin-hanski and its developments.  
%%% Why is the Gilg model different from the modified version of TH97 with seasonal generalist predation?
%%% One hypothesis: It is because of the different numerical response. Looks like this explanation is the more likely...
%%% Testing this here by introducing a Gilg-like numerical response of the stoat into the modified TH97 model with seasonal generalist predation (the version I worked on following discussion with Rachel Taylor, and who to my surprise did not show generalist-generated declines despite the seasonality of generalist predation). 
%%%%%%%%%% Notes on model behaviour ------------------------------------ %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% It looks like generalist predation induces some declines. But more checks might be needed. 
% I came to wonder during debugging: in other cases (typically, nonseasonal gen predation or with the modified TH97) there could be very massive generalist predation >> specialist predation and yet specialist predation could be the driver of the declines if gen predation does not change too much... In that case, due to the timing, I find this doubtful - and I think we mimic the Gilg model. But why we do not find this in the seasonal generalist predation version of the TH97 is still to be decided...
%%%%%%%%% --------------------------------------------------------------- %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc 
clear all
close all

%%%%%% Parameters
global K G H C D t_repro d_high d_low r_s r_w

r_w = 1; % max value of the growth rate in winter
r_s = 0.8; % min summer value - changes a lot the dynamics to have both summer and winter GR, when it comes to the effect of predation rates 
K = 100000; % much more than at Trail island, sounds negligible, but in fact is not. 
G = 80; % 80 max value (30-40, low but plausible)
%%% Note: maybe output this from the empirical data
H = 3; % 2.75 or 5.0, half saturation generalists, maybe between 3 and 1 but this includes the numerical response as well...
C = 1000; % max consumption rate, per year, maybe lower
D = 0.1; % Maybe very low
t_repro = 0.4;% time of stoat repro - centered with that scheduling (winter is between 0.75 and 0.33)
v=2.8; % how many kids per stoat
d_high = 4.0;
d_low = 0.1;
% But maybe more, this seems to enlarge the range of possible 

y0=[0.1 0.00 0];
tstart = 0.0;
tstop = 30.0;

tspan=[tstart tstop];                        % timespan for the numerical solution
options= odeset('Events',@events_gilglike_lemm,'Reltol',1e-6','NonNegative',[1 2]);

%%% A few variables

y1=y0(1);                                   %  Initializing vectors
y2=y0(2);                                   % 
y3=y0(3);

yout = y0.'; tout=0; teout=[]; yeout=[]; ieout=[];       % output of events, teout=time when events occur, yeout=values of 

while (tout(length(tout))<tstop)			% Should work as well with the condition (tout(end)<tstop)		
    %tspan = linspace(tstart,tstop);		
    [t,y,TE,YE,IE] = ode45(@gilglike_lemm_func,[tstart tstop], y0,options);       % ode solver for function = gilg_function 

    length_time=length(t);                                               % length of time (t)
    tout = [tout;t(1:length_time)];                                      % time-out, number of time-steps
    %yout = [yout;t(1:length_time)];                                      % output of y follows the time-steps (why this btw?)
    teout=[teout;TE]; yeout=[yeout;YE]; ieout=[ieout;IE];
    y1=[y1;y(1:length_time,1)];
    y2=[y2;y(1:length_time,2)];
    y3=[y3;y(1:length_time,3)];	% time again

    y0=[y(length_time,1);y(length_time,2)*(1+v);y(length_time,3);];          % values of lemm, stoat, and time at the end (tstop), i.e. initial value for each new step.
    tstart=t(length_time);                                               % length of time-steps (tstop)
	
end;
%%%%%%%%%%%%%% END OF SIMULATION %%%%%%%%%%%%%%%%%%%%%%%% 

%%%%%%%%%%% Computation of %mortality %%%%%%%%%%%%%%%%%%%

%%% When is it winter?
nk=length(y3);
Winter=zeros(1,nk);
for (k=1:nk),
tmod=mod(y3(k),1.0);
if ((tmod>0.33)&&(tmod<0.74)), Winter(k) = 0; else Winter(k) = 1; end; % Step function formulation 
end;

%%% Sould reconstruct the winter variable
MReg =  ((Winter*r_w+(1-Winter)*r_s)').*y1/K; %%% Pop regulation -- see below for the use of ' for transposition
%%%%%%%%%%%%%% Previously generated an error %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% MGen =  (1-Winter).*G*(y1.^3)./(y1.^4+H^4);% g*(1-Winter)*y(1)^2/(y(1)^2+h^2)% in TH97 - ERROR IT SEEMS
%%% Why is this rate max 15 and then when I multiply it by a number either 0 or 1 it becomes 10^5?????
%%%%%%%%%%%%%% Bug solved but keep it in mind %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nsfr = G*(y1.^3)./(y1.^4+H^4); %%% Nonseasonal generalist FR. 
%%% Let's separate the computation and * this with the summer indicator
MGen =  ((1-Winter)').*nsfr; % This bug was because of a transpose!!!
%%% We multiplied 1 x n by n x 1 - Yet we did not get a unique value!? 
%%% WOW!!!!!! Impossible to predict?! That seems like a serious bug...

MSpe = C*y2.*y1./(y1.^2+D^2) ; %%% A bit different than from TH97 (here we had y1^2, above in the FR)

MTot=MReg+MSpe+MGen;
frac_gen = MGen./MTot;
frac_spe = MSpe./MTot; 
frac_reg = MReg./MTot;

%%%%%%%%%%%%%%%%% Plotting %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subplot(311)
plot(tout,y1)
ylabel('Lemming')
subplot(312)
plot(tout,y2)
ylabel('Stoat') 
subplot(313)
plot(tout,frac_reg,'k',tout,frac_gen,'b',tout,frac_spe,'r')% ,'LineWidth',2
%axis([0 tstop 0 1])
ylabel('Variable mortality (fraction)')
xlabel('Time')

figure,
subplot(2,1,1)  
[AX,H1,H2] = plotyy(tout,y1,tout,y2,'semilogy');
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

DeathReg = y1.*MReg;
PredSpe = y1.*MSpe;
PredGen = y1.*MGen;
subplot(2,1,2)
semilogy(tout,DeathReg,'-k',tout,PredSpe,'-r',tout,PredGen,'-b','LineWidth',2);
ylabel('Amount predated')







