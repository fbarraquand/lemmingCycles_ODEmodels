%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Gilg et al. 2003 model, clean integration using Matlab ODE45 %%%%
%%%% 18/01/2011 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%% 19/07/2013 Try to output more properly predation rates for the various predators
%%%% 01/10/2013 Back-calculation of the predation rates fully implemented with the various predator densities as outputs.
%%%% 20/11/2013 Leslie-like numerical response implemented 

clc
clear all
close all

global t_return_owl_fox t_return_skua t_snowmelt t_owl_hatch t_stoat_born t_skua_born t_skua_leaves t_owl_leaves t_fox_leaves r_w r_s v c D d_l d_h b W_o D_o W_l D_l W_f D_f b_o Y_o P_l_const b_f Y_f b1_o Y1_o b1_l Y1_l b1_f Y1_f Delta  Nprime s q

%------------------ Chronology of the year ---------------------%
% October 22 time zero, beginning of the year, tmod=0
% May 1 return of owls and foxes, tmod=(9+30+31+31+28+31+30+1)/365 = 191/365 = 0.523287671
% June 5 predation long-tailed skua, tmod=(9+30+31+31+28+31+30+31+5)/365 = 226/365 = 0.619178082
% June 15 SNOWMELT, fox cubs are born, change in lemming growth,  tmod=(9+30+31+31+28+31+30+31+15)/365 = 236/365 = 0.646575342
% June 25 Owls eggs hatching,  tmod=(9+30+31+31+28+31+30+31+25)/365 = 246/365 = 0.673972603
% July 1 stoat offspring are born,tmod=(9+30+31+31+28+31+30+31+30+1)/365 = 252/365 = 0.690410959
% July 14 skua chicks are born, tmod=(9+30+31+31+28+31+30+31+30+14)/365 = 265/365 = 0.726027397
% August 15 skua leaving, tmod=(9+30+31+31+28+31+30+31+30+31+15)/365 = 297/365 = 0.81369863
% September 30 owl leaving, tmod=(9+30+31+31+28+31+30+31+30+31+31+30)/365 = 297/365 = 0.939726027
% October 22 the fox leaves, tmod=(9+30+31+31+28+31+30+31+30+31+31+30+21)/365 = 365/365 = 1.0
%--------------------------------------------------------------%

t_return_owl_fox = 0.523287671;
t_return_skua = 0.619178082;
t_snowmelt = 0.646575342;
t_owl_hatch = 0.673972603;
t_stoat_born = 0.690410959; 
t_skua_born = 0.726027397;
t_skua_leaves = 0.81369863;
t_owl_leaves = 0.939726027; 
t_fox_leaves = 1.0;

%%%%%%%%% ECOLOGICAL PARAMETERS %%%%%%%%%%%%%%%%%%%%%

%---------- Lemmings Parameters --------------------%
r_w = 6; % up to 6 - blows up then...
r_s = 0.4; %or 0.4 % or 2

%---------- Stoat Parameters -----------------------%
v = 2.8; % 4 normally - should I use 2.8 as before?
c = 1000; % including surplus killing, per year
D = 0.1;
N_crit = D;
d_l = 0.1;
d_h = 4.0;
b = 25;
s = log(1+v);
q = 50.05; % Numbers of lemmings per stoat at carrying capacity -> does not work with large q -> exp growth

%---------- Generalists and nomadic ----------------%
%% Functional responses
W_o = 4.7*365;% so that it is in lemmings per year
D_o = 1.08;
W_l = 4.4*365;
D_l = 2.2;
W_f = 3.8*365;
D_f = 0.13;

%% Numerical responses (adults)
b_o = 0.00366;
Y_o = 2.86;
P_l_const = 0.02; %1.02 produce nice first order cycles with the skua alone
b_f = 0.0008; %0.0016 in JAs code
Y_f = 11.0;

%% Numerical responses (youngs)
b1_o = 0.011;
Y1_o = 4.0;
b1_l = 0.016;
Y1_l = 6.0;
b1_f = 0.0028;
Y1_f = 5.3;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%% NUMERICAL INTEGRATION %%%%%%%%%%%%%%%%%%%%

%---------- Initial conditions ----------------------%

tstart = 0;                                  
tstop = 50;                                
tspan=[tstart tstop];                        % timespan for the numerical solution
y0=[0.1 0.001 0];                            % initial values of state variables vector
%----------------------------------------------------%

%-------- OPTION HANDLE FOR EVENTS (Check Matlab documentation)-%            
options= odeset('Events',@events_gilg_03,'Reltol',1e-3,'NonNegative',[1 2 3]);  % setting the options for interruption of integration and change in predator density at regular intervals with stoat density multiplied by (1+v) 

y1=y0(1);                                   % for prey with mortality inflicted by generalist predators incorporated
y2=y0(2);                                   % for predator (stoat)
y3=y0(3);                                   % for time

yout = y0.'; tout=0; teout=[]; yeout=[]; ieout=[];       % output of events, teout=time when events occur, yeout=values of species abundance at events, ieout=
%generalists = [];
%seasonParams = [];
%PredOwl =[];
%PredSkua =[];
%PredFox =[];
%PredStoat =[];
Np = []; 

while (tout(length(tout))<tstop)			% Should work as well with the condition (tout(end)<tstop)		
    %tspan = linspace(tstart,tstop);		
    [t,y,TE,YE,IE] = ode45(@gilg_func_breakdown_Leslie,[tstart tstop], y0,options);       % ode solver for function = gilg_function 

%%%%%%%%%%%%%% NOTE FOR MYSELF, COMMENTED AND MODIFIED UP TO HERE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    length_time=length(t);                                               % length of time (t)
    tout = [tout;t(1:length_time)];                                      % time-out, number of time-steps
    yout = [yout;t(1:length_time)];                                      % output of y follows the time-steps
    teout=[teout;TE]; yeout=[yeout;YE]; ieout=[ieout;IE];
    y1=[y1;y(1:length_time,1)];
    y2=[y2;y(1:length_time,2)];
    y3=[y3;y(1:length_time,3)];

    Np=[Np;Nprime]; % iterate the density at snowmelt
	% Prepare for next step
    y0=[y(length_time,1);y(length_time,2)*(1);y(length_time,3)];          %%% Not multiplying by 1+v anymore
    test=[y(length_time,2);(1+v);y0(2)];                                   % test for showing value befor stoch value, the stoch value and value after adding stoch value
    tstart=t(length_time);                                               % length of time-steps (tstop)
	
end;
%%%%%%%%%%%%%% END OF SIMULATION %%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%% RECOMPUTE SPECIES-SPECIFIC PREDATION RATES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% NB. These will be obtained as vector - be careful....
%--- Time variable and density at snowmelt
tmod=mod(y3,1.0); % time vector

r=tmod; %initialize r vector for growth rates (giving also winter / summer)
r(tmod<t_snowmelt) = r_w;
r(tmod>t_snowmelt) = r_s;
y1p=tmod;  % density at snowmelt vector - must give the same density for the whole year. 

%% Nprime=lemming density at snowmelt, replaced by N before snowmelt
km = length(tmod)
yp1(1)=y1(1); % init loop
for (ki=1:(km-1))
	if (tmod(ki)<t_snowmelt)
	yp1(ki+1)=y1(ki+1);
	else 
	yp1(ki+1)=yp1(ki);
	end;
end;
plot(y3,yp1,y3,y1,'LineWidth',2)
xlabel('Time')
ylabel('Lemming density variables')
legend('Nprime','N')
print(figure(1),'-dpdf','-r300','DensitiesLemmingsGilgModel')


%%%%% Old, approximate way of back - calculating the predation rates. Gives rough answers, OK for skuas only ---------%
%------------------------------------ Pred Adult skua ----------------------------------------------------------------%
%index_skua = (mod(y3,1.0)>t_return_skua).*(mod(y3,1.0)<t_skua_leaves);
%PredSkua = P_l_const*W_l*index_skua.*(y1.^4)./(y1.^4+D_l^4); %% Adult predation
%------------------------------------ Pred Adult owl -----------------------------------------------------------------%
%index_owl = (mod(y3,1.0)>t_return_owl_fox).*(mod(y3,1.0)<t_owl_leaves);
%PredOwlAnc = b_o*W_l*index_owl.*(y1.^2)./(y1.^2+D_o^2); %% Adult predation, conservative estimates


%%%%%---------------- Proper computing of the total predation rates ----------------------------------------------%%%%%

%%% 1 - Numerical responses of generalists predators %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Use an inefficient yet correct "for loop" - to be sure of what one is computing (especially with the number of youngs)

AO=tmod;% initialization of adult and juvenile variables, AO = adult owls
AF=tmod;% will be replaced by smthg else anyway. Adult foxes, and so forth
AS=tmod;
JO=tmod;
JS=tmod;
JF=tmod;

for (ki=1:km), %%% Loop over the time index of the integrated variables y1,y2,y3
	%%% Adult Owls,  Nprime (=yp1(ki)) is replaced by N before snowmelt (see function)

	if (yp1(ki)>2.0)&&((tmod(ki)>t_return_owl_fox)&&(tmod(ki)<t_owl_leaves)), 
		AO(ki) = (b_o*(yp1(ki)-2.0))/(Y_o+yp1(ki)-4.0); 
	else 
		AO(ki) = 0.0; 
	end;

	%%% Long tailed skua (JA put before a condition (Nprime>2.0) ) 
	if ((tmod(ki)>t_return_skua)&&(tmod(ki)<t_skua_leaves)), 
		AS(ki) = P_l_const; 
	else 
		AS(ki) = 0.0;
	end; 
	%%% Foxes, Nprime replaced by N before snowmelt 
	if (yp1(ki)>2.0)&&((tmod(ki)>t_return_owl_fox)&&(tmod(ki)<t_fox_leaves)), 
		AF(ki) = (b_f*Nprime^2)/(Y_f^2+Nprime^2);  
	else 
		AF(ki) = 0.0;
	end;

%------------------------------------------ Young numbers -------------------------------------------------------------%
if (AO(ki)>0.0000001)&&((tmod(ki)>t_owl_hatch)&&(tmod(ki)<t_owl_leaves)), age = 365*(tmod(ki)-t_owl_hatch); JO(ki) =  ((b1_o*(yp1(ki)-2.0))/(Y1_o+yp1(ki)-4.0)) * 1.0/(1.0+exp(-0.36*(age-9.0))); else JO(ki) = 0.0;end;%%% Young Owls 
if (AS(ki)>0.0000001)&&((tmod(ki)>t_skua_born)&&(tmod(ki)<t_skua_leaves)), age = 365*(tmod(ki)-t_skua_born); JS(ki) =  (b1_l*(yp1(ki)^2)/(yp1(ki)^2+Y1_l^2)) * 1.0/(1.0+exp(-0.464*(age-4.55))); else JS(ki) = 0.0;end;%%% Young Skuas
if (AF(ki)>0.0000001)&&((tmod(ki)>t_snowmelt)&&(tmod(ki)<t_fox_leaves)), age = 365*(tmod(ki)-t_snowmelt); JF(ki) =  (b1_f*(yp1(ki)^2)/(yp1(ki)^2+Y1_f^2)) * 1.0/(1.0+exp(-0.36*(age-9.0))); else JF(ki) = 0.0;end;%%% Young Foxes (actualisation of lemming density at snowmelt)

end;

%%% Multiplication by functional responses %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%------------------------------------ Predation Generalists -----------------------------------------------------------%
PredOwl =(AO+JO).*W_o.*(y1.^2)./(y1.^2+D_o^2);
PredSkua = (AS+JS).*W_l.*(y1.^4)./(y1.^4+D_l^4);%%%% 
PredFox = (AF+JF).*W_f.*(y1.^2)./(y1.^2+D_f^2);

% Sum of variables above
PredGen = PredOwl + PredSkua + PredFox; %% Adult+Young predation

%------------------------------------- Predation Stoat ----------------------------------------------------------------%
PredStoat    = c.*y2.*(y1.^2)./(y1.^2+D^2); 


%%%%%%%%%%%%%% PLOTTING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% (to redo)
figure,
subplot(2,1,1)  
[AX,H1,H2] = plotyy(tout,y1,tout,y2,'semilogy');
ylim=[0 tstop];
ylim=[0 1]; %%% isnt that y2lim?
%title({r_w, r_s});
set(get(AX(1),'Ylabel'),'String','Lemming Density','Color','k')
set(get(AX(2),'Ylabel'),'String','Stoat Density','Color','r')
set(AX,'Xlim',[0 tstop])
set(get(AX(1),'Xlabel'),'String','Time [Years]','Color','k')
set(AX(1),'YColor','k','Ylim',[0 100])
set(AX(2),'YColor','r','Ylim',[0 1])
set(H1,'LineStyle','-','Color','k','LineWidth',2)
set(H2,'LineStyle','-.','Color','r','LineWidth',2)

subplot(2,1,2)
semilogy(tout,PredSkua,'-b',tout,PredStoat,'-r',tout,PredOwl,'-c',tout,PredFox,'-m','LineWidth',2);%,tout,PredSkua+PredStoat+PredOwl+PredFox,'-k'
%semilogy(tout,PredSkua,'-b',tout,PredStoat,'-r',tout,PredSkua+PredStoat,'-k','LineWidth',2);
ylabel('Total amount predated')
title({r_w, r_s});
legend('Skua','Stoat','Owl','Fox','Location','NorthWest')
print(figure(2),'-dpdf','-r300','TimeSeriesLogScale')

figure,
subplot(211)
[AX,H1,H2] = plotyy(tout,y1,tout,y2);
ylim=[0 tstop];
ylim=[0 max(y2)]; %%% isnt that y2lim?
%title({r_w, r_s});
set(get(AX(1),'Ylabel'),'String','Lemming Density','Color','k')
set(get(AX(2),'Ylabel'),'String','Stoat Density','Color','r')
set(AX(1),'YColor','k','Ylim',[0 max(y1)])
set(AX(2),'YColor','r','Ylim',[0 max(y2)])
set(H1,'LineStyle','-','Color','k','LineWidth',2)
set(H2,'LineStyle','-.','Color','r','LineWidth',2)
subplot(212)
MStoat =  PredStoat./y1;
MSkua = PredSkua./y1;
MOwl = PredOwl./y1; 
MFox = PredFox./y1,
Mtot=MOwl+MSkua+MStoat+MFox;
plot(tout,MStoat./Mtot,'r',tout,MSkua./Mtot,'b',tout,MFox./Mtot,'m',tout,MOwl./Mtot,'c','LineWidth',2)
legend('Stoat','Skua','Fox','Owl','Location','NorthWest')
axis([0 tstop 0 1])
ylabel('Fraction of lemming mortality')
xlabel('Time')
print(figure(3),'-dpdf','-r300','TimeSeriesPlusMortality')

figure,
plot(y1(1:end),y2(1:end))
title('Phase-Space plot')
xlabel('Lemming density')
ylabel('Stoat density')
print(figure(4),'-dpdf','-r300','PhaseSpace')

figure,
plot(y3,y2,'r',y3,AO+JO,'c',y3,AS+JS,'b',y3,AF+JF,'m','LineWidth',2)
ylabel('Predator Abundance (Adults + Juveniles)')
xlabel('Time (years)')
legend('Stoat','Owl','Skua','Fox','Location','NorthWest')
print(figure(5),'-dpdf','-r300','PredatorAbundance')

figure,
subplot(211)
plot(y1,c*(y1.^2)./(y1.^2+D^2),'o-');
xlabel('Prey density')
ylabel('Stoat Intake')
subplot(212)
plot(y1,W_l*(y1.^4)./(y1.^4+D_l^4),'o-');
xlabel('Prey density')
ylabel('Skua Intake')
print(figure(6),'-dpdf','-r300','FunctionalResponsesWithUse')

