%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Gilg et al. 2003 model, clean integration using Matlab ODE45 %%%%
%%%% 18/01/2011 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%% 19/07/2013 Try to output more properly predation rates for the various predators
%%%% 01/10/2013 Back-calculation of the predation rates fully implemented with the various predator densities as outputs.
%%%% 14/01/2015 Computation of extremas, to output how much is eaten by generalists. 

clc
clear all
close all

global t_return_owl_fox t_return_skua t_snowmelt t_owl_hatch t_stoat_born t_skua_born t_skua_leaves t_owl_leaves t_fox_leaves r_w r_s v c D d_l d_h b W_o D_o W_l D_l W_f D_f b_o Y_o P_l_const b_f Y_f b1_o Y1_o b1_l Y1_l b1_f Y1_f Delta  Nprime

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
r_w = 1; % 4 in Gilg et al. 2003, up to 6 - blows up then...
r_s = 0.4; %or 0.4 % or 2 or 0.8. In gilg et al. 2003 0, 0.4, 0.8

% Still cyclic for r_w=2, r_s=0.8
%---------- Stoat Parameters -----------------------%
v = 4.0; % 4 normally
c = 1000; % including surplus killing, per year
D = 0.08; % 0.08, 0.1, 0.12 . Lower D increases period. 
N_crit = D;
d_l = 0.1;% 0.0, 0.1, 0.2 % Lowering seems to increase period. 
d_h = 4;% 3.5, 4, 4.5 
b = 25;

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
tstop = 100;                                
tspan=[tstart tstop];                        % timespan for the numerical solution
%y0=[0.1 0.001 0];                            % initial values of state variables vector
%y0=[0.25 0.0025 0];                            % initial values of state variables vector
%y0=[10 0.03 0];                            % initial values of state variables vector

%%% First integration

y0=[4 0.1 0];                            % initial values of state variables vector

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
    [t,y,TE,YE,IE] = ode45(@gilg_func_breakdown,[tstart tstop], y0,options);       % ode solver for function = gilg_function 

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
    y0=[y(length_time,1);y(length_time,2)*(1+v);y(length_time,3)];          % values of prey, pred and time at the end (tstop), i.e. initial value for each new step.
    test=[y(length_time,2);(1+v);y0(2)];                                   % test for showing value befor stoch value, the stoch value and value after adding stoch value
    tstart=t(length_time);                                               % length of time-steps (tstop)
	
end;
%%%%%%%%%%%%%% END OF SIMULATION %%%%%%%%%%%%%%%%%%%%%%
y11=y1;y21=y2;y31=y3;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Second integration
%---------- Initial conditions ----------------------%

tstart = 0;                                  
tstop = 100;                                
tspan=[tstart tstop];                        

y0=[0.01 0.2 0];                            % initial values of state variables vector

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
    [t,y,TE,YE,IE] = ode45(@gilg_func_breakdown,[tstart tstop], y0,options);       % ode solver for function = gilg_function 

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
    y0=[y(length_time,1);y(length_time,2)*(1+v);y(length_time,3)];          % values of prey, pred and time at the end (tstop), i.e. initial value for each new step.
    test=[y(length_time,2);(1+v);y0(2)];                                   % test for showing value befor stoch value, the stoch value and value after adding stoch value
    tstart=t(length_time);                                               % length of time-steps (tstop)
	
end;
%%%%%%%%%%%%%% END OF SIMULATION %%%%%%%%%%%%%%%%%%%%%%
y12=y1;y22=y2;y32=y3;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%% PLOTTING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% (to redo)
figure,
subplot(2,1,1)  
[AX,H1,H2] = plotyy(y31,y11,y31,y21,'semilogy');
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
[AX,H1,H2] = plotyy(y32,y12,y32,y22,'semilogy');
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

print(figure(1),'-dpdf','-r300','TimeSeriesLogScale_N0')

%%% Phase-space for multiple initial conition
figure,
plot(log(y11(1:end)),log(y21(1:end)))
hold on
plot(log(y12(1:end)),log(y22(1:end)),'Color','k')
hold off
title('Phase-Space plot')
xlabel('log(Lemming density)')
ylabel('log(Stoat density)')
print(figure(2),'-dpdf','-r300','PhaseSpace_multipleN0')

