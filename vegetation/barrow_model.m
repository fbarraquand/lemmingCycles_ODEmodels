%%% Barrow model, FB 09/10/2013
%%% Model VIII of Turchin and Batzli Ecology 2001
clc
clear all
close all

global A B G_S G_W K_V K_M U_S U_W u_S u_W tau r_max R alpha

% Parameters 
U_S = 10000;
U_W = 0.0;
u_S = 12;
u_W = 0.0;
tau = 5/6; %before tau it is winter, after tau summer. tau = tau_critical. At t=0, we have V=V*0.10 (reduction in plant biomass)
A = 15;
B = 70;
K_V = 1000;
K_M = 2000;
G_S = 0.44*A;
G_W = 0.63*A;
R = 10.7/A; 
r_max = 6;%or 6,4? % irrelevant now
alpha = 0.1; % 0.1 to 1 . alpha = 0.5 usually

%%% Each year, we have V(t)=0.10*V(t-) at tau = 0 (summer-winter transition)

% --------------------------------------------------------------------------------------------------------------------- %

y0=[K_V/20 K_M/20 0.1 0];                            % initial values of state variables vector: plant - moss - lemming - time 
tstart = 0;
tstop = 50;
tspan=[tstart tstop];                        % timespan for the numerical solution
options= odeset('Events',@events_barrow,'Reltol',1e-6,'NonNegative',[1 2 3]);

y1=y0(1);                                   % 
y2=y0(2);                                   % 
y3=y0(3);
y4=y0(4);

yout = y0.'; tout=0; teout=[]; yeout=[]; ieout=[];       % output of events, teout=time when events occur, yeout=values of 

while (tout(length(tout))<tstop)			% Should work as well with the condition (tout(end)<tstop)		
    %tspan = linspace(tstart,tstop);		
    [t,y,TE,YE,IE] = ode45(@func_barrow,[tstart tstop], y0,options);       % ode solver for function = gilg_function 

    length_time=length(t);                                               % length of time (t)
    tout = [tout;t(1:length_time)];                                      % time-out, number of time-steps
    yout = [yout;t(1:length_time)];                                      % output of y follows the time-steps
    teout=[teout;TE]; yeout=[yeout;YE]; ieout=[ieout;IE];
    y1=[y1;y(1:length_time,1)];
    y2=[y2;y(1:length_time,2)];
    y3=[y3;y(1:length_time,3)];
    y0=[y(length_time,1)*0.10;y(length_time,2);y(length_time,3);y(length_time,4)];          % values of plant, moss, lemm and time at the end (tstop), i.e. initial value for each new step.
    tstart=t(length_time);                                               % length of time-steps (tstop)
	
end;
%%%%%%%%%%%%%% END OF SIMULATION %%%%%%%%%%%%%%%%%%%%%%


%%% Reconstruction of functional responses over time %%%%%%%%%%%%%%%%
FR_plant=A*y1./(y1+alpha*y2+B);
FR_moss=alpha*A*y2./(y1+alpha*y2+B);

percent_moss_in_diet=FR_moss./(FR_moss+FR_plant)
%%%%%%%%%%%%%%--------------------%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%% PLotting %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subplot(411)
plot(tout,y1)
ylabel('Vascular (kg/ha)')

subplot(412)
plot(tout,y2)
ylabel('Moss (kg/ha)')

subplot(413)
plot(tout,y3)
ylabel('Lemmings (inds/ha)')
xlabel('Time')

subplot(414)
plot(tout,percent_moss_in_diet)
ylabel('Fraction moss in diet')
xlabel('Time')

print(figure(1),'-dpdf','-r300','BarrowModel_alpha=01.pdf')

figure,
subplot(311)
ylabel('Plant')
semilogy(tout,y1)

subplot(312)
semilogy(tout,y2)
ylabel('Moss')

subplot(313)
semilogy(tout,y3)
ylabel('Lemmings')
xlabel('Time')

