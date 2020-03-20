%%% SI model for a fictitious parasite 
%%% FB 11/05/2011
%%% Several transmission rates possible, including mass-action, or a mass-action modified by a threshold, or a saturation of the infection rate 
%%% Logistic prey growth produce dampened oscillations in mass-action forms, but these can be amplified by stochasticity... included in a random predation/mortality event every year here

clc
clear all
close all


global r beta N_thresh gamma mu K
r = 8; % lemming max growth rate
K = 500; % lemming carrying capacity
mu = 2; % Death rate of infected individuals 
N_thresh = K/10; % Threshold density at which infection begins
gamma = 8; % 
beta = 0.8;

y0=[10 0.1 0];
tstart = 0.0;
tstop = 50.0;

tspan=[tstart tstop];                        % timespan for the numerical solution
options=odeset('Events',@eventsParasite,'Reltol',1e-3,'NonNegative',[1 2]);

y1 = y0(1);                                   
y2 = y0(2);                                   
y3 = y0(3);                                   % time

yout = y0.';tout = 0; teout = []; yeout = []; ieout = [];  % output of events. 
test = [];

% ----------------------------------------------------------------------- %
while tout(length(tout))<tstop
    [t,y,TE,YE,IE] = ode45(@SI_func,[tstart tstop],y0,options)      % ode solver 
    count = length(t);                                              % length of time (t)
    tout = [tout;t(1:count)];                                       % time-out, number of time-steps
    teout = [teout;TE]; 
    yeout = [yeout;YE]; 
    ieout = [ieout;IE];      
% teout=time when events occur, yeout=values of species abundance at events, 
% ieout = 1 for event occuring
    y1 = [y1; y(1:count,1)];                                         % Lemming
    y2 = [y2; y(1:count,2)];                                         % Stoat
    y3 = [y3; y(1:count,3)];                                         % Time
    u = 0.75*rand; 						     % Perturbation 
    y0 = [y(count,1)*u; y(count,2)*u; y(count,3)];                   % Initial conditions for each timestep
    testadd = [y(count,2);u;y0(2)];                                  % test for Events-value 
    test = [test testadd];
    tstart = t(count);                                              % length of time-steps (tstop)
end;

% Plotting
subplot(311)
plot(tout,y1,teout,yeout(:,1),'-ko','LineWidth',2)
ylabel('Susceptible')
subplot(312)
plot(tout,y2,teout,yeout(:,2),'-ro','LineWidth',2)
ylabel('Infected')
subplot(313)
plot(tout,y1+y2,teout,yeout(:,1)+yeout(:,2),'-o','LineWidth',2)
ylabel('Total')
xlabel('Time')

print(figure(1),'-dpng','-r300','ALovelyParasite')
