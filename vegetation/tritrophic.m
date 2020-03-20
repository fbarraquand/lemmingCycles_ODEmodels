%%% Tritrophic model from Turchin and Batzli Ecology 2001
clc
clear all
close all

global A B C D G K U R s Q

% Parameters 
U = 10000;
A = 15;
B = 70;
K = 2000;

G = 0.6*A;
r_max = 6;
R = r_max/(A*K/(K+B)-G);

C = 600;
D = 6;
s = 1.25;
Q = 40;

% -------------------------------------- %

y0=[K 100 0.01 0];                            % initial values of state variables vector: plant - prey - predator - time 
tstart = 0;
tstop = 30;
tspan=[tstart tstop];                        % timespan for the numerical solution
options= odeset('Reltol',1e-3,'NonNegative',[1 2 3]);

[tout,yout] = ode45(@tri_food_chain , tspan , y0,options);       % ode solver for function = gilg_function 

subplot(311)
plot(tout,yout(:,1))

subplot(312)
plot(tout,yout(:,2))

subplot(313)
plot(tout,yout(:,3))

