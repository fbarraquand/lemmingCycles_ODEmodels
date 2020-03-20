%%% Model 5 of Turchin and Batzli Ecology 2001
clc
clear all
close all

%%% TO DO
% - Try the seasonal version of that model in Matlab with the trick to stop the integration at time tau and replace U
% - Do the same thing with a logistic growth (put some option in the function?)
% - Then extend the seasonal analysis to the variable territory size model + analyze this one with regrowth as well
% Q: With these modifications, are we more likely to observe cycles between 3 and 5 years?

% NB Not sure the regrowth model can always apply in agricultural lands because of mowing... but there is external disturbance anyway in that case (i.e. vegetation appear and disappear, can be modelled mostly by changes in the vegetation carrying capacity)...

global A B G K U_S U_W tau R

% Parameters 
U_S = 10000;
U_W = 0.0;
tau = 0.5;
A = 15;
B = 70;
K = 2000;

G = 0.6*A;
r_max = 6;%6
R = r_max/(A*K/(K+B)-G);

% -------------------------------------- %

y0=[K 100 0.0];                            % initial values of state variables vector: plant - prey - predator - time 
tstart = 0;
tstop = 10;%30

%{ 
% Matlab integration
tspan=[tstart tstop];                        % timespan for the numerical solution
options= odeset('Reltol',1e-3,'NonNegative',[1 2]);
[tout,yout] = ode45(@tri_food_chain , tspan , y0,options);       % ode solver 
% Plotting
subplot(211)
plot(tout,yout(:,1))
subplot(212)
plot(tout,yout(:,2))
%}

tvec = linspace(tstart, tstop, 500);
y = lsode ("func_mod5_oct", [K;100], tvec);

subplot(211)
plot(tvec,y(:,1),tvec(1:50:end),y(1:50:end,1),'o')
subplot(212)
plot(tvec,y(:,2),tvec(1:50:end),y(1:50:end,2),'o')
