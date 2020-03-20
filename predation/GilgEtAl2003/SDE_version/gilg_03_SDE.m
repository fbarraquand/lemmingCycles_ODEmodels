%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Gilg et al. 2003 model, brute force integration %%%%
%%%% 18/01/2011 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Without the stoat this time.... 
%%%% 18/04/2015 Adding stochasticity (perhaps re-adding the stoat after too...)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc
clear all
close all

randn('seed',9)


%% Checks -> (1/dt)*log(N(2:11)./N(1:10)) seems to give the right number yet we don't get the right growth. 
%% Except 0.25*exp(4/2) yields the correct peak at 1.84... Would the other integrations be wrong?
%% The second stagnation is surprising...

%------------------ Chronology ---------------------%
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
%---------------------------------------------------%
t_zero = 0.0;
t_return_owl_fox = 0.523287671;
t_return_skua = 0.619178082;
t_snowmelt = 0.646575342;
t_owl_hatch = 0.673972603;
t_stoat_born = 0.690410959;
t_skua_born = 0.726027397;
t_skua_leaves = 0.81369863;
t_owl_leaves = 0.939726027;
t_fox_leaves = 1.0;
%%%%%%%%%%%%%%%%%%%% PARAMETERS %%%%%%%%%%%%%%%%%%%%%

%---------- Lemmings Parameters --------------------%
r_w = 4.0;
r_s = 0.8; %or 0.4 or 0.8

% carrying capacity in absence pred. when needed
%K = 100; %K=500 makes the pop. escape to the carrying capacity in absence of weasels. Happens stochastically with K=100 now too. 
K=50

%---------- Stoat Parameters -----------------------%
v = 4.0;
c = 1000; % including surplus killing
D = 0.1;
N_crit = D;
d_l = 0.1;
d_h = 4.0;
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
P_l_const = 0.02;
b_f = 0.0008;
Y_f = 11.0;

%% Numerical responses (youngs)
b1_o = 0.011;
Y1_o = 4.0;
b1_l = 0.016;
Y1_l = 6.0;
b1_f = 0.0028;
Y1_f = 5.3;

%---------------- Stochasticity ---------------------------------% 

sigmaE = 0.25; % Intrinsic noise SD (I started with sigmaE=0.1 or 0.5)
% Cycles appear well even without predator for sigmaE=0.1
%%% Would be useful to know exactly how to scale it... W(t=2)-W(t=1)~N(0,1) so  sigmaE*(W(t=2)-W(t=1))~N(0,sigmaE^2)
%%% The infinitesimal SD is in N*sigmaE, so to be multiplied by the amount of time to get the real SD


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initializing % 
r=0.0;

%--------------- Integration parameters -------------------------%
intervals = 50000, %365; %50 -> dt=0.02
dt=1.0/intervals;  %-> NB decrease timestep wrt RK4
T_max=30; %50.0;
N(1) = 0.25;%lemming initial density - works up to 1, blow up at 10
P(1) = 0.00;
%0.001; %stoat initial density

%%% Recording
Po(1)=0; %% Note it is different from P_o, the intermediate variable
Pf(1)=0;
Pl(1)=0;
Nb(1)=0;

%--------------- Loop on time -----------------------------------%
ki=1;
t=0.0;
snowmelt=1;
stoat_born=0;
N1=0.0;

while(t<T_max)
tmod=mod(t,1.0); % Variable between 0 and 1 indicating time of year
%%% case statement to simplify here ???
if (snowmelt==1)&&((tmod>0)&&(tmod<dt)), 
N1=0;% or 0?
stoat_born=0;
snowmelt=0;
end; % set lemming density at snowmelt to zero and Indicator variables for whether snowmelt or the stoat birth has happened this year already

if ((snowmelt==0)&&(tmod>t_snowmelt)),
N1=N(ki); 
snowmelt=1;
end

if (tmod<t_snowmelt), r=r_w; else r=r_s;end; % Which lemming growth rate is used

%%% Lemming growth
%PreyGrowth = r*N(ki); 
PreyGrowth =r*N(ki)*(1-N(ki)/K); %if logistic growth needed

%%% Numerical responses
%%% Adult numbers
if (tmod<t_snowmelt), Nbis = N(ki); else Nbis = N1;end; %% N1 replaced by N before snowmelt

if (Nbis>2.0)&&((tmod>=t_return_owl_fox)&&(tmod<t_owl_leaves)), 
	%P_o = (2*b_o*(Nbis-2.0))/(Y_o+Nbis-4.0); % why times 2 by the way?
	P_o = (b_o*(Nbis-2.0))/(Y_o+Nbis-4.0);
else P_o = 0.0;end;%%% Owls,  N1 replaced by N before snowmelt (introduced variable N2 and if condition)
%P_o=0.0;

if ((tmod>=t_return_skua)&&(tmod<t_skua_leaves)), P_l = P_l_const; else P_l = 0.0;end;%%% Long tailed skua
%
if (Nbis>2.0)&&((tmod>=t_return_owl_fox)&&(tmod<t_fox_leaves)), P_f = (b_f*Nbis^2)/(Y_f^2+Nbis^2);  else P_f = 0.0;end;
%%% Foxes, N1 replaced by N before snowmelt %%% Weird expression in paper - recheck
%P_f=0.0; % avoid foxes


%%% Young numbers
%P_yo=0.0;
if (P_o>0.0000001)&&((tmod>=t_owl_hatch)&&(tmod<t_owl_leaves)), age = 365*t-t_owl_hatch; P_yo =  ((b1_o*(Nbis-2.0))/(Y1_o+Nbis-4.0)) * 1.0/(1.0+exp(-0.36*(age-9.0))); else P_yo = 0.0;end;%%% Young Owls 
if (P_l>0.0000001)&&((tmod>=t_skua_born)&&(tmod<t_skua_leaves)), age = 365*t-t_skua_born; P_yl =  (b1_l*(Nbis^2)/(Nbis^2+Y1_l^2)) * 1.0/(1.0+exp(-0.464*(age-4.55))); else P_yl = 0.0;end;%%% Young Skuas
if (P_f>0.0000001)&&((tmod>=t_snowmelt)&&(tmod<t_fox_leaves)), age = 365*t-t_snowmelt; P_yf =  (b1_f*(Nbis^2)/(Nbis^2+Y1_f^2)) * 1.0/(1.0+exp(-0.36*(age-9.0))); 
else P_yf = 0.0;end;%%% Young Foxes (actualisation of lemming ensity at snowmelt)

%%% Predation Generalists
PredGen	=  0.0;%%% initializing to complete
PredGen = PredGen + P_o*W_o*(N(ki)^2)/(N(ki)^2+D_o^2) + P_f*W_f*(N(ki)^2)/(N(ki)^2+D_f^2) + P_l*W_l*(N(ki)^4)/(N(ki)^4+D_l^4); %% Adult predation
PredGen = PredGen + P_yo*W_o*(N(ki)^2)/(N(ki)^2+D_o^2) + P_yf*W_f*(N(ki)^2)/(N(ki)^2+D_f^2) + P_yl*W_l*(N(ki)^4)/(N(ki)^4+D_l^4); %% Young predation - to modify
%%% Predation stoat
PredSpe    = c*P(ki)*(N(ki)^2)/(N(ki)^2+D^2); 

%%% Stoat decline
Delta = 0.5 + (1.0/pi)*atan(b*(N(ki)-N_crit));
PredGrowth =  -P(ki)*d_h + Delta*P(ki)*(d_h-d_l); 

%----------------------------------------------------%
%%% This is Euler-Maruyama (simplest) scheme for SDE
%%% Multiplication by abundance for environmental stochasticity (instead of sqrt(abundance) for dem. stoch. see S. Engen papers)
% Env. stochasticity
StochCompPrey = sigmaE*N(ki)*sqrt(dt)*randn; % Standard deviation of the Wiener process increment is dt
StochCompPred = sigmaE*P(ki)*sqrt(dt)*randn; % We use the same sigmaE for parsimony, though there's no other reason
% Dem. stochasticity 
%DemogCompPrey = sqrt(abs(PreyGrowth) + PredGen + PredSpe)*sqrt(dt)*randn; % Demographic stochasticity - maybe problematic on density rather than abundance
%%DemogCompPred = sqrt(s*p(ki))*sqrt(dt)*randn; %%% sqrt(abs(PredGrowth)) a bit dirty - better to separate birth death rates, and sum them... but there we mostly have birth
%%% Some problems with ddem. stochasticity in the predator - sqrt(p) with p<1 yields large fluctuations, problems of scaling...

%%% Next step
N(ki+1) = N(ki) + dt*(PreyGrowth - PredGen - PredSpe) + StochCompPrey;%+ DemogCompPrey; % Only env. stoch. on the prey
P(ki+1) = P(ki) + dt*PredGrowth + StochCompPred ; % + DemogCompPred Put only dem. stoch. on the predator, ideally

%%% Add here the increase in Stoat Growth
if (stoat_born==0)&&(tmod>=t_stoat_born), 
	P(ki+1) = P(ki+1)*(1 + v);
	stoat_born=1;
	%disp(N1)
	end;

%%% Recording
Po(ki+1)=P_o;
Pf(ki+1)=P_f;
Pl(ki+1)=P_l;
Nb(ki+1)=Nbis;

%%% Increment time variables
ki=ki+1;
t=t+dt;
end

%--------------- End of loop on time -----------------------------------%

%---------------- Plotting --------------------------------------% 
figure,
subplot(211)
plot(linspace(0,T_max,length(1:ki)),N(1:ki),'-','LineWidth',2)
ylabel('Lemming')
axis tight
subplot(212)
plot(linspace(0,T_max,length(1:ki)),P(1:ki),'-','LineWidth',2)
ylabel('Weasel')
xlabel(' Time (years) ')
axis tight
set(1,'DefaultTextFontSize',20)
set(1,'DefaultAxesFontSize',20)
print(figure(1),'-depsc2','-r300','TS_Gilg')

figure,
subplot(411)
plot(linspace(0,T_max,length(1:ki)),Po(1:ki),'-','LineWidth',1)
ylabel('Owl')
axis tight
subplot(412)
plot(linspace(0,T_max,length(1:ki)),Pl(1:ki),'-','LineWidth',1)
ylabel('Skua')
subplot(413)
plot(linspace(0,T_max,length(1:ki)),Pf(1:ki),'-','LineWidth',1)
ylabel('Fox')
subplot(414)
plot(linspace(0,T_max,length(1:ki)),Nb(1:ki),'-','LineWidth',1)
ylabel('Lemm. crit.')
xlabel(' Time (years) ')
axis tight
set(1,'DefaultTextFontSize',20)
set(1,'DefaultAxesFontSize',20)
print(figure(2),'-depsc2','-r300','TS_Gilg_gene')


figure,
subplot(211)
plot(linspace(0,T_max,length(1:ki)),N(1:ki),'-','LineWidth',1)
ylabel('Prey')
axis tight
subplot(212)
plot(linspace(0,T_max,length(1:ki)),P(1:ki),'-','LineWidth',1)
ylabel('Predator')
xlabel(' Time (years) ')
axis tight
set(2,'DefaultTextFontSize',20)
set(2,'DefaultAxesFontSize',20)
%print(figure(3),'-depsc2','-r300','TS_Gilg')

figure,
subplot(211)
semilogy(linspace(0,T_max,length(1:intervals:ki)),N(1:intervals:ki),'o-','LineWidth',1)
ylabel('Prey')
axis tight
subplot(212)
semilogy(linspace(0,T_max,length(1:intervals:ki)),P(1:intervals:ki),'o-','LineWidth',1)
ylabel('Predator')
xlabel(' Time (years) ')
axis tight
set(3,'DefaultTextFontSize',20)
set(3,'DefaultAxesFontSize',20)
%print(figure(4),'-depsc2','-r300','TS_Gilg_log')

figure,
plot(N(5*intervals:ki),P(5*intervals:ki),'o-','LineWidth',1)
xlabel('Prey')
ylabel('Predator')
axis tight
set(4,'DefaultTextFontSize',20)
set(4,'DefaultAxesFontSize',20)
%print(figure(5),'-depsc2','-r300','Orbits_Gilg')
%----------------------------------------------------------------%



