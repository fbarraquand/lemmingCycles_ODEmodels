clear all
%%%%%%%%%%% SIMULATION SOLUTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global r rw rs v dh dl b Wst bo bf bs Yo Yf byo bys byf Yyo Yys Yyf Wo Do Ws Ds Wf Df Dst Ndot     % parameters that can be varied

%%%%%%% NB: control+c = terminates simulation prosess %%%%%%%%%%%%%


%- PARAMETERS from fig 2. in Gilg et al, 2003 -% 
rw= 4.0;                  % growth rate of lemmings in winter
rs= 0.8;                % growth rate of lemmings in summer. max in summer =0.8 since no young mature, but not all Ad probably reproduce
v= 4;                   % number of weaned young (both sexes) produced per female per year by the stoat
% dh=4.5; dl=0.2; b=25; % parameters usen in fig. 2 in Gilg et al 2003
dh= 4;                  % mortality of stoat at high prey density (may change: 3.5, 4, 4.5)
dl= 0.1;                % mortality of stoat at low prey density (may change: 0.0, 0.1, 0.2)
b= 25;                  % value to obtain steep function in predator mortality function, delta
Wst=1000;               % consumption by stoat (=c)
Dst=0.1; %0.08;         % slope of the type 3 funct resp for the stoat, i.e. half-saturation constant
Ncrit=Dst;              % lemming density at witch stoat mortality is (dh+dl)/2
bo= 0.00366;            % max owl density (per ha)
bf= 0.0016;             % max fox density (per ha)
bs= 0.02;               % density of skua (per ha)
Yo= 2.86;               % Slope of num resp for adult owls
Yf= 11;                 % Slope of num resp for adult fox
byo= 0.011;             % max density for young owls (fledglings per ha)
bys= 0.016;             % max density for young skuas (fledglings per ha)
byf= 0.0028;            % max density for young foxes (weaned young per ha)
Yyo= 4;                 % slope of num resp for young owls
Yys= 6;                 % slope of num resp for young skuas
Yyf= 5.3;               % slope of num resp for young foxes
Wo= 4.7*365;            % max pred rate for the owls (lemming per day)
Do= 1.08;               % slope of type 3 funct resp (e=2) for the owl, i.e. half-saturation constant
Ws= 4.4*365;            % max pred rate for the skuas (lemming per day)
Ds= 2.2;                % slope of type 3 funct resp (e=4) for the skua, i.e. half-saturation constant
Wf= 3.8*365;            % max pred rate for the foxes (lemming per day)
Df= 0.13;               % slope of type 3 funct resp (e=2) for the fox, i.e. half-saturation constant
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%% SIMULATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% outname = 'ts-series';

% Parameter range for simulations: 7 values for each parameter, 10, 20, 30 % increase and decrease from estimated value
% rwval = [2.8 3.2 3.6 4.0 4.4 4.8 5.6];  % growthrate winter of lemming
rsval = [0.28 0.32 0.36 0.4 0.48 0.52 0.56];  % growthrate summer of lemming
% vval = [2.8 3.2 3.6 4 4.4 4.8 5.2];   % numb weaned young per female stoat
% dhval = [2.8 3.2 3.6 4 4.4 4.8 5.2];  % mort of stoat at high prey densities
% dlval = [0.07 0.08 0.09 0.1 0.11 0.12 0.13];  % mort of stoat at low prey densities
% Wstval = [700 800 900 1000 1100 1200 1300];   % yearly consumption for stoat
% boval = [0.002562 0.002928 0.003294 0.00366 0.004026 0.004392 0.004758];  % max owl density
% bfval = [0.00112 0.00128 0.00144 0.0016 0.00176 0.00192 0.00208];  % max fox density
% bsval = [0.014 0.016 0.018 0.02 0.022 0.024 0.026];  % max skua density
% Yoval = [2.002 2.288 2.574 2.86 3.146 3.432 3.718];  % slope of num resp for ad owls 
% Yfval = [7.7 8.8 9.9 11 12.1 13.2 14.3];  % slope of num resp for ad fox 
% byoval = [0.0077 0.0088 0.0099 0.011 0.0121 0.0132 0.0143];    % max density for young owls (fledglings per ha)
% bysval = [0.0112 0.0128 0.0144 0.016 0.076 0.0192 0.0208];     % max density for young skuas (fledglings per ha)
% byfval = [0.00196 0.00224 0.00252 0.0028 0.00308 0.00336 0.00364]; % max density for young foxes (weaned young per ha)
% Yyoval = [2.8 3.2 3.6 4 4.4 4.8 5.2];  %   slope of num resp for young owls
% Yysval = [4.2 4.8 5.4 6 6.6 7.2 7.8];  %   slope of num resp for young skuas
% Yyfval = [3.71 4.24 4.77 5.3 5.83 6.36 6.89];  % slope of num resp for young foxes
% Woval = [3.29*365 3.76*365 4.23*365 4.7*365 5.17*365 5.64*365 6.11*365];   % max pred rate for the owls (lemming per day)
% Doval = [1.08]; % slope of type 3 funct resp (e=2) for the owl, i.e. half-saturation constant
% Wsval = [4.4*365];  % max pred rate for the skuas (lemming per day)
% Dsval = [2.2];                % slope of type 3 funct resp (e=4) for the skua, i.e. half-saturation constant
% Wfval = [2.66*365 3.04*365 3.42*365 3.8*365 4.18*365 4.56*365 4.94*365];   % max pred rate for the foxes (lemming per day)
% Dfval = [0.091 0.104 0.117 0.13 0.143 0.156 0.169];    % slope of type 3 funct resp (e=2) for the fox, i.e. half-saturation constant
Dstval = [0.056 0.064 0.072 0.08 0.088 0.096 0.104];    % half-sat constant of stoat funct resp

% N0=[0.1:0.1:0.5];                               % varying initial conditions of lemming
% P0=[0.001:0.001:0.005];                         % varying initial cond of stoat    

y1=-ones(1, 1);                                 % empty(-1) matrix for prey   
y2=-ones(1, 1);                                 % empty(-1) matrix for stoat   
y3=-ones(1, 1);                                 % empty(-1) matrix for time   
lold=zeros(1, 1);                               % empty matrix for lold (to keep track of time in iterations
% rs=0.8;                                       % value (initial) of growth-rate summer of lemming
% rw=4.0;                                       % value (initial) of growth-rate winter of lemming
corvar=zeros(7,7);  
%- loops -%
for i = 1:1:7                                   % first for-loop for rw
figure(1)                                       % pre-making figure box 
hold on    
for j = 1:1:7                                   % second for-loop for rs

% s=5*(i-1)+j;
s = 7*(i-1)+j;                                  % dummy, counting from 1 to (i*j)

%- initial conditions -%
tstart = 0;                                     % start time
tstop = 50;                                     % end time 
% y0=[N0(i) P0(j) 0];                           % Sim: initial values

%- options-handle -%
options = odeset('Events',@eventsGilg2003,'Reltol',1e-4);        % eventually: 'Jacobian',@Jac,'Vectorized','on');  % setting the options-structure for sertain events happening: Stochasticity in predator
y0 = [0.1 0.001 0];                                               % initial values
y1(1,s) = y0(1);                                                  % making sure that each iteration starts with appropiate initial value for prey
y2(1,s) = y0(2);                                                  % pred iterations start with appropiate initial value
y3(1,s) = y0(3);                                                  % appropiate initial value for time
yout = y0.';tout = 0; teout = []; yeout = []; ieout = [];                 % output of events, teout=time when events occur, yeout=values of species abundance at events, ieout=
lold(1,s) = 2;                                                    % counting timesteps in each simulation

%- Parameters for simulation:
rs = rsval(i);                                                    % first parameter to loop
Dst = Dstval(j);                                                    % second parameter to loop                                             % setting global parameters, i.e. those to be varied
   
while tout(length(tout))<tstop                                      % only iterate for time less than tstop (=20)
    [t,y,TE,YE,IE] = ode45(@gilg2003,[tstart tstop],y0,options);    % ode solver for function = gilg2003 with options statement
    time = length(t);                                                 % length of time (t)
    tout = [tout;t(2:time)];                                          % time-output, number of time-steps
    yout = [yout;t(2:time)];                                          % output of y follows the time-steps
    teout=[teout;TE]; yeout=[yeout;YE]; ieout=[ieout;IE];           % when events happen and values
    y1((lold(1,s)):(lold(1,s)+time-2), s) = [y(2:time,1)];          % setting up matrices for y1 (prey) for each iteration
    y2((lold(1,s)):(lold(1,s)+time-2), s) = [y(2:time,2)];          % setting up matrices for y2 (stoat) for each iteration
    y3((lold(1,s)):(lold(1,s)+time-2), s) = [y(2:time,3)];          % setting up matrices for y3 (time) for each iteration
    y0 = [y(time,1);y(time,2)*v;y(time,3)];                           % vector of initial values of prey, pred and time 
%     test = [y(1:time,2);v;y0(2)];
    tstart = t(time);                                               % length of time-steps (tstop)
    lold(1,s) = time+lold(1,s)-2;                                   % matrix of ??

% Ytls = [y3(:,s) y1(:,1:s) y2(:,1:s)];    % making matrix of timeseries w/1st col=Time, 2nd:10th col=lemming, 11th:19th col=stoat
Time = [y3(:,1:s)]; % matrix of tout of each sim.
Y = [y1(:,1:s)];    % matrix of Lemming density for each sim.
S = [y2(:,1:s)];    % matrix of Stoat density for each sim.
% xlswrite(outname,Y);  % Writing timeseries to exel work-sheet                                             % write to Excel-file

correlation = corrcoef(y1(1:lold(1, s), s),y2(1:lold(1, s), s));     % pairwise correlations of prey vs pred                                       % correlation between pred (stoat) and prey (lemming)
corvar(i,j) = correlation(1,2);                                      % diagonal matrix of correlations 
                                                 
%- Plot of log-transformed timeseries -%
% subplot(3,3,s);
% [AX,H1,H2] = plotyy(y3(1:lold(1, s), s),y1(1:lold(1, s), s),...
% y3(1:lold(1, s), s),y2(1:lold(1, s), s),'semilogy');
% title({rs Dst});                  
% set(get(AX(1),'Ylabel'),'String','Lemming Density','Color','k')
% set(get(AX(2),'Ylabel'),'String','Stoat Density','Color','r')
% set(AX,'Xlim',[0 tstop])
% set(AX(1),'YColor','k','Ylim',[0 100])
% set(AX(2),'YColor','r','Ylim',[0 1])
% set(H1,'LineStyle','-','Color','k','LineWidth',.5)
% set(H2,'LineStyle','-','Color','r','LineWidth',.7)

% - Plot true timeseries (no log) -%
subplot(7,7,s);                                                     % setting up grid of subplots
plotyy(y3(1:lold(:, s), s),y1(1:lold(:, s), s),...
y3(1:lold(:, s), s),y2(1:lold(:, s), s))                            % plotting prey & pred vs time for each iteration/loop
title({rs Dst});   

end;                                                                % end of first for-loop (i)                      
end;                                                                % end of second for-loop (j)
end;                                                                % end for while-statement

%- if exp growth: setting values over 200 to 100 & 10 for prey & pred, respectively
% highy1=find(y1>200);
% y1(highy1)=100
% highy2=find(y2>200);
% y2(highy1)=10;

%- Finding max and min values of timeseries:

Year1 = ceil(y3+0.25)-1;                    % rounding up for max values each year of Lemming (tuned)
Year2 = ceil(y3+0.02)-1;                    % rounding up for max values each year of Stoat (tuned)
Year3 = floor(y3-0.4)+1;                    % rounding down for min values each year of Lemming (tuned)
Year4 = floor(y3-0.4)+1;                    % rounding down for min values each year of Stoat (tuned) 

% count=1:length(Time);
% k = repmat(count.',[1 s]); 

%- empty matrices for max values:
y1max   = zeros(max(tstop),s);
y2max   = zeros(max(tstop),s);
y3maxy1 = zeros(max(tstop),s);
y3maxy2 = zeros(max(tstop),s);

yearmax = zeros(max(tstop),s);
%- empty matrices for min values:
y1min   = zeros(max(tstop),s);
y2min   = zeros(max(tstop),s);
y3miny1 = zeros(max(tstop),s);
y3miny2 = zeros(max(tstop),s);

%- loops for min & max values for each year per timeseries: 
% q=[1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20]; % NB: denne må kjøres for tstop>100, da flere komb. blir NAN etter det. Dvs. exp vekst => dimensjonene blir feil!!!
for p=1:s 
     for i=1:tstop
y1max(i,p)   = max(y1(Year1(:,p)==i,p));
y2max(i,p)   = max(y2(Year2(:,p)==i,p));                % extracting max values each year for Stoat
y3maxy1(i,p) = y3(y1==max(y1(Year1(:,p)==i,p)));        % extracting time values for each max value of Lemming
y3maxy2(i,p) = y3(y2==max(y2(Year2(:,p)==i,p)));        % extracting time values for each max value of Stoat
yearmax(i,p) = max(Year1(Year1(:,p)==i,p));             % every whole year
y1min(i,p)   = min(y1(Year3(:,p)==i,p));                % min values of Lemming
y2min(i,p)   = min(y2(Year4(:,p)==i,p));                % min values of Stoat
y3miny1(i,p) = y3(y1==min(y1(Year3(:,p)==i,p)));        % time value of each min of Lemming
y3miny2(i,p) = y3(y2==min(y2(Year4(:,p)==i,p)));
end;
    end; 
%Ymax = [y1max y2max y3maxy1 y3maxy2 yearmax]; % collecting results for max values
%Ymin = [y1min y2min y3miny1 y3miny2 yearmax]; % collecting results for min values


%- finding the value of amps higher than 4:
% maxy1 = zeros(length(tstop),s);
% for x=1:s
% for i=1:length(tstop)
% [i j]=find(y1max(:,x)>4);
% maxy1(i,x) = y1max(i,x);
% end; 
% end;





%- plotting max amplitude superimposed on timeseries to see accuracy:
for w=1:s  
figure(2)
    subplot(7,7,w)
    plot(y3(:,w),y1(:,w),'-r'), hold on
    plot(y3maxy1(:,w),y1max(:,w),'--d')
figure(3)
    subplot(7,7,w)
    plot(y3(:,w),y2(:,w),'-r'), hold on
    plot(y3maxy2(:,w),y2max(:,w),'--d')
figure(4)
    subplot(7,7,w)
    [AX,H1,H2] = plotyy(y3maxy1(:,w),y1max(:,w),y3maxy2(:,w),y2max(:,w),'semilogy');
    title('Max amplitude of Lemming & Stoat [per Year]');
    set(get(AX(1),'Ylabel'),'String','Lemming Density','Color','k')
    set(get(AX(2),'Ylabel'),'String','Stoat Density','Color','r')
    set(AX,'Xlim',[0 tstop],'XTick', (0:1:tstop))
    set(AX(1),'YColor','k','Ylim',[0 100])
    set(AX(2),'YColor','r','Ylim',[0 1])
    set(H1,'LineStyle','-','Marker','*','Color','k','LineWidth',.7)
    set(H2,'LineStyle','-','Marker','d','Color','r','LineWidth',.7)
figure(5)
    subplot(7,7,w)
    plot(y1max(:,w)),hold on
    plot(y1max(:,w),'r--')
end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%- Plotting ACF & Partial ACF for sim amp series -%
for z=1:s
figure(6)
    subplot(7,7,z)
    plot(pacf2(y1max(:,z),15,1))    % partial ACF for prey amp series
figure(7)
    subplot(7,7,z)
    plot(pacf2(y2max(:,z),15,1))    % partial ACF for stoat amp series
figure(8)
    subplot(7,7,z)
    plot(acf2(y1max(:,z),3))    % ACF for prey amp series
end;

%- Correlations between prey & pred timeseries, both true series and amp series:
corvar2=zeros(7,7);
corvar3=zeros(7,7);
for q=1:s
correlation2 = corrcoef(y1(:,q),y2(:,q));           % correlation based on raw timeseries    
correlation3 = corrcoef(y1max(:,q),y2max(:,q));     % correlation based on max amp timeseries
corvar2(q) = correlation2(1,2);    
corvar3(q) = correlation3(1,2);  % output: element 1.1 =1, 2.1=2, 1.2=4 osv.
end; 

%- Different ACF function:
r=zeros(1,s);
for a=1:s
l = floor(length(y3maxy1(:,1)/2)); 
p = tstop;
r=xcorr(y1max(:,a),p,'coeff'); % autocorrelation function  
d=(-l:l);          % times of delays
figure(9)
subplot(7,7,a);
plot(d,r);
% legend('Autocorrelation');
xlabel('Delay (s)');
ylabel('Correlation coeff.');
line(-tstop/2:tstop/2,0.5)
end;

