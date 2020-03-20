%%% Model of Turchin and Hanski (Am Nat 1997)- seasonal version with breakdown of mortality rates
%%% Barraquand & Henden 2011
%%% Seasonal version
%%% Main thing yet to do: it needs a longer winter period and a smaller summer period...


%%% NB: Threshold formulation of seasonality more prone to generate high generalist predation rates, seasonally
%%% Maybe look with lower specialist predator fertilities - done for the inversed repro schedule, does not seem to change everything
%%% Also, we might compare to cases without generalist predators 
%%% Note: Figs 3 to 4 shows though seasonality of gen. predation has no impact when birth rates are already seasonal, it can have an important effect when they are not. 
%%% Maybe we should look at this with less seasonal birth rates (not zero in winter...)

clc 
clear all
close all

%%%%%% Parameters
global r K G H C D s Q g h d a
r = 6; % twice that = max value
e=1.0;
K=150.0; % between 100 and 300... 
s=1.25; %from 1 to 1.6 (s=0.75 - lower cycles)
C=600.0;
D=6.0; %(very low D/K!) 
Q=40.0; % 
%%%%%%
G=60.0; %140 or 123.0; 10 = value used first
H=15; %13.5;

%--------------- Scaling ----------------------------------------%
g=G/K
h=H/K
d=D/K
a=C/Q

y0=[0.1 0.001 0];
tstart = 0.0;
tstop = 50.0;

tsampling = 501;

%%%% First baseline case
% Generalist predators constant and seasonality in growth rates

tspan=linspace(tstart,tstop,tsampling);                        % timespan for the numerical solution
options= odeset('Reltol',1e-3,'NonNegative',[1 2]);
[tout,yout] = ode45(@th97_seasonal_func, tspan , y0,options);       % ode solver 

MReg =  r*yout(:,1);
MGen =  g*yout(:,1)./(yout(:,1).^2+h^2);% g*(1-Winter)*y(1)^2/(y(1)^2+h^2)% check more abrupt functional response instead -> change the effect of seasonality?
MSpe = a*yout(:,2)./(yout(:,1)+d);
Mtot=MReg+MSpe+MGen;

% Plotting
subplot(311)
plot(tout,yout(:,1),'LineWidth',2)
ylabel('Vole')
subplot(312)
plot(tout,yout(:,2),'LineWidth',2)
ylabel('Weasel')
subplot(313)
plot(tout,MReg./Mtot,'b',tout,MGen./Mtot,'k',tout,MSpe./Mtot,'r','LineWidth',2)
axis([tstart tstop 0 1])
ylabel('Variable mortality (fraction)')
xlabel('Time')
title('Seasonal r and constant G')
% Plot CV, S-axis, max/min (needs annual data)

% Fourier analysis (also annual? or semi...) TO DO...
%{ 
% Matlab example
L=length(yout(:,1));
Fs = 1;                    % Sampling frequency
NFFT = 2^nextpow2(L); % Next power of 2 from length of y
Y = fft(yout(:,1),NFFT)/L;
f = Fs/2*linspace(0,1,NFFT/2+1);

% Plot single-sided amplitude spectrum.
plot(f,2*abs(Y(1:NFFT/2+1))) 
title('Single-Sided Amplitude Spectrum of y(t)')
xlabel('Frequency (Hz)')
ylabel('|Y(f)|')
%}


figure,

%%%% Second case
% Generalist predators variable and seasonality in growth rates

tspan=linspace(tstart,tstop,tsampling);                        % timespan for the numerical solution
options= odeset('Reltol',1e-3,'NonNegative',[1 2]);
[tout,yout] = ode45(@th97_seasonal_gen_func, tspan , y0,options);       % ode solver 

Winter = 1+cos(2*pi*tout);
%{ 
% Threshold version of seasonality
Winter=2*ones(length(tout),1);
 tmod = mod(tout,1.0);
Winter((tmod>0.25)&(tmod<0.75))=0;
%}
MReg =  r*yout(:,1);
MGen =  g*yout(:,1).*(2-Winter)./(yout(:,1).^2+h^2);% g*(1-cos)*y(1)^2/(y(1)^2+h^2)% check more abrupt functional response instead -> change the effect of seasonality?
MSpe = a*yout(:,2)./(yout(:,1)+d);
Mtot=MReg+MSpe+MGen;

% Plotting
subplot(311)
plot(tout,yout(:,1),'LineWidth',2)
ylabel('Vole')
subplot(312)
plot(tout,yout(:,2),'LineWidth',2)
ylabel('Weasel')
subplot(313)
plot(tout,MReg./Mtot,'b',tout,MGen./Mtot,'k',tout,MSpe./Mtot,'r','LineWidth',2)
axis([tstart tstop 0 1])
ylabel('Variable mortality (fraction)')
xlabel('Time')
title('Seasonal r and seasonal G')

figure,

%%%% Third case
% Generalist predators variable and no seasonality in growth rates

tspan=linspace(tstart,tstop,tsampling);                        % timespan for the numerical solution
options= odeset('Reltol',1e-3,'NonNegative',[1 2]);
[tout,yout] = ode45(@th97_seasonal_gen_ng_func, tspan , y0,options);       % ode solver 

Winter = 1+cos(2*pi*tout);
%{ 
% Threshold version of seasonality
Winter=2*ones(length(tout),1);
 tmod = mod(tout,1.0);
Winter((tmod>0.25)&(tmod<0.75))=0;
%}

MReg =  r*yout(:,1);
MGen =  g*yout(:,1).*(2-Winter)./(yout(:,1).^2+h^2);% g*(1-cos)*y(1)^2/(y(1)^2+h^2)% check more abrupt functional response instead -> change the effect of seasonality?
MSpe = a*yout(:,2)./(yout(:,1)+d);
Mtot=MReg+MSpe+MGen;

% Plotting
subplot(311)
plot(tout,yout(:,1),'LineWidth',2)
ylabel('Vole')
subplot(312)
plot(tout,yout(:,2),'LineWidth',2)
ylabel('Weasel')
subplot(313)
plot(tout,MReg./Mtot,'b',tout,MGen./Mtot,'k',tout,MSpe./Mtot,'r','LineWidth',2)
axis([tstart tstop 0 1])
ylabel('Variable mortality (fraction)')
xlabel('Time')
title('Constant r and seasonal G')

%%%% Fourth case
% Generalist predators constant and no seasonality in growth rates
figure,

tspan=linspace(tstart,tstop,tsampling);                        % timespan for the numerical solution
options= odeset('Reltol',1e-3,'NonNegative',[1 2]);
[tout,yout] = ode45(@th97_not_seasonal_func, tspan , y0,options);       % ode solver 

MReg =  r*yout(:,1);
MGen =  g*yout(:,1)./(yout(:,1).^2+h^2);% g*(1-cos)*y(1)^2/(y(1)^2+h^2)% check more abrupt functional response instead -> change the effect of seasonality?
MSpe = a*yout(:,2)./(yout(:,1)+d);
Mtot=MReg+MSpe+MGen;

% Plotting
subplot(311)
plot(tout,yout(:,1),'LineWidth',2)
ylabel('Vole')
subplot(312)
plot(tout,yout(:,2),'LineWidth',2)
ylabel('Weasel')
subplot(313)
plot(tout,MReg./Mtot,'b',tout,MGen./Mtot,'k',tout,MSpe./Mtot,'r','LineWidth',2)
axis([tstart tstop 0 1])
ylabel('Variable mortality (fraction)')
xlabel('Time')
title('Constant r and constant G')


