function dy=gilg2003(t,y)        % function handle for ODE

global r rw rs v dh dl b Wst bo bf bs Yo Yf byo bys byf Yyo Yys Yyf Wo Do Ws Ds Wf Df Dst Ndot K       % parameters that can be varied
% profile on
%- PARAMETERS -% 
% rw= 4;                  % growth rate of lemmings in winter
% rs= 0.8;                % growth rate of lemmings in summer. max in summer =0.8 since no young mature, but not all Ad probably reproduce
% v= 4;                   % number of weaned young (both sexes) produced per female per year by the stoat
% dh=4.5; dl=0.2; b=25; % parameters usen in fig. 2 in Gilg et al 2003
% dh= 4.4;                % mortality of stoat at high prey density (may change: 3.5, 4, 4.5)
% dl= 0.2;                % mortality of stoat at low prey density (may change: 0.0, 0.1, 0.2)
% b= 25;                  % value to obtain steep function in predator mortality function, delta
% Wst=1000/365;           % consumption by stoat (=c)
% Dst=0.08;               % slope of the type 3 funct resp for the stoat, i.e. half-saturation constant
Ncrit=Dst;                % lemming density at witch stoat mortality is (dh+dl)/2
c=Wst*365;
% bo= 0.00366;            % max owl density (per ha)
% bf= 0.0016;             % max fox density (per ha), NB: endret i en revidert utgave av artikelen (0.0008 i org. artikkelen)
% bs= 0.02;               % density of skua (per ha)
% Yo= 2.86;               % Slope of num resp for adult owls
% Yf= 11;                 % Slope of num resp for adult fox
% byo= 0.011;             % max density for young owls (fledglings per ha)
% bys= 0.016;             % max density for young skuas (fledglings per ha)
% byf= 0.0028;            % max density for young foxes (weaned young per ha)
% Yyo= 4;                 % slope of num resp for young owls
% Yys= 6;                 % slope of num resp for young skuas
% Yyf= 5.3;               % slope of num resp for young foxes
% Wo= 4.7;                % max pred rate for the owls (lemming per day)
% Do= 1.08;               % slope of type 3 funct resp (e=2) for the owl, i.e. half-saturation constant
% Ws= 4.4;                % max pred rate for the skuas (lemming per day)
% Ds= 2.2;                % slope of type 3 funct resp (e=4) for the skua, i.e. half-saturation constant
% Wf= 3.8;                % max pred rate for the foxes (lemming per day)
% Df= 0.13;               % slope of type 3 funct resp (e=2) for the fox, i.e. half-saturation constant

                          % These generalists are not present before mai 1st and even then only when Ndot > 2.
%%%%%%%%%%%%%%%%%%%%%%%%%%% END %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%- SEASONAL INDEX FOR PRAY (summer/winter) AND PREDATORS (present/absent) -%

%- Lemming -%
if mod(y(3),1)>(5*30+14)/365 && mod(y(3),1)<(9*30)/365 
     lsummer.inx=1;                                                           % summer between may and september 
else lsummer.inx=0;                                                         % otherwise winter
end;
r=(rs*lsummer.inx + rw*abs(lsummer.inx-1));                                   % index for when to use rs and rw for lemming growth
%- Owls -%
if mod(y(3),1)>(4*30)/365 && mod(y(3),1)<(9*30)/365 
     aopres.inx=1;                                                            % adult owls present, given dates
else aopres.inx=0;                                                          % otherwise not present
end;
if mod(y(3),1)>(5*30+24)/365 && mod(y(3),1)<(9*30)/365 
     yopres.inx=1;                                                            % juvenile owls present, given dates
else yopres.inx=0;                                                          % otherwise not present
end;
%- Skua -%
if mod(y(3),1)>(5*30+4)/365 && mod(y(3),1)<(7*30+15)/365 
     aspres.inx=1;                                                            % adult skua present, given dates
 else aspres.inx=0;                                                          % otherwise not present
end;
if mod(y(3),1)>(6*30+13)/365 && mod(y(3),1)<(7*30+15)/365 
    yspres.inx=1;                                                            % juvenile skua present, given dates
 else yspres.inx=0;                                                          % otherwise not present
end;
%- Fox -%
if mod(y(3),1)>(4*30)/365 && mod(y(3),1)<(9*30+22)/365 
    afpres.inx=1;                                                           % adult fox present, given dates
else afpres.inx=0;                                                          % otherwise not present
end;
if mod(y(3),1)>(5*30+14)/365 && mod(y(3),1)<(9*30+22)/365 
    yfpres.inx=1;                                                            % juvenile fox present, given dates
else yfpres.inx=0;                                                          % otherwise not present   
end;     
% NB: these have to be incorporated in the numerical response equations for generalist predators   
%%%%%%%%%%%%%%%%%%%%%%%%%%%% END %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%- DAYS-functions for young generalist predators -%
Oyage= 365*(mod(y(3),1)-(5*30+24)/365);                                    % days from 25 of june for young Owls
Syage= 365*(mod(y(3),1)-(6*30+13)/365);                                    % days from 14 of july for young Skua
Fyage= 365*(mod(y(3),1)-(5*30+14)/365);                                    % days from 15 of june for young Fox
%%%%%%%%%%%%%% END %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%- Incorporating Ndot as a constant from snowmelt -%
if mod(y(3),1)<(5*30+15)/365 
    Ndot= y(1);
else Ndot=Ndot ;                                                            % Ndot, is a constant of y(1) at snowmelt (june 15) thrue the entire season
end;
%%%%%%%%%%%%%%% END %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%- NUMERICAL RESPONSES (gilg et al 2003_online material -%
if ((aopres.inx >= 1) && (Ndot>2))
% Pao = (2*bo)*(Ndot-2)/(Yo+Ndot-4);                                           % num resp of adult Owls (Ndot replaced with N before snowmelt, june 15)
Pao = bo*(Ndot-2)/(Yo+Ndot-4);         %endret til denne i revidert utgave
else Pao= 0;
end;
if ((yopres.inx > 0) && (Pao > 0))
Pyo = ((byo*(Ndot-2))/(Yyo+(Ndot-4)))*(1/(1+exp(-0.36*(Oyage-9))));           % num resp of young owls, as adult eqvivalents (multipied with "adult" factor), from june 25 to sept 30 when Po>0)  
else Pyo = 0;
end;
if ((aspres.inx >= 1) )    %&& (Ndot>=2))   %%% skuas are supposed to NOT be associated with lemmings at snowmelt!!!!!!!
Pas = bs;                                                                    % num resp of adult skua
else Pas = 0;
end;
if ((yspres.inx >= 1) && (Pas>0))
Pys = (bys*Ndot^2)/(Yys^2+Ndot^2)*(1/(1+exp(-0.464*(Syage-4.55))));           % num resp of young skua, from july 14 to august 15
else Pys = 0;
end;
if ((afpres.inx >= 1) && (Ndot>2))
% Paf = (0.04+(bf-0.04))*Ndot^2/(Yf^2+Ndot^2);                                 % num resp adult arctic fox, from may 1 to october 22, with Ne replaced with N before snowmelt (june 15)
Paf = (0.0004+(bf-0.0004))*Ndot^2/(Yf^2+Ndot^2); % endret til denne i revidert utgave
else Paf = 0;
end;
if ((yfpres.inx >= 1) && (Paf>0))
Pyf = (byf*Ndot^2)/(Yyf^2+Ndot^2)*(1/(1+exp(-0.464*(Fyage-4.55))));           % num resp of young fox, from june 15
else Pyf = 0;
end; 
%%%%%%%%%%%%%%%%%%%%%%%%%%% END %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%- FUNCTIONAL RESPONSES (eq. S4b in gilg et al 2003_online material -%
DCRowl= Wo*y(1)^2/(Do^2+y(1)^2);                                            % funct resp of Owl (Ad & Juv merged together)
DCRskua= Ws*y(1)^4/(Ds^4+y(1)^4);                                           % funct resp of Skua (Ad & Juv merged together)
DCRfox= Wf*y(1)^2/(Df^2+y(1)^2);                                            % funct resp of Fox (Ad & Juv merged together)
DCRstoat= Wst*y(1)^2/(y(1)^2+ Dst^2);                                       % daily lemming consumprion rate by one stoat?
%%%%%%%%%%%%%%%%%%%%%%%%% END %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%- PREY MORTALITY INFLICTED BY GENERALIST PREDATORS -%
mowl= ((Pao+Pyo)*DCRowl)/y(1);                                              % mortality inflicted by Owls (ad & Juv merged) 
%myowl= (Pyo/DCRowl)/y(1);                                                   % mortality inflicted by Young Owls
mskua= ((Pas+Pys)*DCRskua)/y(1);                                            % mortality inflicted by Skuas (ad & Juv merged)
%myskua= (Pys/DCRskua)/y(1);                                                 % mortality inflicted by Young Skuas
mfox= ((Paf+Pyf)*DCRfox)/y(1);                                              % mortality inflicted by Fox (ad & Juv merged)
%myfox= (Pyf/DCRfox)/y(1);                                                   % mortality inflicted by Young Fox
mstoat= (y(2)*DCRstoat)/y(1);
%%%%%%%%%%%%%%%%%%%%%%%%%% END %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%psi=sin(2*pi*y(3));                                        % sin function of time y(3)
%omega=(sign(psi)/2)*abs(psi)^a+0.5;                        % smooth-season function omega
delta=0.5+atan(b*(y(1)-Ncrit))/pi;                          % predator mortality function, delta

%- PREY MORTALITY INFLICTED BY STOAT -%
m2stoat = (dh-(dh-dl)*delta);                               % Stoat mortality rate
%%%%%%%%%%%%%%%%%%% END %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


dy = zeros(3,1);                                            % creating a zero_matrix for timeseries

%- EQUATION SET: PREY w/ GENERALIST PREDATORS & SPECIALIST PREDATOR -%
 
dy(1) = (r-mowl-mfox-mskua-mstoat)*y(1);                    % the prey equation with mortality inflicted by generalist predators
dy(2) = -y(2)*m2stoat;                                        % stoat equation
dy(3) = 1;                                                  % the TIME-function, derivative of 1 leads to continous time
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 



