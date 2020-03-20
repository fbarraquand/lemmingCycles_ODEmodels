%%%% Function specifying the Gilg et al. 2003 Model for clean numerical integration with Matlab ode45
%%%% F. Barraquand and J.A. Henden. Last updated 01/10/2013
%%%% 20/11/2013 Leslie-like numerical response implemented 

function dydt=gilg_func_breakdown_Leslie(t,y)        % function handle for ODE
%%% List of global variables whose values are specified within the main script
global t_return_owl_fox t_return_skua t_snowmelt t_owl_hatch t_stoat_born t_skua_born t_skua_leaves t_owl_leaves t_fox_leaves r_w r_s v c D d_l d_h b W_o D_o W_l D_l W_f D_f b_o Y_o P_l_const b_f Y_f b1_o Y1_o b1_l Y1_l b1_f Y1_f Delta Nprime s q

%%% Relationship between parameters
N_crit = D;
%%%%%%%% STATE VARIABLES %%%%%%%%
%%% y(1) = N = lemming density
%%% y(2) = P = stoat density
%%% y(3) = time vector
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tmod=mod(y(3),1.0); % Variable between 0 and 1 indicating time of year

%%%%%% Determination of the the period of the year (winter or summer) %%%%%%%
%-------------------------------------------- WINTER VS SUMMER ------------------------------------------------%
if (tmod<t_snowmelt), 
	snowmelt = 0;
	r=r_w; %% winter growth rate of change for lemmings
	Nprime = y(1); %% Nprime=lemming density at snowmelt, replaced by N before snowmelt
else 
	snowmelt = 1; %% snowmelt has happened
	r=r_s;
	%Nprime = Nprime; %% Nprime=lemming density at snowmelt, keep it that way, passed as a global variable (not sure this is very clean...)
end; 
PreyGrowth = r*y(1); % exponential growth *(1-y(1)/K); if logistic growth needed, dont forget to parametrize K

%%% Numerical responses of generalists predators %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%------------------------------------------- Adult numbers -----------------------------------------------------------%
if (Nprime>2.0)&&((tmod>t_return_owl_fox)&&(tmod<t_owl_leaves)), P_o = (b_o*(Nprime-2.0))/(Y_o+Nprime-4.0); else P_o = 0.0; end;%%% Owls,  Nprime replaced by N before snowmelt (introduced variable N2 and if condition)
%P_o = 0.0;
if ((tmod>t_return_skua)&&(tmod<t_skua_leaves)), P_l = P_l_const; else P_l = 0.0;end;%%% Long tailed skua -> Should I put (Nprime>2.0)&& before like JA ??????
%P_l = 0.0;
if (Nprime>2.0)&&((tmod>t_return_owl_fox)&&(tmod<t_fox_leaves)), P_f = (b_f*Nprime^2)/(Y_f^2+Nprime^2);  else P_f = 0.0;end;%%% Foxes, Nprime replaced by N before snowmelt %%% Weird expression in paper - recheck
%P_f=0.0; % Put the foxes at zero

%------------------------------------------ Young numbers -------------------------------------------------------------%
if (P_o>0.0000001)&&((tmod>t_owl_hatch)&&(tmod<t_owl_leaves)), age = 365*(t-t_owl_hatch); P_yo =  ((b1_o*(Nprime-2.0))/(Y1_o+Nprime-4.0)) * 1.0/(1.0+exp(-0.36*(age-9.0))); else P_yo = 0.0;end;%%% Young Owls 
if (P_l>0.0000001)&&((tmod>t_skua_born)&&(tmod<t_skua_leaves)), age = 365*(t-t_skua_born); P_yl =  (b1_l*(Nprime^2)/(Nprime^2+Y1_l^2)) * 1.0/(1.0+exp(-0.464*(age-4.55))); else P_yl = 0.0;end;%%% Young Skuas
if (P_f>0.0000001)&&((tmod>t_snowmelt)&&(tmod<t_fox_leaves)), age = 365*(t-t_snowmelt); P_yf =  (b1_f*(Nprime^2)/(Nprime^2+Y1_f^2)) * 1.0/(1.0+exp(-0.36*(age-9.0))); else P_yf = 0.0;end;%%% Young Foxes (actualisation of lemming density at snowmelt)
%P_yf = 0.0;% to check the effect of juvenile predation
%P_yo = 0.0;
%P_yl = 0.0;
%%% Functional responses %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%------------------------------------ Predation Generalists -----------------------------------------------------------%
PredOwlNow =(P_o+P_yo)*W_o*(y(1)^2)/(y(1)^2+D_o^2);
PredSkuaNow = (P_l+P_yl)*W_l*(y(1)^4)/(y(1)^4+D_l^4);
PredFoxNow = (P_f+P_yf)*W_f*(y(1)^2)/(y(1)^2+D_f^2);
% Persistent variables above
PredGen = PredOwlNow + PredSkuaNow + PredFoxNow; %% Adult+Young predation

%------------------------------------- Predation Stoat ----------------------------------------------------------------%
PredSpe    = c*y(2)*(y(1)^2)/(y(1)^2+D^2); 

%---------------------------------- Stoat mortality processes ----------------------------------------------------------%
Delta = 0.5 + atan(b*(y(1)-N_crit))/pi; % Delta = 1 Good conditions, 0 bad conditions
%PredGrowth =  -y(2)*d_h + Delta*y(2)*(d_h-d_l); 
if ((r==r_s)&&(Nprime>=N_crit)),
PredGrowth =  s*y(2)*(1-q*y(2)/y(1));
else
PredGrowth =  -s*y(2)*q*y(2)/y(1);
end;
%-y(2)*(d_h - Delta*(d_h-d_l)); % Perhaps could be multiplied by Delta
%%% Integration %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dydt = zeros(3,1);                                            % initializing matrix
dydt(1) = PreyGrowth  - PredGen - PredSpe;                     % Changes in lemming density 
dydt(2) = PredGrowth;                                         % Changes in stoat density
dydt(3) = 1;                                                  % Matlab rule - see documentation Ode45
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


