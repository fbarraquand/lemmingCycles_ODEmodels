%%% Model of Turchin and Hanski (Am Nat) 
%%% Barraquand & Henden 2011
%%% Function for ODE integration

%%% Note ------------------------------------------------------------------- %%%
%%% Predation by generalist seasonal -> change from the original model
%%% Or we could add a parameter summer season length -> Maybe more parsimonious

function dydt = th97_seasonal_func(t,y)

global r K G H C D s Q g h d a

%%% Note: time zero is january 1st, therefore the indicator variable Winter = 2 at midwinter and 0 and midsummer (1 on average)
Winter = (1+cos(2*pi*t)); 
% tmod = mod(t,1.0);
% if (tmod>0.25)&&(tmod<0.75), Winter = 0; else Winter = 1; end;
% if (tmod>0.33)&&(tmod<0.74), Winter = 0; else Winter = 1; end; % Step function formulation 
%Winter = 0.5*(1-cos(pi*cos(pi*t))); % Other possibility, here winter is longer although we have no smoothing parameter to tune the duration
dydt = zeros(3,1);   
M1 = r*y(1);
M2 = g*y(1)/(y(1)^2+h^2); % g*(1-Winter)*y(1)^2/(y(1)^2+h^2)% check more abrupt functional response instead -> change the effect of seasonality?
M3 = a*y(2)/(y(1)+d);
dydt(1) = (r*(2-Winter) - M1 - M2 - M3)*y(1);  % Now 2-Winter as the max of this variable is 2 now (mean=1), so that 1-cos is in the parenthesis
%%% This is not perfect of course as this assumes winter/summer are of equal length.
dydt(2) = s*(2-Winter)*y(2)-s*y(2)^2/y(1);
dydt(3) = 1;




