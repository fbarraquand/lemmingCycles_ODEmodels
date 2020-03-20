%%% Model of Turchin and Hanski (Am Nat) adapted to lemmings so that we can model the Traill island population with it
%%% Barraquand & Henden 2011
%%% Function for ODE integration

%%% Note ------------------------------------------------------------------- %%%
%%% Predation by generalist seasonal -> change from the original model
%%% Or we could add a parameter summer season length -> Maybe more parsimonious

function dydt = th_lemm_seasonal_func(t,y)

global r_min r_max r K G H C D s Q

%%% Note: time zero is january 1st, therefore the indicator variable Winter = 1 at midwinter and 0 and midsummer (maybe we can change that)
Winter = 0.5*(1+cos(2*pi*t)); %
% tmod = mod(t,1.0);
% if (tmod>0.25)&&(tmod<0.75), Winter = 0; else Winter = 1; end; % Here it means the parameters are specified as max, not average values (for averages, there should be Winter between 0 and 2)
 %if (tmod>0.33)&&(tmod<0.74), Winter = 0; else Winter = 1; end; % Step function formulation 
%Winter = 0.5*(1-cos(pi*cos(pi*t))); % Other possibility, here winter is longer although we have no smoothing parameter to tune the duration
dydt = zeros(3,1);   
r = (r_max-r_min)*Winter + r_min;
dydt(1) = r*y(1) - r_max*y(1)^2/K - G*(1-Winter)*y(1)^4/(y(1)^4+H^4) - C*y(1)^2*y(2)/(y(1)^2+D^2);  % - C*y(1)*y(2)/(y(1)+D); % type II functional response instead / in the regulation term we could have r instead of r_max
%%% There is less lemming growth in summer (zero in midsummer here), and more generalist predation in summer than in winter where it is null. 
%%% This is not perfect of course as this assumes winter/summer are of equal length.
dydt(2) = s*(1-Winter)*y(2)-s*Q*y(2)^2/y(1);
dydt(3) = 1;




