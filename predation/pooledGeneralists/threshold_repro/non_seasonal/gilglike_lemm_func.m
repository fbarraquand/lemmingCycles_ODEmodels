%%% Mix between Turchin and Hanski (Am Nat 1997) seasonal generalist predation and Gilg model 
%%% Modified from dydt = th_lemm_starv_func(t,y)
%%% Note ------------------------------------------------------------------- %%%
%%% Threshold reproduction for the stoat 
%%% We need to add a parameter summer season length -> Maybe more parsimonious?

function dydt = gilglike_lemm_func(t,y)

global r_s r_w K G H C D t_repro d_high d_low

%%% Note: time zero is january 1st, therefore the indicator variable Winter = 1 at midwinter and 0 and midsummer (maybe we can change that)
 Winter = 0.5*(1+cos(2*pi*t)); 
 b=25;
 Delta = 0.5 + (1.0/pi)*atan(b*(y(1)-D));
% tmod = mod(y(3),1.0);
%if (tmod>0.25)&&(tmod<0.75), Winter = 0; else Winter = 1; end;
%if ((tmod>0.33)&&(tmod<0.74)), Winter = 0; else Winter = 1; end; % Step function formulation 
%Winter = 0.5*(1-cos(pi*cos(pi*t))); % Other possibility, here winter is longer although we have no smoothing parameter to tune the duration
dydt = zeros(3,1);   
dydt(1) = r_s*y(1) + (r_w-r_s)*y(1)*Winter - G*(1-Winter)*y(1)^4/(y(1)^4+H^4) - C*y(1)^2*y(2)/(y(1)^2+D^2);  % - C*y(1)*y(2)/(y(1)+D); % type I functional response instead 
%%% Regulation term - (Winter*r_w+(1-Winter)*r_s)*y(1)^2/K , shut off now
%%% There is less lemming growth in summer (zero in midsummer here), and more generalist predation in summer than in winter where it is null. 
%%% This is not perfect of course as this assumes winter/summer are of equal length.
%%% dydt(2) = s*(1-Winter)*y(2)-s*Q*y(2)^2/y(1);
%%% Changed stoat dynamics -> that of the Gilg model literally
dydt(2) = -y(2)*(d_high - Delta*(d_high-d_low));
dydt(3) = 1;




