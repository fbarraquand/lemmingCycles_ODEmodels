%%%% State event finder for the timing of reproduction of the stoat %%%%%
function [value,isterminal,direction]=events_barrow(t,y)

global tau

value = mod(y(3),1.0)-tau;      % value goes through zero at the summer-winter transition
isterminal= 1;             		 % stop computation at event
direction= 1;                            % 1 means that event only occur when value is increasing when passing zero
%- END State event finder -%
