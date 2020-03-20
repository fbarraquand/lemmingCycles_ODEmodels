%%%% State event finder for the timing of reproduction of the stoat %%%%%
function [value,isterminal,direction]=events_gilg_03(t,y)

global t_stoat_born

value = mod(y(3),1.0)-t_stoat_born;      % value goes through zero at the moment of birth
isterminal= 1;             		 % stop computation at event
direction= 1;                            % 1 means that event only occur when value is increasing when passing zero
%- END State event finder -%
