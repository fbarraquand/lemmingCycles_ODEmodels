%%%% State event finder for the timing of reproduction of the stoat %%%%%
function [value,isterminal,direction]=events_gilglike_lemm(t,y)

global t_repro

value = mod(y(3),1.0)-t_repro;           % value goes through zero at the time of stoat reproduction
isterminal= 1;             		 % stop computation at event
direction= 1;                            % 1 means that event only occur when value is increasing when passing zero
%- END State event finder -%
