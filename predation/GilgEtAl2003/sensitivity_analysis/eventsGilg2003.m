%- State event finder for reprod of stoat (v=4)  -%
function [value,isterminal,direction]=eventsGilg2003(t,y)

value= mod(y(3),1)-(5*30+15)/365;       % value declines through zero at onset of winter, 
isterminal= 1;             % stop computation at event
direction= 1;              % 1 means that event only occur when value is increasing when passing zero
%- END State event finder -%