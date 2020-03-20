function [value,isterminal,direction]=eventsParasite(t,y)

value = mod(y(3),1.0)-0.3;   % stops in spring
isterminal = 1;             % stop computation at event
direction = 1;              % 1 means that event only occur when value is increasing when passing zero

end
