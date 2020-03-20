function dydt=tri_food_chain(t,y)   % function handle for ODE

global A B G K U R

dydt = zeros(4,1);   
dydt(1) = U*(1 - y(1)/K) - A*y(1)*y(2)/(y(1)+B);
dydt(2) = R*y(2)*(A*y(1)/(y(1)+B) - G) - C*y(2)*y(3)/(y(2)+D);
dydt(3) = s*y(3)*(1-Q*y(3)/y(2));
dydt(4) = 1;

