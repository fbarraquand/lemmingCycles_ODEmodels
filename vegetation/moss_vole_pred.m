function dydt=moss_vole_pred(t,y)   % function handle for ODE

global A B C D G K_M u R s Q

dydt = zeros(4,1);   
dydt(1) = u*y(1)*(1 - y(1)/K_M) - A*y(1)*y(2)/(y(1)+B);
dydt(2) = R*y(2)*(A*y(1)/(y(1)+B) - G) - C*y(2)*y(3)/(y(2)+D);
dydt(3) = s*y(3)*(1-Q*y(3)/y(2));
dydt(4) = 1;

