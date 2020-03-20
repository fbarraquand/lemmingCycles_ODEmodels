function dydt=func_mod5(t,y)   % function handle for ODE - model 5 of Turchin and Batzli 2001

global A B G K U_S U_W tau R

dydt = zeros(3,1);   
if mod(y(3),1)>tau, U=U_S; else U = U_W;end; 
dydt(1) = U*(1 - y(1)/K) - A*y(1)*y(2)/(y(1)+B);
dydt(2) = R*y(2)*(A*y(1)/(y(1)+B) - G); %- C*y(2)*y(3)/(y(2)+D);
dydt(3) = 1;
