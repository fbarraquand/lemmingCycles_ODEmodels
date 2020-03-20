function dydt=func_barrow(t,y)   % function handle for ODE - model 5 of Turchin and Batzli 2001

global A B G_S G_W K_V K_M U_S U_W u_S u_W tau r_max R alpha

dydt = zeros(3,1);   
if mod(y(4),1)>tau, U=U_S;u=u_S;G=G_S; else U = U_W;u=u_W;G=G_W; end; % tau=5/6, before = winter, after = summer (til t=0). 
%%% R = r_max/(A-G); approx the same result
dydt(1) = U*(1 - y(1)/K_V) - A*y(1)*y(3)/(y(1)+alpha*y(2)+B);
dydt(2) = u*y(2)*(1 - y(2)/K_M) - alpha*A*y(2)*y(3)/(y(1)+alpha*y(2)+B);
dydt(3) = R*y(3)*(A*(y(1)+alpha*y(2))/(y(1)+alpha*y(2)+B) - G); %- C*y(2)*y(3)/(y(2)+D);
dydt(4) = 1;
