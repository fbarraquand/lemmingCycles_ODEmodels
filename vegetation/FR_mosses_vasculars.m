%%% Theoretical prediction

%%% Parameters
%%% Looks like this could be a type II ref w.r.t brown, w.r.t. collared this is likelty to be a synchrony effect...
n1=linspace(0.1,1000,50); % vasculars
n2=linspace(0.1,2000,50); % mosses
alpha = 0.1; % 0.1 to 1
A = 15;
B = 70;

%%% Functional response 
[x,y]=meshgrid(n1,n2); % x vasculars
Ne_vasculars=A*x./(x+alpha*y+B); 
Ne_moss = alpha*A*y./(x+alpha*y+B); 

percent_coll_theor = Ne_moss./(Ne_moss+Ne_vasculars);
mesh(n1,n2,percent_coll_theor)
mesh(n1,n2,0.5*ones(50,50))
xlabel('Vasculars density')
ylabel('Mosses density')
zlabel('% mosses')
title('Theoretical FR Turchin Batzli')
hold on
%colormap hot
mesh(n1,n2,percent_coll_theor)
hidden off
hold off


