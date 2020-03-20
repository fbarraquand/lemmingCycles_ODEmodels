function dydt = SI_func(t,y)

global r beta N_thresh gamma mu K

N = y(1)+y(2); % Total density
%InfectionRate = beta*y(2); % mass action or 'density-dependent transmission (McCallum et al TREE 2001)'
InfectionRate = beta*y(2)*N^gamma/(N_thresh^gamma+N^gamma); % Infection coefficient that has a treshold density
%%% The rationale is that mass action happen at large densities, but otherwise we have infections within population clusters that do not affect the whole population, and the percentage of connected clusters increase steeply after a threshold large-scale density. At low densities the rate of infection is very very low. The threshold is modelled through a sigmoid function, and otherwise this follows the rule of mass-action.
%k = 10;
%InfectionRate = k*log(1+beta*y(2)/k); % Negative binomial, aggregated infections
%InfectionRate = beta*10*y(2)/(50+y(2)); % Saturating infection rate -> at some point adding more infective individuals does not increase the likehood of getting infected. Unlikely I think. 
%InfectionRate = beta*10*y(2)/(y(1)+y(2)); % Depend on the proportion of infected individual (connection to the ratio-dependent predator-prey model)
%k0 = 1;% 0.1 % Here the spatial dispersion parameter k of the negative binomial is increasing with N (more well-mixed at high densities)
%InfectionRate = k0*N*log(1+beta*y(2)/(k0*N)); % Modified negative binomial, aggregated infections

dydt = zeros(3,1);   
dydt(1) = r*(1-y(1)/K)*y(1) - InfectionRate*y(1);
dydt(2) = InfectionRate*y(1) - mu*y(2);
dydt(3) = 1;

