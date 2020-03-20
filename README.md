## Lemming population cycles models in Matlab

This repository contains Matlab codes for numerical integration of (mostly) ODE models describing lemming population dynamics. Namely, we consider models implementing different cycling mechanisms: 

- ``predation`` models, inspired by the case study of Gilg et al. (2003) at Traill Island, Greenland. 
- ``vegetation`` models, from Turchin & Batzli (2001) at Point Barrow, Alaska. 
- ``host-parasite`` models. 

*Nota bene* Although it is convenient for modellers to implement these mechanisms in separate models, there is of course no guarantee that they do not all interact together in real life. They might interact with stage structure too, which is neglected in most models here. 

### Practical notes 

- Most codes have been written in 2011 - 2014 by Frédéric Barraquand based on previous codes of the Traill island model by John-André Henden. We initially thought to use these for a publication that has not materialized (yet), hence the public release. 
- The codes use Matlab's ``ode45`` algorithm for integration. Although satisfactory in many cases, it can be inappropriate for some stiff problems. In other cases, more brute force methods like Runge-Kutta methods will be sufficient; to add noise (stochastic differential equations), Euler is best. We prefer to avoid discontinuities of functions with respect to time, which has guided some model simplifications.
- Send feeback to ``frederic.barraquand@u-bordeaux.fr`` 




