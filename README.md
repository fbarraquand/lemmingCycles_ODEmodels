[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4271834.svg)](https://doi.org/10.5281/zenodo.4271834)

## Models for lemming population cycles (in Matlab)

This repository contains Matlab codes for numerical integration of (mostly) ODE models describing lemming population dynamics: 

- predator-prey models in the ``predation`` folder. These include the model of Gilg et al. (2003) for the lemming population of Traill Island, Greenland, as well as simpler variants. 
- consumer-plant models in the ``vegetation`` folder, after the paper by Turchin & Batzli (2001) which focuses explaning fluctuations at Point Barrow, Alaska. 
- ``host-parasite`` models. 

*Nota bene* Although it is convenient for modellers to implement these mechanisms in separate models, there is of course no guarantee that they do not all interact together in reality. They might interact with each other -- and with stage structure too, which is neglected in most models here. 

### Practical notes 

- Most codes have been initially written in 2011 - 2014 by Frédéric Barraquand (CNRS; then at University of Tromsø) based on previous codes of the Traill island model by John-André Henden (University of Tromsø). We initially thought to use these for a publication that has not materialized (yet), hence the public release. Our follow-up study on the Traill island model involved Frédéric Barraquand, John-André Henden, Olivier Gilg, Rolf A. Ims and Nigel G. Yoccoz. [Summary here](https://github.com/fbarraquand/lemmingCycles_ODEmodels/tree/master/predation/Traill_model_description.pdf). 
- The codes use Matlab's ``ode45`` algorithm for integration. Although satisfactory in many cases, it can be inappropriate for some stiff problems. In other cases, Runge-Kutta methods will be sufficient; to add noise (stochastic differential equations), Euler is usually best. We prefer in general to avoid discontinuities of functions with respect to time, which has guided some model simplifications.
- Please send feeback to ``frederic.barraquand@u-bordeaux.fr`` 
- If using these codes, please cite the associated Zenodo DOI as well as the original publications (depending on which model you are using):
     + Gilg, O., Hanski, I., & Sittler, B. (2003). Cyclic dynamics in a simple vertebrate predator-prey community. Science, 302(5646), 866-868
     + Turchin, P., & Batzli, G. O. (2001). Availability of food and the population dynamics of arvicoline rodents. Ecology, 82(6), 1521-1534.




