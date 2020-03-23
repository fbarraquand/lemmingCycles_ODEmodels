## Models of lemming population cycles driven by predation 

The models presented here are based on the case study of lemming population cycles at Traill island, east Greenland (although they might apply with some modifications to other places as well). The original model for this case study has been published in Gilg, O., Hanski, I., & Sittler, B. (2003). Cyclic dynamics in a simple vertebrate predator-prey community. Science, 302(5646), 866-868. This case study has one prey (the collared lemmings) and four predators, only one of which (the stoat) is modelled dynamically. The three `generalist' predators are the snowy owl, the long-tailed skua, and the Arctic fox. 

We consider a number of variants to the original model explained below. 

A major deviation from previous literature on rodent population dynamics pertains to the "inverted reproduction schedule" of lemmings: maximal rates of intrinsic growth for lemmings are in winter (they reproduce well under the snow), while population growth is (here) lower in summer. This allows seasonal generalist predation, occurring in the summer when the avian predators arrive, to have a large impact on lemming dynamics; it functions as seasonal, fluctuation-enhancing perturbation. In that sense, the models presented here differ from the classics like Turchin & Hanski American Naturalist 1997, in which (1) generalist predation is not seasonal and (2) rodent reproduce better in summer. 

We present three types of models: 
- The Gilg et al. (2003) model plus a variant with a different, Leslie-like numerical response (which produces nonetheless similar results). 
- ``Pooled generalists`` models, where the four generalist predators are aggregated in one compartment (their behaviour resembles then that of the long-tailed skua, which is the most abundant predator in Greenland). These are parsimonious simplifications of the original Traill model. 
- In the folder ``comparison_TH97``, we simplify modify the Turchin & Hanski (1997) model by (1) adding the possibility for generalist predation to be seasonal as well and (2) inverting the reproduction schedule (better reproduction in winter). This allows to visualize the dynamical implifications of seasonal generalist predation and the reproduction schedule. 

### What diagnostics can you find in these folders? 
- time series
- mortality rates inflicted by the various predators. For instance, an inverted reproduction schedule (with more reproduction in winter) combined with a large seasonal generalist predation in summer tend to make most of the summer mortality (and hence, the population declines) attributable to generalist predators instead of the specialist predator. 
