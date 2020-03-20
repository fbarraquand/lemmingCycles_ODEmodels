## Variants of the classic Turchin-Hanski (1997) model with different seasonality scenarios

The original model is presented in 

Turchin, P., & Hanski, I. (1997). An empirically based model for latitudinal gradient in vole population dynamics. The American Naturalist, 149(5), 842-874.

In the ``with_and_without_seasonal_generalists`` subrepo we consider different forms of seasonality.
* With/without seasonality in the maximum (intrinsic) growth rate r (considered by TH97)
* With/without seasonality in the quantity of generalist predation G (which TH97 did not consider)

In the ``inverted_repro_schedule`` subrepo we consider what happens when -- unlike the original publication -- we allow rodents to be more efficient at reproducing in winter rather than in summer. This corresponds relatively well to the lemmings and allows comparison to the model of Gilg et al. (2003). Probably the maximum growth rate parameters would need to be changed a little to detect whether this model could show cycles where the generalist predators can be responsible for most of the mortality during the summer declines. 
