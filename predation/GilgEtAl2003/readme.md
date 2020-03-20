## Gilg et al. (2003) Traill island model and variants

In this repository, we provide 
* the ``original_model``. The model provide cycles of amplitude and periodicity similar to the data (the reader should note that despite this, while the data has a one-year increase and two-year decline, the model shows a multi-year increase and one-year sharp decline). 
* ``Leslie_numResp`` implements a Leslie-style numerical response instead of having the young stoats bursting out in the spring. 
* ``without_skuas`` and ``without_stoats`` examine what happens when removing skuas (the lemming population explodes) and removing stoats (cycles of lower periodicity can remain, especially if increasing generalist predation rates). 
* ``sensitivity_analysis`` is another version of the code by JAH that allows to loop over parameter values. 
* ``ASS`` looks for alternative stable states that are frequent in such predation models (we did not find any so far, although the trajectories can be slightly different along the attractor, i.e., there is very likely a weak chaos). 
