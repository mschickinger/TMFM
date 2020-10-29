This folder contains code for two important corrections of the TMFM data:
1) accounting for missed events due to limited temporal resolution as well as finite duration of the data acquisition.
2) identifying outlier data points that can arise from molecular heterogeneity.

The intended order of running scripts is:
run_postHMM.m
main_population_Corr.m

The concept is outlined in brief within the SI of the following publication:
M. Schickinger, M. Zacharias, and H. Dietz, “Tethered multifluorophore motion reveals equilibrium transition kinetics of single DNA double helices,” Proceedings of the National Academy of Sciences, vol. 115, no. 32, pp. E7512–E7521, 2018.

A more detailed description of the concept and its foundation is given in my Dissertation at the TUM:
Matthias Schickinger, “Untersuchung der Hybridisierungskinetik kurzer DNA-Doppelstraenge mit der neuartigen Einzelmolekuel-Methode Tethered Multi-Fluorophore Motion”, 2020.
(in German language)

I would like to point out that correction for missed events was implemented by myself in Matlab, but was adapted from the method published in:
J. Stigler and M. Rief, “Hidden Markov analysis of trajectories in single-molecule experiments and the effects of missed events,” ChemPhysChem, vol. 13, no. 4, pp. 1079–1086, 2012.