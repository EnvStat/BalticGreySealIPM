-comparison.R: compare the population trajectory and growth rate of Caroll's and our model
-Caroll_original: unchanged code of Caroll et al
-Caroll_modified: Our modification to run the compariosn

Changes to Caroll et al.'s code:

1_Ordered_scripts.R,2_Historical_data.R, Generation_time_calculator.R: NO modification (except name of parent folder)
3_Explicit_population_modelling.R,4_Fertility_survival_estimation.R,5_Undisturbed_population_modelling_MOD.R: Added third case for K= 120000


Baisc_model.R:  added E, E_sk as variables of the PopMod() function, otherwise it was not running (the function does not end up using it though)

Parameters.R: 
	- Added 3rd scenario for K=120,000
	- Note: values evalutatd earlier are copied here. Reducing the population/ K by *0.59 comes from only modelling females, who are assumed to be 59% of the population
	- F1,F2,F3: fertility corresponding to the 3 carrying capacity scenarios. F is multiplied by the effect of  density dependence, which is copied in here, but actually computed in Script 3 as Q_m_1, Q_m_2, Q_m_3 (So Fi=F*Q_m_i for i=1,2,3)
	- Added the value for H_sk, originally H_sk here was declared as a vector of length 5, which caused the simulation to not run, as it needs to be a vector of length 46. H_sk is originally calculated from Script 5, with H_sk=c(age_count_h_pl$Proportion)

6_Scenarios_MOD.R: Added 3rd scenario (K=120,000) & Reduced to the cases which were presented by Caroll et al., as we only modelled these, i.e. fertility not reduced, no good year-bad year scenarios, hunting bias towards females is 0.43 (historical), hunting bias hunting pressure is 0, 2400 or 3600

7_Risk_averages.R: 
	- Added 3rd scenario
	- Removed unused scenarios
	- Added simulation with our model
	- Added comparison plot

