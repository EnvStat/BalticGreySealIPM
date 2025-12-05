Data from Source_data was processed in R and prepared in the format where the "fit_run.R" R script can read it in to fit the stan model. R scripts used for processing can be found in https://github.com/jpvanhat/milena/tree/main/data_processing


"Aerial_count.csv"
	- Year: 2003-2024
	- Total_count: Total yearly aerial survey count
	- Processed with "aerial_count.R"

"Bycatch_samples_FI.csv"
	- Number of bycaught seals in the sample from each demographic group in each year
	- 2005-2024
	- Processed with "samples_FI.R"


"Bycatch_samples_SW.csv"
	- Number of bycaught seals in the sample from each demographic group in each year
	- 2002-2024
	- Processed with "samples_SW.R"

"Herring.csv"
	- 2001 - 2024
	- BP_GoF, GoB: mean weight of herring over age 5 in the the Baltic Proper+Gulf of Finland (BP_GoF) and in the Gulf of Bothnia (GoB)  
	- - BP_GoF.norm, GoB.norm: normalized mean weight of herring over age 5 in the the Baltic Proper+Gulf of Finland (BP_GoF) and in the Gulf of Bothnia (GoB) 
	
"Hunting_bag_quota_FI_SW.csv"
	- 2002-2024
	- Hunting_bag_FI: Reported Finish hunting bag 
	- Hunting_bag_SW: Reported Swedish hunting bag 
	- Quota_FI: Yearly quota in Finland
	- Quota_SW: Yearly quota in Sweden
	- Processe with "hunting_bag.R"
	
"Hunting_samples_FI.csv"
	- Number of hunted seals in the sample from each demographic group in each year
	- 2002-2024
	- Processed with "samples_FI.R"

"Hunting_samples_SW.csv"
	- Number of hunted seals in the sample from each demographic group in each year
	- 2002-2024
	- Processed with "samples_SW.R"

"Pregnancy_samples_SW.csv"
	- Year:2002-2024
	- Pregnant: Number of pregnant seals out of all the sampled seals Sweden that could be pregnant (i.e. adult females sampled after July, when the implantation of the embiro happens)
	- Total: Number of sampled seals in Sweden that could be pregnant 
	- Processed with "samples_SW.R"
	
"Reproductive_signs_FI.csv"
	- Prsence/abscense data of reproductive signs (placental scar and CA) on sampled seals in Finland that coul have reproductive signs (i.e. adult females that died before August, when the scars haven't yet faded)
	- Year: 2002-2024
	- 1: no placental scar or CA present
	- 2: scar present, CA not
	- 3: CA present, scar not
	- 4: both scar and CA present
	
	
