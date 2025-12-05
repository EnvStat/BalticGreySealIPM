#Daire Carroll University of Gothenburg 2022
#Grey seal PVA script, estimated parameters for ppualiton growth in Baltic Grey Seals based on fitting

t = 100 #number of cycles representing 100 years

A = 47 #maximum age

logistic = function(y0,K,mu,x){
  (K * y0)/(y0 + (K - y0) * exp(-mu * x))
}

F26plus = 0.4166667
F = c(0,0,logistic(0.0004946038,0.8594168,1.572681,c(3:26)),rep((F26plus),(A)-26))#Each age class is assigned fertility of next cycle

K1 = 100000* 0.59
K2 = 88000* 0.59#carrying capacities
K3 = 120000 * 0.59   ##Added by Vanko et al.

P01 = 55997.36 * 0.59
P02 = 55736.80 * 0.59 #starting populations in 2020
P03 = 55997.36 * 0.59  ##Added by Vanko et al.

P2009_1 = 31163.51 * 0.59
P2009_2 = 30968.65 * 0.59 #starting populations in 2009
P2009_3 = 31422.64 * 0.59  ##Added by Vanko et al.

F1 = F*1.041543
F2 = F*1.053192 #estimated intrinsic fertilities
F3=F* 1.026829  ##Added by Vanko et al.

S1 = c(0.662       ,rep(0.932     , 3), rep(0.95     , (A-4)))
S2 = c(0.735       ,rep(0.932     , 3), rep(0.95     , (A-4))) #estimated survivals (disturbed population) FROM #4
S3 = c(0.5679371       ,rep(0.932     , 3), rep(0.95     , (A-4))) #estimated survivals (disturbed population)

H_sk = c(0.21329365, 0.07440476, 0.06448413, 0.05555556, 0.07142857) #hunting skew 
H_sk= c(0.2132936508, 0.0744047619, 0.0644841270 ,0.0555555556, 0.0714285714 ,0.0585317460, 0.0505952381, 0.0416666667 ,0.0406746032,
       0.0486111111, 0.0347222222, 0.0307539683, 0.0297619048, 0.0228174603 ,0.0178571429, 0.0297619048 ,0.0198412698 ,0.0178571429,
       0.0198412698, 0.0079365079 ,0.0089285714, 0.0099206349, 0.0079365079, 0.0049603175 ,0.0049603175, 0.0069444444 ,0.0019841270,
       0.0009920635, 0.0009920635, 0.0000000000 ,0.0009920635 ,0.0000000000 ,0.0000000000, 0.0000000000 ,0.0000000000 ,0.0000000000,
       0.0000000000 ,0.0000000000 ,0.0000000000, 0.0000000000 ,0.0000000000, 0.0000000000 ,0.0000000000, 0.0000000000, 0.0000000000,
       0.0000000000 ,0.0009920635)  ##Added by Vanko et al., c(age_count_h_pl$Proportion)from Script 5, needed for the msimulation to run



E_sk = c(0.53741497, 0.08163265, 0.11224490, 0.11904762, 0.04081633) #entanglement skew


R = 0.053 #random variation term

Sex_Ratio = 0.59 #Females as proportion of population
Hunt_Ratio = 0.43 #Females in the hunt

S_1_Hunt = c(0.709     , rep(0.952    , 3), rep(0.962   , (A-4)))
S_2_Hunt = c(0.784     , rep(0.953     , 3), rep(0.962    , (A-4)))#estimated survivals (undisturbed population)
S_3_Hunt = c(0.613     , rep(0.951     , 3), rep(0.961    , (A-4))) ##Added by Vanko et al.
