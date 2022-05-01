---
title: Simulating and Running an Economic Model in R
date: 2022-01-27 22:07:18
author: Helen Huier Ma
categories: Economics
tags:
 - R
 - Econometrics
---
#Load library
install.packages(sfcr)
library(sfcr)
library(tidyverse)

#Given below Canada's interprovincial trade model
rm(list = Is())
reg_eqs <- sfcr_set(
  
  Y_AB ~ C_AB + G_AB + EQZ_AB + X_AB - M_AB, 
  Y_ON ~ C_ON + G_ON + EQZ_ON + X_ON - M_ON, 
  Y_QC ~ C_QC + G_QC + EQZ_QC + X_QC - M_QC, 
  Y ~ Y_AB + Y_ON + Y_QC,
  
  M_AB_ON ~ mu_AB_ON * Y_AB, 
  M_AB_QC ~ mu_AB_QC * Y_AB, 
  M_AB ~ M_AB_ON + M_AB_QC,
  
  M_ON_AB ~ mu_ON_AB * Y_ON, 
  M_ON_QC ~ mu_ON_QC * Y_ON, 
  M_ON ~ M_ON_AB + M_ON_QC,
  
  M_QC_AB ~ mu_QC_AB * Y_QC, 
  M_QC_ON ~ mu_QC_ON * Y_QC, 
  M_QC ~ M_QC_ON + M_QC_AB,
  
  X_AB_ON ~ M_ON_AB, 
  X_AB_QC ~ M_QC_AB, 
  X_AB ~ X_AB_ON + X_AB_QC,
  
  X_ON_AB ~ M_AB_ON, 
  X_ON_QC ~ M_QC_ON, 
  X_ON ~ X_ON_AB + X_ON_QC,
  
  X_QC_AB ~ M_AB_QC, 
  X_QC_ON ~ M_ON_QC, 
  X_QC ~ X_QC_AB + X_QC_ON,
  
  YD_AB ~ Y_AB - TX_AB + r[-1] * Bh_AB[-1], 
  YD_ON ~ Y_ON - TX_ON + r[-1] * Bh_ON[-1], 
  YD_QC ~ Y_QC - TX_QC + r[-1] * Bh_QC[-1],
  
  TX_AB ~ theta * ( Y_AB + r[-1] * Bh_AB[-1] ), 
  TX_ON ~ theta * ( Y_ON + r[-1] * Bh_ON[-1] ), 
  TX_QC ~ theta * ( Y_QC + r[-1] * Bh_QC[-1] ),
  
  V_AB ~ V_AB[-1] + ( YD_AB - C_AB ), 
  V_ON ~ V_ON[-1] + ( YD_ON - C_ON ), 
  V_QC ~ V_QC[-1] + ( YD_QC - C_QC ),
  
  C_AB ~ alpha1_AB * YD_AB + alpha2_AB * V_AB[-1], 
  C_ON ~ alpha1_ON * YD_ON + alpha2_ON * V_ON[-1], 
  C_QC ~ alpha1_QC * YD_QC + alpha2_QC * V_QC[-1],
  
  Hh_AB ~ V_AB - Bh_AB, 
  Hh_ON ~ V_ON - Bh_ON, 
  Hh_QC ~ V_QC - Bh_QC,
  
  Bh_AB ~ V_AB * ( lambda0_AB + lambda1_AB * r - lambda2_AB * ( YD_AB/V_AB ) ), 
  Bh_ON ~ V_ON * ( lambda0_ON + lambda1_ON * r - lambda2_ON * ( YD_ON/V_ON ) ), 
  Bh_QC ~ V_QC * ( lambda0_QC + lambda1_QC * r - lambda2_QC * ( YD_QC/V_QC ) ),
  
  TX ~ TX_AB + TX_ON + TX_QC, 
  G ~ G_AB + G_ON + G_QC, 
  Bh ~ Bh_AB + Bh_ON + Bh_QC, 
  Hh ~ Hh_AB + Hh_ON + Hh_QC,
  
  Bs ~ Bs[-1] + ( G + EQZ_CA + r[-1] * Bs[-1] ) - ( TX + r[-1] * Bcb[-1] ), 
  Hs ~ Hs[-1] + Bcb - Bcb[-1], 
  Bcb ~ Bs - Bh,
  
  AVGY ~ Y/3,
  
  zeta_AB ~ if ( Y_AB[-1]-AVGY[-1] <0 ) {1} else {0}, 
  zeta_ON ~ if ( Y_ON[-1]-AVGY[-1] <0 ) {1} else {0}, 
  zeta_QC ~ if ( Y_QC[-1]-AVGY[-1] <0 ) {1} else {0},
  
  EQZ_AB ~ -zeta_AB*rho_eqz*(Y_AB[-1]-AVGY[-1]), 
  EQZ_ON ~ -zeta_ON*rho_eqz*(Y_ON[-1]-AVGY[-1]), 
  EQZ_QC ~ -zeta_QC*rho_eqz*(Y_QC[-1]-AVGY[-1]), 
  EQZ_CA ~ EQZ_AB + EQZ_ON + EQZ_QC,
  
  redondant ~ Hs - Hh
)

reg_ext <- sfcr_set( 
  r ~ 0.025,
  
  G_AB ~ 30, 
  G_ON ~ 30, 
  G_QC ~ 30, 
  rho_eqz ~ 0,
  
  mu_AB_ON ~ 0.1, 
  mu_AB_QC ~ 0.1,
  
  mu_ON_AB ~ 0.1, 
  mu_ON_QC ~ 0.1,
  
  mu_QC_AB ~ 0.1, 
  mu_QC_ON ~ 0.1,
  
  alpha1_AB ~ 0.70, 
  alpha1_ON ~ 0.70, 
  alpha1_QC ~ 0.70,
  
  alpha2_AB ~ 0.30, 
  alpha2_ON ~ 0.30, 
  alpha2_QC ~ 0.30,
  
  lambda0_AB ~ 0.7, 
  lambda0_ON ~ 0.7, 
  lambda0_QC ~ 0.7,
  
  lambda1_AB ~ 0.08, 
  lambda1_ON ~ 0.08, 
  lambda1_QC ~ 0.08,
  
  lambda2_AB ~ 0.01, 
  lambda2_ON ~ 0.01, 
  lambda2_QC ~ 0.01,
  
  theta ~ 0.25
)

#To simulate the above model 100 times

nsimul=100
reg=sfcr_baseline(
  equations=reg_eqs,
  external=reg_ext,
  periods=nsimul,
  method = "Broyden")
## above regression in short term: reg=sfcr_baseline(reg_eqs, reg_ext, nsimul)

#To find each the consumption level of ON, QC, and AB at stationary state
cat('C_AB at stationary state=',reg$C_AB[nsimul])
##C_AB at stationary state = 94.92669
cat('C_ON at stationary state=',reg$C_ON[nsimul])
##C_ON at stationary state = 94.92669
cat('C_QC at stationary state=',reg$C_QC[nsimul])
##C_QC at stationary state = 94.92669

#b)
#To find Canada's GDP at stationary state
cat('national GDP at stationary state=',reg$Y[nsimul])
##national GDP at stationary state = 374.7801

#c)
# Assume that the propensity (mu) of Quebec and Ontario to import goods and services by from Alberta increases from 0.1 to 0.125, as federal government spending in Ontario increases from 30 to 35, and that federal government spending in Quebec and Alberta drop from 30 to 27.5
##Create shock using start=3 and end=100
shkmu=sfcr_shock(
variables=sfcr_set(
  mu_ON_AB~0.125,
  mu_QC_AB~0.125,
  G_ON~35,
  G_AB~27.5,
  G_QC~27.5
  ),
  start = 3,
  end = 100
)
reg_shkmu=sfcr_scenario(baseline=reg,scenario=shkmu,periods = 100)

##Recall that we have already simulated the model 100 times at the beginning. Now calculate the % change between the provincial GDP at the initial stationary state and the new stationary state (after shock)
init_Y_AB=reg$Y_AB[nsimul]
init_Y_ON=reg$Y_ON[nsimul]
init_Y_QC=reg$Y_QC[nsimul]
new_Y_AB=reg_shkmu$Y_AB[nsimul]
new_Y_ON=reg_shkmu$Y_ON[nsimul]
new_Y_QC=reg_shkmu$Y_QC[nsimul]
cat('percent_change_AB=',(new_Y_AB-init_Y_AB)/init_Y_AB*100)
##percent_change_AB= 5.306341
cat('percent_change_ON=',(new_Y_ON-init_Y_ON)/init_Y_ON*100)
##percent_change_ON= 2.658365
cat('percent_change_QC=',(new_Y_QC-init_Y_QC)/init_Y_QC*100)
##percent_change_QC= -7.964688

#Based on your answer to c), which province saw its GDP increase the most compared to
#in the initial steady state? How do you explain this result, given the increase in spending by the
#federal government in Ontario and lower spending in Quebec and Alberta?

#d)Alberta sees its GDP increase the most. Given that Alberta and Quebec lower government spending by 2.5 while Ontario increase government spending by 5. Because both Ontario and Quebec are increasing their imports from Alberta simultaneously; in other words, Alberta has a higher exportation without increasing its imports from other provinces.So when other condition remains the same, the GDP of Ontario and Quebec tend to be lower According to this equation: Y ~ C + G + EQZ + X - M

#e)Calculate the percentage change between Canada's GDP at the initial steady state that you got in b) and the GDP at the new steady state. How do you explain this result?
cat('Y difference btw initial and new steady state=',reg_shkmu$Y[nsimul]-reg$Y[nsimul])
#Y difference btw initial and new steady state= 2.335871e-05
cat('Y percent change btw initial and new steady state=',(reg_shkmu$Y[nsimul]-reg$Y[nsimul])/reg$Y[nsimul]*100)
#Y percent change btw initial and new steady state= 6.232644e-06
#Canada's GDP in new steady state is slightly higher than in the initial state. Given the equation Y=Y_AB+Y_ON+Y_QC, and in c) we know that Y_AB increases 5.306341%, Y_ON increases 2.658365% while Y_QC decreases -7.964688%. The % increase in Y_AB and Y_ON almost overcome the % decrease in Y_QC. Therefore, Canada's GDP at initial steady state is very close to that in new steady state. The % change between them is minor. 

#f)The federal government's equalization program compensates provinces whose provincial GDP is lower than the average GDP of the Canadian provinces. The eqz (rho_eqz) parameter determines the generosity of the equalization program. So far, the basic calibration assumes that the program equalization does not compensate the poorer provinces, since eqz = 0.
#Given the new steady state obtained in c), identify the province (s) that would be eligible to receive an equalization payment!
cat('diff_AB_Average=',reg_shkmu$Y_AB[nsimul]-(reg_shkmu$Y[nsimul]/3))
#diff_AB_Average= 6.629029
cat('diff_ON_Average=',reg_shkmu$Y_ON[nsimul]-(reg_shkmu$Y[nsimul]/3))
#diff_ON_Average= 3.321
cat('diff_QC_Average=',reg_shkmu$Y_QC[nsimul]-(reg_shkmu$Y[nsimul]/3))
#diff_QC_Average= -9.950029
#Quebec would be eligible to receive an equalization payement, because its GDP is lower than the average GDP of the Canadaian provinces.

#g)Repeat the exercise in c), FURTHER assuming that the federal government decides that eqz = 0.5. What is the GDP of each province at the new steady state in this case, and identify the provinces that receive an equalization payment
nsimul=100
reg=sfcr_baseline(reg_eqs, reg_ext, nsimul)
shkeqz=sfcr_shock(
  variables=sfcr_set(
    mu_ON_AB~0.125,
    mu_QC_AB~0.125,
    G_ON~35,
    G_AB~27.5,
    G_QC~27.5,
    rho_eqz~0.5
  ),
  start = 3,
  end = 100
)
reg_shkeqz=sfcr_scenario(baseline=reg,scenario=shkeqz,periods = 100)

#le PIB de chaque province au nouvel état stationnaire
cat('Y_AB_eqz=',reg_shkeqz$Y_AB[nsimul])
##Y_AB_eqz= 134.4947
cat('Y_ON_eqz=',reg_shkeqz$Y_ON[nsimul])
##Y_ON_eqz= 130.5989
cat('Y_QC_eqz=',reg_shkeqz$Y_QC[nsimul])
##Y_QC_eqz= 122.9741

#Identify which province will get the equalized payment
cat('diff_AB_Average=',reg_shkeqz$Y_AB[nsimul]-(reg_shkeqz$Y[nsimul]/3))
##diff_AB_Average= 5.138824
cat('diff_ON_Average=',reg_shkeqz$Y_ON[nsimul]-(reg_shkeqz$Y[nsimul]/3))
##diff_ON_Average= 1.242993
cat('diff_QC_Average=',reg_shkeqz$Y_QC[nsimul]-(reg_shkeqz$Y[nsimul]/3))
##diff_QC_Average= -6.381817
##Quebec will receive an equalization payment in this case, because its provincal GDP is lower than Canada's average GDP.

#h) Compare the provincial GDP to the steady state obtained in g), with an equalization program, to the steady state provincial GDP obtained in c) without an equalization program. Identify the provinces whose GDP increases following the implementation of the equalization program?
cat('Y_AB difference between g) and c)=',reg_shkeqz$Y_AB[nsimul]-reg_shkmu$Y_AB[nsimul])
#Y_AB difference between g) and c)= 2.939014
cat('Y_ON difference between g) and c)=',reg_shkeqz$Y_ON[nsimul]-reg_shkmu$Y_ON[nsimul])
#Y_ON difference between g) and c)= 2.351211
cat('Y_QC difference between g) and c)=',reg_shkeqz$Y_QC[nsimul]-reg_shkmu$Y_QC[nsimul])
#Y_QC difference between g) and c)= 7.997429
#The GDP of all three provinces increases following the implementation of the equalization program, but Quebec's GDP increases the most.

#i) 
cat('equalization_payment_AB in h)=',reg_shkeqz$EQZ_AB[nsimul])
#equalization_payment_AB in h)= 0
cat('equalization_payment_ON in h)=',reg_shkeqz$EQZ_ON[nsimul])
#equalization_payment_ON in h)= 0
cat('equalization_payment_QC in h)=', reg_shkeqz$EQZ_QC[nsimul])
#equalization_payment_QC in h)= 3.190909
#It is possible that a province that does not receive an equalization payment still sees its GDP increase. For example, we know that in g) Alberta and Ontario will not receive an equlization payment because their provincial GDP is higher than Canada's average. However, in h), with an equalization program neither Alberta nor Ontario receive any equalization payment. Alberta's GDP increases 2.939014 and Ontario's GDP increases 2.351211. 

#j)Does the equalization program succeed in reducing the gaps between the GDP of each province and the average GDP of the provinces of Canada?
cat('diff_AB_Average=',reg_shkmu$Y_AB[nsimul]-(reg_shkmu$Y[nsimul]/3))
#diff_AB_Average= 6.629029
cat('diff_AB_Average_eqz=',reg_shkeqz$Y_AB[nsimul]-(reg_shkeqz$Y[nsimul]/3))
##diff_AB_Average_eqz= 5.138824

cat('diff_ON_Average=',reg_shkmu$Y_ON[nsimul]-(reg_shkmu$Y[nsimul]/3))
#diff_ON_Average= 3.321
cat('diff_ON_Average_eqz=',reg_shkeqz$Y_ON[nsimul]-(reg_shkeqz$Y[nsimul]/3))
##diff_ON_Average_eqz= 1.242993

cat('diff_QC_Average=',reg_shkmu$Y_QC[nsimul]-(reg_shkmu$Y[nsimul]/3))
#diff_QC_Average= -9.950029
cat('diff_QC_Average_eqz=',reg_shkeqz$Y_QC[nsimul]-(reg_shkeqz$Y[nsimul]/3))
##diff_QC_Average_eqz= -6.381817
#The equalization program succeed in reducing the gaps between the GDP of each province and the average GDP of the provinces of Canada.For example, with an equalization program, Quebec's GDP is 6.38 lower than the average. But it was 9.95 lower than the average when without the equalization payment.

Z<-matrix(c(400,0,100,100,275,125,250,180,525),3,3)
X<-matrix(c(950,600,1200),3,1)
950-(400+0+100)
#[1] 450
600-(100+275+125)
#[1] 100
1200-(250+180+525)
#[1] 245
f<-matrix(c(450,100,245),3,1)

#b)
400+100+250
#[1] 750
0+275+280
#[1] 555
100+125+525
#[1] 750
950-750
#[1] 200
600-555
#[1] 45
1200-750
#[1] 450
print(v<-matrix(c(200,45,450),1,3))
#[,1] [,2] [,3]
#[1,]  200   45  450