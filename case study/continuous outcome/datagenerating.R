set.seed(2024)

###############################################################
#####                                                     #####                                                                       
#####                      MAJIC PV                       #####
#####                                                     ##### 
###############################################################

PV_treatment <- sample(rep(c(0,1),c(87,93)))

PV_age_0 <- round(rnorm(87,mean=66,sd=13))
PV_age_0[PV_age_0<28] <- 28
PV_age_0[which.min(PV_age_0)] <- 28
PV_age_0[PV_age_0>85] <- 85
PV_age_0[which.max(PV_age_0)] <- 85

PV_age_1 <- round(rnorm(93,mean=67,sd=13))
PV_age_1[PV_age_1<34] <- 34
PV_age_1[which.min(PV_age_1)] <- 34
PV_age_1[PV_age_1>88] <- 88
PV_age_1[which.max(PV_age_1)] <- 88

PV_sex_0 <- sample(rep(c(1,0),c(38,49)))
PV_sex_1 <- sample(rep(c(1,0),c(37,56)))

PV_haemoglobin_0 <- round(rnorm(87,mean=136,sd=17))
PV_haemoglobin_0[PV_haemoglobin_0<65] <- 65
PV_haemoglobin_0[which.min(PV_haemoglobin_0)] <- 65
PV_haemoglobin_0[PV_haemoglobin_0>163] <- 163
PV_haemoglobin_0[which.max(PV_haemoglobin_0)] <- 163


PV_haemoglobin_1 <- round(rnorm(93,mean=136,sd=17))
PV_haemoglobin_1[PV_haemoglobin_1<85] <- 85
PV_haemoglobin_1[which.min(PV_haemoglobin_1)] <- 85
PV_haemoglobin_1[PV_haemoglobin_1>173] <- 173
PV_haemoglobin_1[which.max(PV_haemoglobin_1)] <- 173

PV_therapy_0 <- sample(c(rep(1:2,c(40,30)),sample(3:6,17,replace = T)))
PV_therapy_1 <- sample(c(rep(1:2,c(55,30)),sample(3:4,8,replace = T)))

PV_thrombosis_0 <- sample(rep(c(0,1),c(49,38)))
PV_thrombosis_1 <- sample(rep(c(0,1),c(67,26)))

PV_intolerant_0 <- sample(rep(c(0,1),c(50,37)))
PV_intolerant_1 <- sample(rep(c(0,1),c(50,43)))

PV_intoresis_0 <- rep(0,87)
PV_intoresis_0[PV_intolerant_0==0] <- rep(c(0,1),c(23,27))
PV_intoresis_1 <- rep(0,93)
PV_intoresis_1[PV_intolerant_1==0] <- rep(c(0,1),c(31,19))

PV_splenomegaly_0 <- sample(rep(c(0,1),c(65,22)))
PV_splenomegaly_1 <- sample(rep(c(0,1),c(70,23)))

PV_splenectomy_0 <- rep(0,87)
PV_splenectomy_0[PV_splenomegaly_0==0] <- sample(rep(c(0,1),c(60,5)))
PV_splenectomy_1 <- rep(0,93)
PV_splenectomy_1[PV_splenomegaly_1==0] <- sample(rep(c(0,1),c(65,5)))

PV_WBC_0 <- rnorm(87,mean=9,sd=3)
PV_WBC_0[PV_WBC_0<2] <- 2
PV_WBC_0[which.min(PV_WBC_0)] <- 2
PV_WBC_0[PV_WBC_0>37] <- 37
PV_WBC_0[which.max(PV_WBC_0)] <- 37

PV_WBC_1 <- rnorm(93,mean=9,sd=4)
PV_WBC_1[PV_WBC_1<2] <- 2
PV_WBC_1[which.min(PV_WBC_1)] <- 2
PV_WBC_1[PV_WBC_1>73] <- 73
PV_WBC_1[which.max(PV_WBC_1)] <- 73


PV_platelets_0 <- rnorm(87,mean=356,sd=220)
PV_platelets_0[PV_platelets_0<99] <- 99
PV_platelets_0[which.min(PV_platelets_0)] <- 99
PV_platelets_0[PV_platelets_0>1420] <- 1420
PV_platelets_0[which.max(PV_platelets_0)] <- 1420

PV_platelets_1 <- rnorm(93,mean=401,sd=220)
PV_platelets_1[PV_platelets_1<61] <- 61
PV_platelets_1[which.min(PV_platelets_1)] <- 61
PV_platelets_1[PV_platelets_1>1546] <- 1546
PV_platelets_1[which.max(PV_platelets_1)] <- 1546

PV_JAK2V617F_0 <- sample(rep(c(0,1),c(2,85)))
PV_JAK2V617F_1 <- sample(rep(c(0,1),c(4,89)))

PV_0 <- data.frame(age=PV_age_0,
                   sex=PV_sex_0,
                   haemoglobin=PV_haemoglobin_0,
                   therapy=PV_therapy_0,
                   thrombosis=PV_thrombosis_0,
                   intolerant=PV_intolerant_0,
                   intoresis=PV_intoresis_0,
                   splemomegaly=PV_splenomegaly_0,
                   splenectomy=PV_splenectomy_0,
                   WBC=PV_WBC_0,
                   platelets=PV_platelets_0,
                   JAK2V617F=PV_JAK2V617F_0)

PV_1 <- data.frame(age=PV_age_1,
                   sex=PV_sex_1,
                   haemoglobin=PV_haemoglobin_1,
                   therapy=PV_therapy_1,
                   thrombosis=PV_thrombosis_1,
                   intolerant=PV_intolerant_1,
                   intoresis=PV_intoresis_1,
                   splemomegaly=PV_splenomegaly_1,
                   splenectomy=PV_splenectomy_1,
                   WBC=PV_WBC_1,
                   platelets=PV_platelets_1,
                   JAK2V617F=PV_JAK2V617F_1)

PV <- rbind(PV_0,PV_1)
PV[PV_treatment==1,] <- PV_1
PV[PV_treatment==0,] <- PV_0
PV$treatment <- PV_treatment

PV_centered <- data.frame(scale(PV,center = T,scale=F))

PV_intercept <- log(63/(180-63))

PV_score <- PV_intercept+PV_centered$treatment*log(2.03)+PV_centered$age*log(1.01)+PV_centered$haemoglobin*log(1.02)+
  PV_centered$therapy*log(0.77)+PV_centered$sex*log(0.97)+PV_centered$thrombosis*log(0.63)+PV_centered$intolerant*log(0.94)+
  PV_centered$intoresis*log(0.7)+PV_centered$splemomegaly*log(0.13)+PV_centered$splenectomy*log(1.26)

PV_error <- rnorm(n=180,mean=0,sd=sqrt(var(PV_score)/10))

PV_score <- PV_score + PV_error

###############################################################
#####                                                     #####                                                                       
#####                      MAJIC ET                       #####
#####                                                     ##### 
###############################################################

ET_treatment <- sample(rep(c(0,1),c(52,58)))

ET_age_0 <- round(rnorm(52,mean=65.6,sd=13.5))
ET_age_0[ET_age_0<37] <- 37
ET_age_0[which.min(ET_age_0)] <- 37
ET_age_0[ET_age_0>85] <- 85
ET_age_0[which.max(ET_age_0)] <- 85

ET_age_1 <- round(rnorm(58,mean=62.9,sd=12.3))
ET_age_1[ET_age_1<34] <- 34
ET_age_1[which.min(ET_age_1)] <- 34
ET_age_1[ET_age_1>90] <- 90
ET_age_1[which.max(ET_age_1)] <- 90

ET_sex_0 <- sample(rep(c(1,0),c(30,22)))
ET_sex_1 <- sample(rep(c(1,0),c(36,22)))

ET_haemoglobin_0 <- round(rnorm(52,mean=126,sd=17))
ET_haemoglobin_0[ET_haemoglobin_0<90] <- 90
ET_haemoglobin_0[which.min(ET_haemoglobin_0)] <- 90
ET_haemoglobin_0[ET_haemoglobin_0>160] <- 160
ET_haemoglobin_0[which.max(ET_haemoglobin_0)] <- 160


ET_haemoglobin_1 <- round(rnorm(58,mean=119,sd=17))
ET_haemoglobin_1[ET_haemoglobin_1<87] <- 87
ET_haemoglobin_1[which.min(ET_haemoglobin_1)] <- 87
ET_haemoglobin_1[ET_haemoglobin_1>152] <- 152
ET_haemoglobin_1[which.max(ET_haemoglobin_1)] <- 152

ET_therapy_0 <- sample(rep(1:6,c(15,20,8,5,2,2)))
ET_therapy_1 <- sample(rep(c(1:5,9),c(14,24,12,5,2,1)))


ET_intolerant_0 <- sample(rep(c(0,1),c(25,27)))
ET_intolerant_1 <- sample(rep(c(0,1),c(28,30)))

ET_splenomegaly_0 <- sample(rep(c(0,1),c(43,9)))
ET_splenomegaly_1 <- sample(rep(c(0,1),c(44,14)))

ET_splenectomy_0 <- rep(0,52)
ET_splenectomy_0[ET_splenomegaly_0==0] <- sample(rep(c(0,1),c(41,2)))
ET_splenectomy_1 <- rep(0,58)
ET_splenectomy_1[ET_splenomegaly_1==0] <- sample(rep(c(0,1),c(41,3)))

ET_WBC_0 <- rnorm(52,mean=6.8,sd=2.7)
ET_WBC_0[ET_WBC_0<2.8] <- 2.8
ET_WBC_0[which.min(ET_WBC_0)] <- 2.8
ET_WBC_0[ET_WBC_0>15.2] <- 15.2
ET_WBC_0[which.max(ET_WBC_0)] <- 15.2

ET_WBC_1 <- rnorm(58,mean=7.5,sd=4.8)
ET_WBC_1[ET_WBC_1<1.7] <- 1.7
ET_WBC_1[which.min(ET_WBC_1)] <- 1.7
ET_WBC_1[ET_WBC_1>29.8] <- 29.8
ET_WBC_1[which.max(ET_WBC_1)] <- 29.8


ET_platelets_0 <- rnorm(52,mean=573,sd=227.1)
ET_platelets_0[ET_platelets_0<166] <- 166
ET_platelets_0[which.min(ET_platelets_0)] <- 166
ET_platelets_0[ET_platelets_0>1406] <- 1406
ET_platelets_0[which.max(ET_platelets_0)] <- 1406

ET_platelets_1 <- rnorm(58,mean=545.4,sd=215.3)
ET_platelets_1[ET_platelets_1<89] <- 89
ET_platelets_1[which.min(ET_platelets_1)] <- 89
ET_platelets_1[ET_platelets_1>1139] <- 1139
ET_platelets_1[which.max(ET_platelets_1)] <- 1139

ET_JAK2V617F_0 <- sample(rep(c(0,1),c(26,26)))
ET_JAK2V617F_1 <- sample(rep(c(0,1),c(30,28)))

ET_CALR_0 <- rep(0,52)
ET_CALR_0[ET_JAK2V617F_0==0] <- sample(rep(c(0,1),c(12,14)))
ET_CALR_1 <- rep(0,58)
ET_CALR_1[ET_JAK2V617F_1==0] <- sample(rep(c(0,1),c(10,20)))


ET_0 <- data.frame(age=ET_age_0,
                   sex=ET_sex_0,
                   haemoglobin=ET_haemoglobin_0,
                   therapy=ET_therapy_0,
                   intolerant=ET_intolerant_0,
                   splemomegaly=ET_splenomegaly_0,
                   splenectomy=ET_splenectomy_0,
                   WBC=ET_WBC_0,
                   platelets=ET_platelets_0,
                   JAK2V617F=ET_JAK2V617F_0,
                   CALR=ET_CALR_0)

ET_1 <- data.frame(age=ET_age_1,
                   sex=ET_sex_1,
                   haemoglobin=ET_haemoglobin_1,
                   therapy=ET_therapy_1,
                   intolerant=ET_intolerant_1,
                   splemomegaly=ET_splenomegaly_1,
                   splenectomy=ET_splenectomy_1,
                   WBC=ET_WBC_1,
                   platelets=ET_platelets_1,
                   JAK2V617F=ET_JAK2V617F_1,
                   CALR=ET_CALR_1)

ET <- rbind(ET_0,ET_1)
ET[ET_treatment==1,] <- ET_1
ET[ET_treatment==0,] <- ET_0
ET$treatment <- ET_treatment

ET$platelets_median <- ET$platelets>quantile(ET$platelets,1/3) & ET$platelets<=quantile(ET$platelets,2/3)
ET$platelets_upper <- ET$platelets>quantile(ET$platelets,2/3)
ET$allele_neither <- 1-ET$JAK2V617F-ET$CALR

ET_centered <- data.frame(scale(ET,center = T,scale=F))

ET_intercept <- log(50/(110-50))

ET_score <- ET_intercept+ET_centered$treatment*log(1.14)+ET_centered$intolerant*log(1.17)+ET_centered$WBC*log(1.45)+
  ET_centered$platelets_median*log(0.58)+ET_centered$platelets_upper*log(0.44)+
  ET_centered$haemoglobin*log(1.82)+ET_centered$CALR*log(0.9)+ET_centered$allele_neither*log(0.79)

ET_error <- rnorm(n=110,mean=0,sd=sqrt(var(ET_score)/10))

ET_score <- ET_score + ET_error
