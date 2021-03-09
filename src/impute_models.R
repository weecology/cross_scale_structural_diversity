library(data.table)
library(tidyverse)
library(brms)
library(GGally)

#get data from the FHM to get predictions of tree crown size for the FIA
NEON_FHM = fread("./indir/FHM_NEON.csv")

PHENOLOGY = fread("./indir/FHM_NEON_phylogeny.csv")
colnames(NEON_FHM)
NEON_FHM = NEON_FHM %>% filter(CROWN_DIAMETER_WIDE>6.5, 
                               DBH>5, TREE_HEIGHT > 2)

NEON_FHM = NEON_FHM %>% filter(taxonID %in% colnames(PHENOLOGY)) 
NEON_FHM = data.frame(NEON_FHM)
set.seed(1987)
# clean features that are too noisy
NEON_FHM = NEON_FHM %>%
  filter(DBH < 90) %>%
  filter(TREE_HEIGHT < 300) %>%
  filter(DBI < 120000) %>%
  filter(BA < 600) %>%
  filter(DENSITY < 2500)


# clean features that are too noisy
NEON_FHM$DBH = sqrt(NEON_FHM$DBH)
NEON_FHM$DBI = log(NEON_FHM$DBI)
NEON_FHM$DENSITY = log(NEON_FHM$DENSITY)
NEON_FHM$CROWN_DIAMETER_WIDE = log(NEON_FHM$CROWN_DIAMETER_WIDE)
#NEON_FHM$ninetyCrownDiameter = log(NEON_FHM$ninetyCrownDiameter)

#calculate stature and height variance
NEON_FHM = NEON_FHM %>% group_by(P3ID, POINT_NBR, YEAR) %>% 
  mutate(STATURE = max(TREE_HEIGHT, na.rm = T),
         H_VAR = sd(TREE_HEIGHT, na.rm = T))

#correct H_VAR to 0 for plots with only 1 tree height
NEON_FHM = data.frame(NEON_FHM)
#scale all numeric variables
vars_to_scale = c("DBH","TREE_HEIGHT", "DBI", "DENSITY","BA",
                  "CROWN_DIAMETER_WIDE", "STATURE", "H_VAR")
FHM = NEON_FHM %>% filter(dataset == "FHM")
NEON = NEON_FHM %>% filter(dataset == "NEON")

climate_neon = fread("./indir/climate/NEON.csv") %>%select(-one_of("geometry"))
NEON = left_join(NEON, climate_neon, by = c("TREE_NBR" = "individualID"))
NEONENV = NEON %>% filter(!is.na(ppt_11))
ENV = prcomp(NEONENV[colnames(climate_neon)[-1]])[1:10]
climate_fhm = fread("./indir/climate/FHM.csv")%>%select(-one_of("geometry"))
FHM = left_join(FHM, climate_fhm, by = "TREE_NBR")
FHMENV = FHM %>% filter(!is.na(ppt_11))

ENV.fhm = predict(ENV)
NEONENV = cbind.data.frame(NEONENV[,1:49], data.frame(ENV["x"]))
                           
#NEON_FHM_[NEON_FHM_$dataset=="NEON", "TREE_HEIGHT"] = NEON_FHM_[NEON_FHM_$dataset=="NEON", "TREE_HEIGHT"] *2. 
#check which var is not overlapping
ggplot(NEON_FHM, aes(x=(CROWN_DIAMETER_WIDE), fill = dataset))+geom_histogram(alpha = 0.5)
vars_to_scale = c(vars_to_scale, paste("x.PC", 1:10, sep=""))
scaleList = scale(NEONENV[vars_to_scale])
scaleList <- list(scale = attr(scaleList, "scaled:scale"),
                      center = attr(scaleList, "scaled:center"))

NEON[vars_to_scale] = scale(NEON[vars_to_scale])
# scale FHM in NEON space
FHM[vars_to_scale] = scale(FHM[vars_to_scale], center = scaleList$center, 
                           scale = scaleList$scale)
NEON$siteID = factor(NEON$siteID)
NEON$taxonID = factor(NEON$taxonID)
NEON["domainID"] = substr(NEON$PK_TREE, 10,12)
NEON= NEON %>% filter(!is.na(NA_L1CODE))
train = NEON %>% group_by(taxonID) %>% dplyr::sample_frac(0.8)
test =  NEON %>% data.frame %>% filter(!TREE_NBR %in% train$TREE_NBR)

crown_hat.neon <-  brm(CROWN_DIAMETER_WIDE ~ s(DBH) + s(TREE_HEIGHT)
                       # population level features
                       + s(x.PC1)+ s(x.PC2)+ s(x.PC3)+ s(x.PC4)+ s(x.PC5)
                       + s(x.PC6)+ s(x.PC7)+ s(x.PC8)+ s(x.PC9)+ s(x.PC10)
                       + s(H_VAR) 
                       + s(STATURE)
                       + s(DBI) + s(BA) + s(DENSITY)
                       #interaction terms
                       + TREE_HEIGHT:DBH
                       + DENSITY:DBH + DENSITY:TREE_HEIGHT
                       + BA:DBH + BA:TREE_HEIGHT
                       + DBI:DBH + DBI:TREE_HEIGHT
                       # scale nesting
                       + ((1 + DBH + TREE_HEIGHT 
                           #+ BA + DBI + DENSITY
                       ) | NA_L1NAME / domainID ) 
                       #taxon level features
                       + ((1 + (DBH) + (TREE_HEIGHT)
                           + (x.PC1)+ (x.PC2)+ (x.PC3)+ (x.PC4)+ (x.PC5)
                           + (x.PC6)+ (x.PC7)+ (x.PC8)+ (x.PC9)+ (x.PC10)
                           # population level features
                           + (H_VAR) 
                           + (STATURE)
                           + (DBI) + (BA) + (DENSITY)
                           #interaction terms
                           + TREE_HEIGHT:DBH
                           + DENSITY:DBH + DENSITY:TREE_HEIGHT
                           + BA:DBH + BA:TREE_HEIGHT
                           + DBI:DBH + DBI:TREE_HEIGHT
                       ) |gr(taxonID, cov = PHENOLOGY, cor=F))
                       # landuse level feature
                       + ((1 + DBH + TREE_HEIGHT
                          #+ STATURE + BA + DBI
                       ) |gr(canopyPosition,cor=F))
                       ,
                       # define data and grouping correlation structures
                       data = (train), family = student(),
                       data2 = list(PHENOLOGY = PHENOLOGY), 
                       #other parameters
                       chains=2, cores =6, backend = "cmdstanr", threads = threading(3))

test = test %>% select(H_VAR, CROWN_DIAMETER_WIDE, DBH, TREE_HEIGHT, STATURE, 
                       x.PC1, x.PC2, x.PC3, x.PC4, x.PC5, x.PC6, x.PC7, x.PC8, x.PC9, x.PC10,
                       DBI, BA, DENSITY, taxonID, canopyPosition, NA_L1NAME, domainID)
test = test_fhm %>% filter(!is.na(H_VAR))
test = test[complete.cases(test),]
#test = test %>% filter(CROWN_DIAMETER_WIDE > -2)
bayes_R2(crown_hat.neon, newdata = test, allow_new_levels =T)
write_rds(crown_hat.neon, "./out/models/crown_hat.neon_env.rds")

# append domainID to FHM data points
domains = fread("./indir/ecoregions_domain.csv")
test_fhm = FHM %>% select(STATE, COUNTY, H_VAR, CROWN_DIAMETER_WIDE, DBH, TREE_HEIGHT, STATURE, 
                               DBI, BA, DENSITY, taxonID, canopyPosition, NA_L1NAME)
test_fhm = inner_join(test_fhm, domains) 
test_fhm$canopyPosition = "Partially shaded"
test_fhm = test_fhm[complete.cases(test_fhm),]

#test NEON model on NEON test and FHM
bayes_R2(crown_hat.neon, newdata = test_fhm, allow_new_levels = T)
bayes_R2(crown_hat.neon, newdata = test, allow_new_levels = T)


test_predict = predict(crown_hat.neon, newdata = test, allow_new_levels = T)
test_predict = data.frame(test_predict)
predictions = cbind.data.frame(test_predict, test)
columns_to_rescale = c("Estimate", "Est.Error","Q2.5","Q97.5", "CROWN_DIAMETER_WIDE")
predictions[,columns_to_rescale] = predictions[,columns_to_rescale] * scaleList$scale["CROWN_DIAMETER_WIDE"] + 
  scaleList$center["CROWN_DIAMETER_WIDE"]
ggplot(predictions, aes( x = exp(CROWN_DIAMETER_WIDE)*0.3048, y = exp(Estimate)*0.3048)) + geom_point() + theme_bw()+
  geom_abline(slope = 1) +# ylim(0,40)+#facet_wrap(NA_L2NAME~., scales = "free")+
  geom_ribbon(aes(ymin=exp(Q2.5)*0.3048, ymax=exp(Q97.5)*0.3048), alpha=0.4) 

predictions$delta =  exp(predictions$CROWN_DIAMETER_WIDE) -  exp(predictions$Estimate) 
summary(predictions$delta*0.3048)
quantile(predictions$delta*0.3048, 0.025)
quantile(predictions$delta*0.3048, 0.975)
offset = predictions %>% filter(abs(delta*0.3048) > 4)
