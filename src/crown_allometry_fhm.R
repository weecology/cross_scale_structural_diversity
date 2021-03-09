library(data.table)
library(tidyverse)
library(brms)
NEON_FHM = fread("./indir/FHM_NEON.csv")
#FHM = fread("./indir/FHM_STRUCTURE_PLOT.csv")
PHENOLOGY = fread("./indir/FHM_NEON_phylogeny.csv")

NEON_FHM$siteID %>% summary
NEON_FHM$siteID[is.na(NEON_FHM$siteID)] = "UNKN"

# NEON_FHM[NEON_FHM$dataset == "FHM", "stemDiameter"] = NEON_FHM[NEON_FHM$dataset == "FHM", "stemDiameter"] * 2.54
# NEON_FHM[NEON_FHM$dataset == "FHM", "ninetyCrownDiameter"] = NEON_FHM[NEON_FHM$dataset == "FHM", "ninetyCrownDiameter"] * 0.3048
# NEON_FHM[NEON_FHM$dataset == "FHM", "maxCrownDiameter"] = NEON_FHM[NEON_FHM$dataset == "FHM", "maxCrownDiameter"] * 0.3048
# NEON_FHM[NEON_FHM$dataset == "FHM", "height"] = NEON_FHM[NEON_FHM$dataset == "FHM", "height"] * 0.3048

NEON_FHM = NEON_FHM %>% filter((ninetyCrownDiameter)>0, (maxCrownDiameter)>0, 
                               (stemDiameter)>5, height > 2)

NEON_FHM = NEON_FHM %>% filter(taxonID %in% colnames(PHENOLOGY)) 
NEON_FHM = data.frame(NEON_FHM)
NEON_FHM["stemDiameter"] = sqrt(NEON_FHM["stemDiameter"])
is_numeric = c("stemDiameter","height","maxCrownDiameter", 
               "ninetyCrownDiameter", "DENSITY", "BASAL_AREA",
               "PLOT_SIZE", "ASPECT", "ELEVATION")
NEON_FHM = NEON_FHM %>% filter(ELEVATION < 10000)
NEON_FHM[NEON_FHM$dataset == "FHM", "ELEVATION"] = NEON_FHM[NEON_FHM$dataset == "FHM", "ELEVATION"] * 0.3048
# NEON_FHM[NEON_FHM$dataset == "FHM", "ASPECT"] = 
#   NEON_FHM[NEON_FHM$dataset == "FHM", "ASPECT"] * 0.3048
scaled = NEON_FHM
scaled[is_numeric] = scale(scaled[is_numeric])
set.seed(1987)


NEON_FHM = NEON_FHM %>% filter(dataset == "FHM")
train_plts = NEON_FHM %>% group_by(UID, taxonID) %>% sample_frac(0.2)

train =  NEON_FHM %>% filter(plotID %in% train_plts$plotID, subplotID %in% train_plts$subplotID,
                             (eventID %in% train_plts$eventID))
test =  NEON_FHM %>% filter(!UID %in% train$UID)

train_ = train %>% filter(taxonID %in% colnames(PHENOLOGY)) 

mv_structure <- brm(mvbind((stemDiameter), (maxCrownDiameter), (height)) ~ 
                      # population level features
                      + s(density) + s(basal_area)
                    # scale nesting
                    + (height| NA_L1CODE/ NA_L2CODE/ US_L3CODE/ US_L4CODE) 
                    #taxon level features
                    + ((1 )|gr(taxonID, cov = PHENOLOGY)) #1 + height + stemDiameter + height:stemDiameter
                    # landuse level feature
                    + ((1 )|gr(nlcdClass)),
                    # define data and grouping correlation structures
                    data = train_, family = student(),
                    data2 = list(PHENOLOGY = PHENOLOGY), 
                    # define priors
                    #prior = set_prior("horseshoe(3)", class = "b"),
                    #c(prior(lkj(1), class = "rescor"),
                    #prior(lkj(1), class = "cor"),
                    #prior_string("normal(0,10)", resp = "stemDiameter"),
                    #prior_string("normal(0,10)", resp = "ninetyCrownDiameter"),
                    #prior_string("normal(0,10)",  resp = "maxCrownDiameter"),
                    #prior_string("normal(0,10)",  resp = "height"),
                    
                    #set_recor=T,
                    #other parameters
                    chains=4, cores = 4)

bayes_R2(mv_structure)
bayes_R2(mv_structure, newdata=test, allow_new_levels=T)

write_rds(mv_structure, "./out/models/mv_structure.rds")

# plot predictions
test_predict = predict(crown_radius.fhm, newdata = test, allow_new_levels=T)
test_predict = data.frame(test_predict)
predictions = cbind.data.frame(test_predict, test)
ggplot(predictions, aes( x = maxCrownDiameter, y = Estimate)) + geom_point() + theme_bw()+
  geom_abline(slope = 1) +facet_wrap(NA_L2NAME~.)


library(bayesplot)
color_scheme_set("red")
posterior <- as.array(crown_radius.fhm)
summary(posterior)
dimnames(posterior)
boundaries = list()
for(ii in 1:dim(posterior)[3]){
  tmp = posterior[,,ii]
  tmp = unlist(tmp)
  lb = quantile(tmp, 0.20)
  cnt = median(tmp)
  ub = quantile(tmp, 0.80)
  boundaries[[ii]] = c(lb, cnt, ub)
}
boundaries = do.call(rbind.data.frame, boundaries)
colnames(boundaries) = c("Q20", "Q50", "Q80")
cross_zero = apply(boundaries, 1, function(x) x["Q20"]*x["Q80"] > 0)
#cross_zero = apply(boundaries, 1, function(x)abs(x[2])> 0.1)
cross_zero[length(cross_zero)] = F
dimnames(posterior)$parameters[cross_zero]
good_preds = cbind.data.frame(dimnames(posterior)$parameters[cross_zero],
                              boundaries[cross_zero,])
plt = mcmc_intervals(posterior[,,cross_zero])
plt



crown_radius <- brm(maxCrownDiameter ~ 
                      # population level features
                      s(stemDiameter) + s(height) + s(ELEVATION) #+ s(ASPECT)
                    + s(DENSITY) #+ s(BASAL_AREA)  
                    #interaction terms
                    + height:stemDiameter
                    # + BASAL_AREA:stemDiameter + BASAL_AREA:height 
                    + DENSITY:stemDiameter + DENSITY:height 
                    + ELEVATION:stemDiameter + ELEVATION:height 
                    #+ ASPECT:stemDiameter + ASPECT:height 
                    # scale nesting
                    + (1+height:stemDiameter | NA_L1CODE/ NA_L2CODE/ US_L3CODE/ US_L4CODE / PID) 
                    #taxon level features
                    + ((1+height:stemDiameter)|gr(taxonID, cov = PHENOLOGY)) #1 + height + stemDiameter + height:stemDiameter
                    # landuse level feature
                    + ((1+height:stemDiameter)|gr(nlcdClass)) 
                    + ((1+height:stemDiameter)|gr(BIO_PLOT_SOIL_DEPTH)), 
                    #+ ((1+height:stemDiameter)|gr(BIO_PLOT_DISTURBANCE)),
                    # define data and grouping correlation structures
                    data = train_, family = student(),
                    data2 = list(PHENOLOGY = PHENOLOGY), 
                    # define priors
                    #other parameters
                    chains=4, cores = 4)

crown_radius.u1

bayes_R2(crown_radius)
bayes_R2(crown_radius, newdata=test, allow_new_levels=T)

write_rds(crown_radius, "./out/models/fhm_crown_radius.rds")

# plot predictions
test_predict = predict(crown_radius, newdata = test, allow_new_levels=T)
test_predict = data.frame(test_predict)
predictions = cbind.data.frame(test_predict, test)
ggplot(predictions, aes( x = maxCrownDiameter, y = Estimate)) + geom_point() + theme_bw()+
  geom_abline(slope = 1) +facet_wrap(NA_L2NAME~.)

predictions$delta = abs(predictions$maxCrownDiameter - predictions$Estimate)

very_off = predictions$delta > 4
bayes_R2(crown_radius, newdata=test[!very_off,], allow_new_levels=T)

ggplot(predictions[!very_off,], aes( x = maxCrownDiameter, y = Estimate)) + geom_point() + theme_bw()+
  geom_abline(slope = 1) +facet_wrap(NA_L2NAME~.)

cor(predictions$Estimate, predictions$maxCrownDiameter)^2
