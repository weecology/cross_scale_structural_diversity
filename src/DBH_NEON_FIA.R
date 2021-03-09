library(data.table)
library(tidyverse)
library(brms)
library(GGally)

#get data from the FHM to get predictions of tree crown size for the FIA
NEON_FHM = fread("./indir/FHM_NEON.csv")
NEON_FIELD = fread("./indir/NEON_DF.csv") %>%
  filter(str_detect(plantStatus,"Live")) 
  
canopy_position = fread("./indir/indir/neon_vst_data_022021.csv", select=c("individualID", "canopyPosition")) %>%
  filter(!is.na(canopyPosition)) %>% 
  group_by(individualID) %>% slice(1)

axis1 = NEON_FIELD$right - NEON_FIELD$left
axis2= NEON_FIELD$top - NEON_FIELD$bottom
NEON_FIELD$CROWN_DIAMETER_90_hat = lapply(1:length(axis2), function(x)min(axis1[x], axis2[x])) %>%
  unlist
NEON_FIELD$CROWN_DIAMETER_WIDE_hat = lapply(1:length(axis2), function(x)max(axis1[x], axis2[x])) %>%
  unlist
NEON_FIELD = NEON_FIELD %>%
  select(individualID, domainID, siteID, taxonID, stemDiameter, height, 
         maxCrownDiameter, ninetyCrownDiameter,
         height.y, CROWN_DIAMETER_WIDE_hat,CROWN_DIAMETER_90_hat, 
         DENSITY, BA, DBI, plotID, nlcdClass)

NEON_FIELD = left_join(NEON_FIELD, canopy_position)

PHENOLOGY = fread("./indir/FHM_NEON_phylogeny.csv")
colnames(NEON_FIELD)[c(1,5:9, 15, 16)] = c("PK_TREE", "DBH", "TREE_HEIGHT", 
                                           "CROWN_DIAMETER_WIDE", "CROWN_DIAMETER_90", 
                                           "HEIGHT_HAT", "P3ID", "FOREST_TYPE_DESCR")
NEON_FIELD = NEON_FIELD %>% group_by(P3ID) %>% 
  mutate(STATURE = max(HEIGHT_HAT, na.rm = T))
NEON_FIELD = NEON_FIELD %>% group_by(P3ID) %>% 
  mutate(PLOT_MAX_WIDTH = max(CROWN_DIAMETER_WIDE_hat, na.rm = T))

NEON_FIELD = NEON_FIELD %>% group_by(P3ID) %>% 
  mutate(H_VAR = sd(HEIGHT_HAT, na.rm = T))

NEON_FIELD = clean_typos_taxonID(NEON_FIELD)
NEON_FIELD = NEON_FIELD %>% filter(taxonID %in% colnames(PHENOLOGY)) 
NEON_FIELD = data.frame(NEON_FIELD)

set.seed(1987)
# clean features that are too noisy
NEON_FIELD = NEON_FIELD %>%
  filter(DBH < 90) %>%
  filter(TREE_HEIGHT < 300) %>%
  filter(DENSITY < 2500) 

#check pairs of features
NEON_FIELD$CROWN_DIAMETER_WIDE_hat = NEON_FIELD$CROWN_DIAMETER_WIDE_hat* 3.28084
NEON_FIELD$CROWN_DIAMETER_90_hat = NEON_FIELD$CROWN_DIAMETER_90_hat* 3.28084

#NEON_FIELD = NEON_FIELD %>% filter(!canopyPosition %in% c("Full shade")) 

FHM = NEON_FIELD %>% 
  #filter(!is.na(canopyPosition)) %>%
  filter(!is.na(HEIGHT_HAT)) %>%
  filter(!is.na(DBH)) %>%
  filter(!is.na(CROWN_DIAMETER_WIDE_hat)) 

#FHM[FHM$canopyPosition != "Open grown", "canopyPosition"] = "Not open grown"

screen_features = FHM[c(5:14,17,19)]
(screen_features)
ggpairs(screen_features, aes(alpha = 0.4))

#transform all vars and standardize
# cor(sqrt(train$DBH), log(train$CROWN_DIAMETER_WIDE))
# cor((train$STAND_AGE), log(train$CROWN_DIAMETER_WIDE))
# cor(log(train$DBI), log(train$CROWN_DIAMETER_WIDE))
# cor((train$BA), log(train$CROWN_DIAMETER_WIDE))
# cor(sqrt(train$DENSITY), log(train$CROWN_DIAMETER_WIDE))

FHM$DBH = sqrt(FHM$DBH)
FHM$DBI = log(FHM$DBI)
FHM$DENSITY = log(FHM$DENSITY)
FHM$CROWN_DIAMETER_WIDE = log(FHM$CROWN_DIAMETER_WIDE)
FHM$CROWN_DIAMETER_90 = log(FHM$CROWN_DIAMETER_90)
FHM$CROWN_DIAMETER_WIDE_hat = log(FHM$CROWN_DIAMETER_WIDE_hat)
FHM$CROWN_DIAMETER_90_hat = log(FHM$CROWN_DIAMETER_90_hat)

#scale all numeric variables
vars_to_scale = c("DBH","TREE_HEIGHT", "DBI", "DENSITY","BA", "HEIGHT_HAT",
                  "CROWN_DIAMETER_WIDE", "CROWN_DIAMETER_90", 
                  "CROWN_DIAMETER_WIDE_hat", "CROWN_DIAMETER_90_hat", "STATURE", "PLOT_MAX_WIDTH", "H_VAR")

FHM= data.frame(FHM)
FHM[is.na(FHM$H_VAR),"H_VAR"] = 0
# Save scaled attibutes:
scaleList = scale(FHM[vars_to_scale])
scaleList <- list(scale = attr(scaleList, "scaled:scale"),
                  center = attr(scaleList, "scaled:center"))

FHM[vars_to_scale] = scale(FHM[vars_to_scale])

train_plts = FHM %>% group_by(P3ID) %>% sample_frac(0.1)
train =  FHM %>% filter(P3ID %in% train_plts$P3ID)
train = train %>% filter(taxonID %in% colnames(PHENOLOGY))
test =  FHM %>% filter(!P3ID %in% train$P3ID)

train = FHM %>% group_by(siteID, taxonID) %>% sample_frac(0.8)
test =  FHM %>% filter(!PK_TREE %in% train$PK_TREE)

# mv_structure <- brm(mvbind((stemDiameter), (maxCrownDiameter), (height)) ~ 
#                       # population level features
#                     + s(density) + s(basal_area)
#                     # scale nesting
#                     + (height| NA_L1CODE/ NA_L2CODE/ US_L3CODE/ US_L4CODE / UID) 
#                     #taxon level features
#                     + ((1 )|gr(taxonID, cov = PHENOLOGY)) #1 + height + stemDiameter + height:stemDiameter
#                     # landuse level feature
#                     + ((1 )|gr(nlcdClass)),
#                     # define data and grouping correlation structures
#                     data = train_, family = student(),
#                     data2 = list(PHENOLOGY = PHENOLOGY), 
#                     # define priors
#                     #prior = set_prior("horseshoe(3)", class = "b"),
#                     #c(prior(lkj(1), class = "rescor"),
#                     #prior(lkj(1), class = "cor"),
#                     #prior_string("normal(0,10)", resp = "stemDiameter"),
#                     #prior_string("normal(0,10)", resp = "ninetyCrownDiameter"),
#                     #prior_string("normal(0,10)",  resp = "maxCrownDiameter"),
#                     #prior_string("normal(0,10)",  resp = "height"),
#                     
#                     #set_recor=T,
#                     #other parameters
#                     chains=4, cores = 4)
# 
# bayes_R2(mv_structure)
# bayes_R2(crown_radius, newdata=test, allow_new_levels=T)

# crown_dbh.hat3 <- brm(DBH ~ s(CROWN_DIAMETER_WIDE_hat) + s(HEIGHT_HAT)
#                       # population level features
#                       + s(H_VAR) 
#                       #interaction terms
#                       + HEIGHT_HAT:CROWN_DIAMETER_WIDE_hat
#                       + H_VAR:CROWN_DIAMETER_WIDE_hat + H_VAR:HEIGHT_HAT
#                       # scale nesting
#                       + ((1 + HEIGHT_HAT:CROWN_DIAMETER_WIDE_hat) | domainID / siteID) 
#                       #taxon level features
#                       + ((1 + HEIGHT_HAT:CROWN_DIAMETER_WIDE_hat) |gr(taxonID, cov = PHENOLOGY)) #1 + height + stemDiameter + height:stemDiameter
#                       # landuse level feature
#                       #+ ((1 + HEIGHT_HAT:CROWN_DIAMETER_WIDE_hat) |gr(LAND_USE_CLASS_DESCR))
#                       + ((1 + HEIGHT_HAT:CROWN_DIAMETER_WIDE_hat) |gr(FOREST_TYPE_DESCR)),
#                       # define data and grouping correlation structures
#                       data = train, family = student(),
#                       data2 = list(PHENOLOGY = PHENOLOGY), 
#                       # define priors
#                       #prior = set_prior("horseshoe(1)"),
#                       #c(prior(lkj(1), class = "rescor"),
#                       #prior(lkj(1), class = "cor"),
#                       #prior_string("normal(0,10)", resp = "stemDiameter"),
#                       #prior_string("normal(0,10)", resp = "ninetyCrownDiameter"),
#                       #prior_string("normal(0,10)",  resp = "maxCrownDiameter"),
#                       #prior_string("normal(0,10)",  resp = "height"),
#                       #set_recor=T,
#                       #other parameters
#                       chains=4, cores = 4)

train = data.frame(train) 

crown_dbh.hat <-  brm(DBH ~ s(CROWN_DIAMETER_WIDE_hat) + s(HEIGHT_HAT)
                       # population level features
                       + s(H_VAR) 
                       + s(STATURE)
                       + s(PLOT_MAX_WIDTH)
                       #interaction terms
                       + HEIGHT_HAT:CROWN_DIAMETER_WIDE_hat
                       + DENSITY:CROWN_DIAMETER_WIDE_hat + DENSITY:HEIGHT_HAT
                       + H_VAR:CROWN_DIAMETER_WIDE_hat + H_VAR:HEIGHT_HAT
                       + PLOT_MAX_WIDTH:CROWN_DIAMETER_WIDE_hat + PLOT_MAX_WIDTH:HEIGHT_HAT
                       # scale nesting
                       + ((1 + CROWN_DIAMETER_WIDE_hat + HEIGHT_HAT 
                           + STATURE 
                           #+ PLOT_MAX_WIDTH 
                           #+ HEIGHT_HAT:CROWN_DIAMETER_WIDE_hat  
                           #+ STATURE:CROWN_DIAMETER_WIDE_hat 
                           #+ STATURE:HEIGHT_HAT
                       ) | domainID / siteID) 
                       #taxon level features
                       + ((1 + CROWN_DIAMETER_WIDE_hat 
                           + HEIGHT_HAT  
                           + STATURE
                           + PLOT_MAX_WIDTH
                           + H_VAR 
                           + HEIGHT_HAT:CROWN_DIAMETER_WIDE_hat  
                           + PLOT_MAX_WIDTH:HEIGHT_HAT  
                           + STATURE:HEIGHT_HAT
                           + H_VAR:HEIGHT_HAT
                       ) |gr(taxonID, cov = PHENOLOGY)) 
                       # landuse level feature
                       # + ((1 + CROWN_DIAMETER_WIDE_hat 
                       #     + HEIGHT_HAT 
                       #     #+ HEIGHT_HAT:CROWN_DIAMETER_WIDE_hat 
                       # ) |gr(canopyPosition))
                       ,
                       # define data and grouping correlation structures
                       data = train, family = asym_laplace(),
                       data2 = list(PHENOLOGY = PHENOLOGY), 
                       #other parameters
                       chains=2, cores =3, backend = "cmdstanr", threads = threading(2))

train$FOREST_TYPE_DESCR %>% table
test = test[complete.cases(test),]
test_ = test %>% filter(canopyPosition %in% c("Full sun", "Open grown", "Partially shaded"))
test_[test_$canopyPosition != "Open grown", "canopyPosition"] = "Not open grown"

#remove clear outliers
outliers= c("NEON.PLA.D12.YELL.00471", "9	NEON.PLA.D10.RMNP.03697", "NEON.PLA.D10.RMNP.03821")
test_ = test_ %>% filter(!PK_TREE %in% outliers)

#test_ = test %>% filter(canopyPosition %in% c("Full sun", "Open grown", "Partially shaded"))
test = test[complete.cases(test),]
test_ = test %>% filter(domainID %in% train$domainID) %>%
  filter(siteID %in% train$siteID) %>%
  filter(taxonID %in% train$taxonID) %>%
  filter(FOREST_TYPE_DESCR %in% train$FOREST_TYPE_DESCR)
bayes_R2(crown_dbh.hat.weibull, newdata = test, allow_new_levels = T)
write_rds(crown_dbh.hat.weibull, "./out/models/deepforest_student_family.rds")
test_known = test %>% filter(taxonID %in% unique(train$taxonID))
test_predict = predict(deepforest_student_family, newdata = test, allow_new_levels=T)
test_predict = data.frame(test_predict)
predictions = cbind.data.frame(test_predict, test)

#retrotransform scaled features
predictions$Estimate =  predictions$Estimate *
  scaleList$scale["DBH"] +
  scaleList$center["DBH"]
predictions$DBH =  predictions$DBH *
  scaleList$scale["DBH"] +
  scaleList$center["DBH"]
predictions$TREE_HEIGHT =  predictions$TREE_HEIGHT * 
  scaleList$scale["TREE_HEIGHT"] + 
  scaleList$center["TREE_HEIGHT"]
predictions$HEIGHT_HAT =  predictions$HEIGHT_HAT * 
  scaleList$scale["HEIGHT_HAT"] + 
  scaleList$center["HEIGHT_HAT"]
predictions$CROWN_DIAMETER_WIDE =  predictions$CROWN_DIAMETER_WIDE * 
  scaleList$scale["CROWN_DIAMETER_WIDE"] + 
  scaleList$center["CROWN_DIAMETER_WIDE"]

predictions$Q2.5 =  predictions$Q2.5 * 
  scaleList$scale["DBH"] + 
  scaleList$center["DBH"]

predictions$Q97.5 =  predictions$Q97.5 * 
  scaleList$scale["DBH"] + 
  scaleList$center["DBH"]

predictions$CROWN_DIAMETER_90 =  predictions$CROWN_DIAMETER_90 * 
  scaleList$scale["CROWN_DIAMETER_90"] + 
  scaleList$center["CROWN_DIAMETER_90"]
predictions$CROWN_DIAMETER_WIDE_hat =  predictions$CROWN_DIAMETER_WIDE_hat *
  scaleList$scale["CROWN_DIAMETER_WIDE_hat"] + 
  scaleList$center["CROWN_DIAMETER_WIDE_hat"]
predictions$CROWN_DIAMETER_90_hat =  predictions$CROWN_DIAMETER_90_hat *
  scaleList$scale["CROWN_DIAMETER_90_hat"] + 
  scaleList$center["CROWN_DIAMETER_90_hat"]
predictions$DENSITY =  predictions$DENSITY * 
  scaleList$scale["DENSITY"] + 
  scaleList$center["DENSITY"]

predictions$delta = abs(predictions$DBH^2 - predictions$Estimate^2)
predictions$TREE_HEIGHT =  predictions$TREE_HEIGHT/3.280841666667
pp = predictions %>% filter(siteID == "OSBS")
# ggplot(pp[-c(457, 615),], aes( x = DBH^2, y = Estimate^2)) + geom_point() + theme_bw()+
#   geom_abline(slope = 1) +facet_wrap(siteID~.)
ggplot(predictions, aes( x = DBH^2*2.54, y = Estimate^2*2.54)) + 
  geom_point(alpha = 0.5) + theme_bw()+ #xlim(0,100) + ylim(0,100)+
  #facet_wrap(siteID~.) +
  geom_abline(slope = 1)  + #facet_wrap(siteID~.) +
  theme(legend.position = "none")+
  geom_ribbon(aes(ymin=Q2.5^2*2.54, ymax=Q97.5^2*2.54), alpha=0.4) 

# conditions <- make_conditions(crown_dbh.hat5, "taxonID")
# ce = conditional_effects(crown_dbh.hat5, 
#                          effects = c("HEIGHT_HAT:CROWN_DIAMETER_WIDE_hat"),
#                          conditions = conditions)
# 
# ce


# evaluate if there is a pattern behind bad trees
offset = predictions %>% filter(delta > 5)
offset$Estimate = offset$Estimate^2 *2.54
offset$DBH = offset$DBH^2 *2.54

offset = offset %>% select(PK_TREE, taxonID, siteID, domainID, Estimate, DBH, HEIGHT_HAT, 
                           TREE_HEIGHT, CROWN_DIAMETER_WIDE_hat, CROWN_DIAMETER_WIDE,
                           FOREST_TYPE_DESCR, DENSITY, STATURE, PLOT_MAX_WIDTH, delta)
ggplot(offset, aes(y = delta*2.54, x = abs(CROWN_DIAMETER_WIDE_hat-CROWN_DIAMETER_WIDE), 
                   size = (abs(HEIGHT_HAT - TREE_HEIGHT)))) + theme_bw() + #geom_smooth(method="lm")+
  #scale_color_gradient2()+
  geom_point()+theme(legend.position = "bottom")#  + facet_wrap(siteID~.)
