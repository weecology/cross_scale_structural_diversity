
#transform all vars and standardize
cor(sqrt(train$DBH), log(train$CROWN_DIAMETER_WIDE))
cor((train$STAND_AGE), log(train$CROWN_DIAMETER_WIDE))
cor(log(train$DBI), log(train$CROWN_DIAMETER_WIDE))
cor((train$BA), log(train$CROWN_DIAMETER_WIDE))
cor(sqrt(train$DENSITY), log(train$CROWN_DIAMETER_WIDE))

FHM$DBH = sqrt(FHM$DBH)
FHM$DBI = log(FHM$DBI)
FHM$DENSITY = log(FHM$DENSITY)
FHM$CROWN_DIAMETER_WIDE = log(FHM$CROWN_DIAMETER_WIDE)
FHM$CROWN_DIAMETER_90 = log(FHM$CROWN_DIAMETER_90)

#scale all numeric variables
vars_to_scale = c("DBH", "TREE_HEIGHT", "STAND_AGE", "CROWN_RATIO" , "DBI", "DENSITY","BA", 
                  "CROWN_DIAMETER_WIDE", "CROWN_DIAMETER_90")

FHM[vars_to_scale] = scale(FHM[vars_to_scale])

train_plts = FHM %>% group_by(P3ID) %>% sample_frac(0.15)
train =  FHM %>% filter(P3ID %in% train_plts$P3ID)
test =  FHM %>% filter(!UID %in% train$UID)


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

crown_radius <- brm(CROWN_DIAMETER_WIDE ~ s(DBH) + s(TREE_HEIGHT)
                    # population level features
                    + s(BA) + s(DBI) + s(DENSITY) + s(STAND_AGE)
                    #interaction terms
                    + TREE_HEIGHT:DBH 
                    + BA:DBH + BA:TREE_HEIGHT 
                    + DBI:DBH + DBI:TREE_HEIGHT 
                    + DENSITY:DBH + DENSITY:TREE_HEIGHT
                    + STAND_AGE:DBH + STAND_AGE:TREE_HEIGHT
                    # scale nesting
                    + ((1 + TREE_HEIGHT:DBH) | NA_L1CODE/ NA_L2CODE/ US_L3CODE/ US_L4CODE / P3ID) 
                    #taxon level features
                    + ((1 + TREE_HEIGHT:DBH) |gr(taxonID, cov = PHENOLOGY)) #1 + height + stemDiameter + height:stemDiameter
                    # landuse level feature
                    + ((1 + TREE_HEIGHT:DBH) |gr(LAND_USE_CLASS_DESCR))
                    + ((1 + TREE_HEIGHT:DBH) |gr(FOREST_TYPE_DESCR)),
                    # define data and grouping correlation structures
                    data = train, family = student(),
                    data2 = list(PHENOLOGY = PHENOLOGY), 
                    # define priors
                    #prior = set_prior("horseshoe(1)"),
                    #c(prior(lkj(1), class = "rescor"),
                    #prior(lkj(1), class = "cor"),
                    #prior_string("normal(0,10)", resp = "stemDiameter"),
                    #prior_string("normal(0,10)", resp = "ninetyCrownDiameter"),
                    #prior_string("normal(0,10)",  resp = "maxCrownDiameter"),
                    #prior_string("normal(0,10)",  resp = "height"),
                    #set_recor=T,
                    #other parameters
                    chains=4, cores = 4)

crown_radius.all <- brm((maxCrownDiameter) ~ 
                          # population level features
                          s(BASAL_AREA) + s(sqrt(stemDiameter)) + s(height) +s(ELEVATION)
                        #interaction terms
                        + sqrt(stemDiameter):BASAL_AREA + BASAL_AREA:height + sqrt(stemDiameter):height
                        # scale nesting
                        + ((1 +(ELEVATION) +  sqrt(stemDiameter):height) | dataset) 
                        #taxon level features
                        + ((1+sqrt(stemDiameter):height)|gr(taxonID, cov = PHENOLOGY)) #1 + height + stemDiameter + height:stemDiameter
                        # landuse level feature
                        + ((1+(ELEVATION)+sqrt(stemDiameter):height)|gr(nlcdClass)),
                        # define data and grouping correlation structures
                        data = train_, family = student(),
                        data2 = list(PHENOLOGY = PHENOLOGY), 
                        # define priors
                        #prior = set_prior("horseshoe(1)"),
                        #c(prior(lkj(1), class = "rescor"),
                        #prior(lkj(1), class = "cor"),
                        #prior_string("normal(0,10)", resp = "stemDiameter"),
                        #prior_string("normal(0,10)", resp = "ninetyCrownDiameter"),
                        #prior_string("normal(0,10)",  resp = "maxCrownDiameter"),
                        #prior_string("normal(0,10)",  resp = "height"),
                        
                        #set_recor=T,
                        #other parameters
                        chains=4, cores = 4)

train$FOREST_TYPE_DESCR %>% table

bayes_R2(crown_radius, newdata = test, allow_new_levels = T)
bayes_R2(crown_radius_fhm, newdata = test, allow_new_levels = T)
bayes_R2(fhm_crown_radius, newdata = test, allow_new_levels = T)

write_rds(crown_radius, "./out/models/crown_radius_fhm.rds")
test$canopyPosition = "Partially shaded"
test$H_VAR = 0
colnames(test)[3:4] =c("stemDiameter", "height" )
test_ = test_neon %>% filter(CROWN_DIAMETER_WIDE > -2)
test_ = data.frame(test_)
test_ = test_ %>% select(colnames(test_fhm))
test_ = test_[complete.cases(test_),]