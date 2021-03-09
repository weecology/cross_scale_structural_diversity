#the idea here is to train on NEON, test on FHM, adn apply on FIA
# features to derive: STATURE, delta H, BA, DBI, Density
NEON =  fread("indir/indir/neon_vst_data_022021.csv")
NEON  = NEON %>% select(individualID, eventID, siteID, domainID, plotID, subplotID, height, maxCrownDiameter, ninetyCrownDiameter, 
                        stemDiameter, taxonID, plantStatus, growthForm, canopyPosition, elevation, nlcdClass)

#since the other two datasets will be in imperial units, turn dbh, height, and diameter
NEON$stemDiameter = NEON$stemDiameter * 0.393701
NEON$height = NEON$height *3.280841666667
NEON$maxCrownDiameter = NEON$maxCrownDiameter *3.280841666667
NEON$ninetyCrownDiameter = NEON$ninetyCrownDiameter *3.280841666667
summary(NEON)

NEON = NEON %>% filter(growthForm %like% "tree")%>% 
  group_by(plotID, subplotID, eventID) %>% 
  mutate(DBI = n()* 10*(mean(stemDiameter^2, na.rm=T)/10)^1.605) %>%
  #number of trees per hectar
  mutate(DENSITY = n()*25) %>%
  mutate(BA = sum(10*0.005454*stemDiameter^2, na.rm = T))


NEON = NEON %>% filter((maxCrownDiameter > 5),#) %>%
                       (height > 2),
                       (stemDiameter > 0),
                       plantStatus %like% "Live") %>%
  group_by(individualID) %>% slice(1)
summary(NEON)

NEON = NEON %>% group_by(plotID, subplotID, eventID) %>% 
  mutate(STATURE = max(height, na.rm = T),
         H_VAR = sd(height, na.rm = T))

NEON = clean_typos_taxonID(NEON)
NEON = NEON %>% filter(taxonID %in% colnames(PHENOLOGY)) 
NEON = data.frame(NEON)

set.seed(1987)
# clean features that are too noisy
NEON$stemDiameter = sqrt(NEON$stemDiameter)
NEON$DBI = log(NEON$DBI)
NEON$DENSITY = log(NEON$DENSITY)
NEON$maxCrownDiameter = log(NEON$maxCrownDiameter)
NEON$ninetyCrownDiameter = log(NEON$ninetyCrownDiameter)
#scale all numeric variables
vars_to_scale = c("stemDiameter","height", "DBI", "DENSITY","BA",
                  "maxCrownDiameter", "ninetyCrownDiameter", 
                   "STATURE", "H_VAR")

NEON= data.frame(NEON)
NEON[is.na(NEON$H_VAR),"H_VAR"] = 0
# Save scaled attibutes:
scaleListNEON = scale(NEON[vars_to_scale])
scaleListNEON <- list(scale = attr(scaleListNEON, "scaled:scale"),
                  center = attr(scaleListNEON, "scaled:center"))

NEON[vars_to_scale] = scale(NEON[vars_to_scale])
train = NEON %>% group_by(siteID, taxonID) %>% sample_frac(0.8)
test =  NEON %>% filter(!individualID %in% train$individualID)


crown_hat.neon2 <-  brm(maxCrownDiameter ~ s(stemDiameter) + s(height)
                      # population level features
                      + s(H_VAR) 
                      + s(STATURE)
                      + s(DBI) + s(BA) + s(DENSITY)
                      #interaction terms
                      + height:stemDiameter
                      + DENSITY:stemDiameter + DENSITY:height
                      + BA:stemDiameter + BA:height
                      + DBI:stemDiameter + DBI:height
                      + H_VAR:stemDiameter + H_VAR:height
                      + STATURE:stemDiameter + STATURE:height
                      # scale nesting
                      + ((1 + stemDiameter + height 
                          + STATURE + BA
                      ) | gr(domainID, cor=F)) 
                      #taxon level features
                      + ((1 + (stemDiameter) + (height)
                          # population level features
                          + (H_VAR) 
                          + (STATURE)
                          + (DBI) + (BA) + (DENSITY)
                          #interaction terms
                          + height:stemDiameter
                          + DENSITY:stemDiameter + DENSITY:height
                          + BA:stemDiameter + BA:height
                          + DBI:stemDiameter + DBI:height
                          + H_VAR:stemDiameter + H_VAR:height
                          + STATURE:stemDiameter + STATURE:height
                      ) |gr(taxonID, cov = PHENOLOGY, cor=F))
                      # landuse level feature
                      + ((1 + stemDiameter + height 
                          + STATURE + BA
                      ) |gr(canopyPosition,cor=F))
                      ,
                      # define data and grouping correlation structures
                      data = train, family = student(),
                      data2 = list(PHENOLOGY = PHENOLOGY), 
                      #other parameters
                      chains=2, cores =3, backend = "cmdstanr", threads = threading(2))

bayes_R2(crown_hat.neon2, newdata = test, allow_new_levels = T)
test_predict = predict(crown_hat.neon, newdata = test, allow_new_levels=T)
test_predict = data.frame(test_predict)
predictions = cbind.data.frame(test_predict, test)

ggplot(predictions, aes( x = maxCrownDiameter, y = Estimate, color = taxonID)) + 
  geom_point(alpha = 0.5) + theme_bw()+ #xlim(0,100) + ylim(0,100)+
  #facet_wrap(siteID~.) +
  geom_abline(slope = 1)  + #facet_wrap(siteID~.) +
  theme(legend.position = "none")

write_rds(crown_hat.neon, "./out/models/crown_size_from_structure.rds")

#retrotransform scaled features
predictions$Estimate =  predictions$Estimate *
  scaleList$scale["maxCrownDiameter"] +
  scaleList$center["maxCrownDiameter"]
predictions$maxCrownDiameter =  predictions$maxCrownDiameter *
  scaleList$scale["maxCrownDiameter"] +
  scaleList$center["maxCrownDiameter"]
predictions$stemDiameter =  predictions$stemDiameter * 
  scaleList$scale["stemDiameter"] + 
  scaleList$center["stemDiameter"]
predictions$height =  predictions$height * 
  scaleList$scale["height"] + 
  scaleList$center["height"]

