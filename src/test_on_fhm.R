test_fhm_ = test_fhm_[complete.cases(test_fhm_),]
test_fhm_$STATURE = test_fhm_$height
test_fhm_$H_VAR = mean(train$H_VAR)
test_fhm_$canopyPosition = "Partially shaded"

bayes_R2(crown_hat.neon, newdata = data.frame(test_fhm_), allow_new_levels = T)
train$canopyPosition %>% table
test_fhm$nlcdClass %>% table



# get NEON domains from Ecoregions and apply to FHM
colnames(NEON_FHM)[c(8, 5, 6,18,13,12, 11, 16,21:23)] = colnames(NEON)[c(1,5:10, 16:19)]
domains = fread("./indir/ecoregions_domain.csv") %>% select(domainID, US_L4CODE) %>% unique
NEON_FHM = inner_join(NEON_FHM, domains)
# clean features that are too noisy
NEON_FHM$stemDiameter = sqrt(NEON_FHM$stemDiameter)
NEON_FHM$DBI = log(NEON_FHM$DBI)
NEON_FHM$DENSITY = log(NEON_FHM$DENSITY)
NEON_FHM$maxCrownDiameter = log(NEON_FHM$maxCrownDiameter)
NEON_FHM$ninetyCrownDiameter = log(NEON_FHM$ninetyCrownDiameter)


#scale all numeric variables
vars_to_scale = c("stemDiameter","height", "DBI", "DENSITY","BA",
                  "maxCrownDiameter")

scaleListFHM = scale(NEON_FHM[vars_to_scale])
scaleListFHM <- list(scale = attr(scaleListFHM, "scaled:scale"),
                     center = attr(scaleListFHM, "scaled:center"))


NEON_FHM$STATURE = NEON_FHM$height
NEON_FHM$H_VAR = mean(NEON$H_VAR)
NEON_FHM$canopyPosition = "Partially shaded"
test_fhm = NEON_FHM %>% select("stemDiameter","height", "DBI", "DENSITY","BA",
                                 "maxCrownDiameter", "taxonID", "domainID")
NEON_FHM = NEON_FHM[complete.cases(NEON_FHM),]
bayes_R2(crown_size_from_structure, newdata = data.frame(NEON_FHM), allow_new_levels = T)
