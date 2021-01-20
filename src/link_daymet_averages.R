#load FIA data 
aa = data.table::fread("../../FIA_climate/DAYMET/climate_features_county.csv")
coords = data.table::fread("../../FIA_climate/my_sites.csv")
colnames(aa)[9] = "ID"
dd = aa$ID %>% table %>% data.frame
dd = dd %>% filter(Freq == 12)
coor_ = coords %>% filter(ID %in% unique(dd[,1]))
dd = inner_join(coords, aa)
write_csv(dd, "../data/climate_monthly_averages_2015.csv")
sum(is.na(dd$LAT))

ee = dd$ID %>% table %>% data.frame
ee %>% filter(Freq >12)
