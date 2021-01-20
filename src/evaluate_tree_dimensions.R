library(bien)
library(data.table)
library(tidyverse)

get_diameter = function(dt, max_min){
  side1 = abs(dt[["top"]] - dt[["bottom"]])
  side2 = abs(dt[["left"]] - dt[["right"]])
  if(max_min == "max"){
    diam = pmax(side1, side2)
  }else if(max_min == "min"){
    diam = pmin(side1, side2)
  }else{
    warning("no function used to define with diameter to be retrieved")
  }
}

DFor = fread("./indir/deepForest_full_uniques.csv",
             select = c("individualID", "left","bottom","right","top",
                        "score","label","height","area"))
colnames(DFor)[8] = "df_height"
vst = fread("./indir/vst_dec2020_3m.csv")

DFor[["df_max_diameter"]] = get_diameter(DFor, "max")
DFor[["df_min_diameter"]] = get_diameter(DFor, "min")


jStructure = inner_join(vst, DFor) %>%
  select(individualID, plotID, subplotID, siteID,
         taxonID, nlcdClass, stemDiameter, height, maxCrownDiameter,
         ninetyCrownDiameter,df_max_diameter,df_min_diameter,
         score, df_height, area)

#jStructure = jStructure %>% filter(score > 0.4)
ggplot(data = jStructure, aes(x = maxCrownDiameter, y = df_max_diameter)) +
  geom_point()+geom_abline(slope = 1) + theme_bw()

ggplot(data = jStructure, aes(x = log(maxCrownDiameter), y = log(df_max_diameter))) +
  geom_point()+geom_abline(slope = 1) + theme_bw()

ggplot(data = jStructure, aes(x = ninetyCrownDiameter, y = df_min_diameter)) +
  geom_point()+geom_abline(slope = 1) + theme_bw()

ggplot(data = jStructure, aes(x = log(ninetyCrownDiameter), y = log(df_min_diameter))) +
  geom_point()+geom_abline(slope = 1) + theme_bw()


# Subplot average
jStructure = inner_join(vst, DFor) %>%
  select(individualID, plotID, subplotID, siteID,
         taxonID, nlcdClass, stemDiameter, height, maxCrownDiameter,
         ninetyCrownDiameter,df_max_diameter,df_min_diameter,
         score, df_height, area)

jStructure = jStructure %>% filter(abs(df_height - height) < 4)

jAllo = jStructure %>% select(stemDiameter, maxCrownDiameter,
          ninetyCrownDiameter,df_max_diameter,df_min_diameter)
jAllo = melt(jAllo, id.vars = "stemDiameter")
# train a dbh - cr allometry from field and df
ggplot(data = jAllo, aes(x = sqrt(value), y = log(stemDiameter), color = variable)) +
  geom_point(alpha = 0.3, size = 0.3)+ theme_bw() + geom_smooth(method = "lm") + facet_wrap(.~variable)+
  theme(legend.position = "bottom")

jave_Structure = jStructure %>%
  group_by(plotID, subplotID) %>%
  summarize_if(is.numeric, mean)

jave_Allo = jave_Structure %>% ungroup %>%
  select(stemDiameter, maxCrownDiameter,
                              ninetyCrownDiameter,df_max_diameter,df_min_diameter)
jave_Allo = melt(jave_Allo, id.vars = "stemDiameter")

# train a dbh - cr allometry from field and df
ggplot(data = jave_Allo, aes(x = sqrt(value), y = log(stemDiameter), color = variable)) +
  geom_point(alpha = 0.3, size = 0.3)+ theme_bw() + geom_smooth(method = "lm") + facet_wrap(.~variable)+
  theme(legend.position = "bottom")



