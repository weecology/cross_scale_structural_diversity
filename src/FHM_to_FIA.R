# create dbh and crown allometries from FHM to be plugged in FIA data
library(data.table)
library(tidyverse)
FHM = fread("./indir/FHM_STRUCTURE_PLOT.csv")
summary(FHM)
PHENOLOGY = fread("./indir/FHM_NEON_phylogeny.csv")
