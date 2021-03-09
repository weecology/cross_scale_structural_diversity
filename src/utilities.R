clean_typos_taxonID <- function(species_classification_data){
  #make it a function but for testing ok replicating upper loop
  species_classification_data$taxonID[species_classification_data$taxonID == "ABLAL"] = "ABAM"
  species_classification_data$taxonID[species_classification_data$taxonID == "ACFA"] = "ACSA2"
  species_classification_data$taxonID[species_classification_data$taxonID == "ACNEN"] = "ACNE2"
  species_classification_data$taxonID[species_classification_data$taxonID == "ACSAS"] = "ACSA3"
  species_classification_data$taxonID[species_classification_data$taxonID == "ACRUR"] = "ACRU"
  species_classification_data$taxonID[species_classification_data$taxonID == "ARVIM"] = "OTHER"
  species_classification_data$taxonID[species_classification_data$taxonID == "BETUL"] = "BEPA"
  species_classification_data$taxonID[species_classification_data$taxonID == "BEPAP"] = "BEPA"
  species_classification_data$taxonID[species_classification_data$taxonID == "BOURR"] = "OTHER"
  species_classification_data$taxonID[species_classification_data$taxonID == "BOURRSPP"] = "OTHER"
  species_classification_data$taxonID[species_classification_data$taxonID == "BUBU"] = "OTHER"
  species_classification_data$taxonID[species_classification_data$taxonID == "CAOV3"] = "CAOV2"
  species_classification_data$taxonID[species_classification_data$taxonID == "CELE2"] = "OTHER"
  species_classification_data$taxonID[species_classification_data$taxonID == "DIOSP"] = "OTHER"
  species_classification_data$taxonID[species_classification_data$taxonID == "DIVI5"] = "OTHER"
  species_classification_data$taxonID[species_classification_data$taxonID == "FRAXI"] = "FRAPE"
  species_classification_data$taxonID[species_classification_data$taxonID == "HALES"] = "OTHER"
  species_classification_data$taxonID[species_classification_data$taxonID == "ILAN"] = "OTHER"
  species_classification_data$taxonID[species_classification_data$taxonID == "JUGLA"] = "JUNI"
  species_classification_data$taxonID[species_classification_data$taxonID == "JUOS"] = "JUVI"
  species_classification_data$taxonID[species_classification_data$taxonID == "LOMA6"] = "OTHER"
  species_classification_data$taxonID[species_classification_data$taxonID == "AGNO"] = "OTHER"
  species_classification_data$taxonID[species_classification_data$taxonID == "MEPO5"] = "OTHER"
  species_classification_data$taxonID[species_classification_data$taxonID == "NYBI"] = "OTHER"
  species_classification_data$taxonID[species_classification_data$taxonID == "NYSY"] = "OTHER"
  species_classification_data$taxonID[species_classification_data$taxonID == "OXYDE"] = "OTHER"
  species_classification_data$taxonID[species_classification_data$taxonID == "PARK12"] = "OTHER"
  species_classification_data$taxonID[species_classification_data$taxonID == "PATO2"] = "OTHER"
  species_classification_data$taxonID[species_classification_data$taxonID == "PICEA"] = "PIAB"
  species_classification_data$taxonID[species_classification_data$taxonID == "PINUS"] = "PIPA2"
  species_classification_data$taxonID[species_classification_data$taxonID == "PIPOS"] = "PIPA2"
  species_classification_data$taxonID[species_classification_data$taxonID == "POPUL"] = "POTR5"
  species_classification_data$taxonID[species_classification_data$taxonID == "PRAM"] = "PRAV"
  species_classification_data$taxonID[species_classification_data$taxonID == "PRSES"] = "PRAV"
  species_classification_data$taxonID[species_classification_data$taxonID == "PRVE"] = "PRAV"
  species_classification_data$taxonID[species_classification_data$taxonID == "QUERC"] = "QUAL"
  species_classification_data$taxonID[species_classification_data$taxonID == "QUHE"] = "QUAL"
  species_classification_data$taxonID[species_classification_data$taxonID == "QUHE2"] = "QUAL"
  species_classification_data$taxonID[species_classification_data$taxonID == "QUMA13"] = "QUAL"
  species_classification_data$taxonID[species_classification_data$taxonID == "SASSA"] = "OTHER"
  species_classification_data$taxonID[species_classification_data$taxonID == "SWMA2"] = "OTHER"
  species_classification_data$taxonID[species_classification_data$taxonID == "TABR2"] = "OTHER"
  species_classification_data$taxonID[species_classification_data$taxonID == "TEDA"] = "OTHER"
  species_classification_data$taxonID[species_classification_data$taxonID == "ULMUS"] = "OTHER"
  species_classification_data$taxonID[species_classification_data$taxonID == "ULPU"] = "OTHER"
  species_classification_data$taxonID[species_classification_data$taxonID == "VIBUR"] = "OTHER"
  species_classification_data$taxonID[species_classification_data$taxonID == "PICOL"] = "PICO"
  species_classification_data$taxonID[species_classification_data$taxonID == "PSMEM"] = "PSME"

  
  return(species_classification_data)
}


plot_errors_relations_to_features <- function(offset){
  errors_relations = melt(offset, id.vars="delta")
  ggplot(errors_relations, aes(x=value, y = delta))+
    geom_point()+facet_wrap(variable~., scales = "free")
}