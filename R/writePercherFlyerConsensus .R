
#raw_pf_data_location="D:/GitProjects/odomatic-wings-analysis/data/percher_flyer.csv"
#save_location="D:/GitProjects/dragonfly-wings-bodies-na/data/flight_styles.csv"
writePercherFlyerConsensus <- function(save_location, raw_pf_data_location="D:/GitProjects/odomatic-wings-analysis/data/percher_flyer.csv"){
  library(dplyr)
  # add percher flyer
  pf <- read.csv(raw_pf_data_location)
  colnames(pf) <- c("Taxon","John_Percher","John_Flyer","John_Intermediate","John_Reference",
                    "Jess_Percher","Jess_Flyer","Jess_Intermediate","Jess_Reference")
  pf <- pf[-1,]
  consensus <- c()
  for(i in 1:nrow(pf)){
    row <- pf[i,]
    if(row$John_Percher == "x" & row$Jess_Percher == "X"){
      consensus <- c(consensus,"percher")
    }
    else if(row$John_Flyer == "x" & row$Jess_Flyer == "X"){
      consensus <- c(consensus,"flyer")
    }
    else if(row$John_Intermediate == "x" & row$Jess_Intermediate == "X"){
      consensus <- c(consensus,"intermediate")
    }
    else{
      consensus <- c(consensus,NA)
    }
  }
  library(stringr)
  pf$species <- paste(str_split_fixed(pf$Taxon," ", 3)[,1],str_split_fixed(pf$Taxon," ", 3)[,2])
  pf$flight_type <- consensus  
  pf <- pf %>% select(c(species,flight_type))
  
  write.csv(pf, save_location,row.names = FALSE)
} 