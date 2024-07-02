#----------------------------------------------------------#
#
#       Latitudinal analysis of phylogenetic dispersion
#
#               Surface sample pollen data 
#          
#        Estimate phylogenetic relatedness (MPD, MNTD) ----
#----------------------------------------------------------#
#-----------------------------------------------#
# Load configuration ----
#-----------------------------------------------#
source("R/00_Config_file.R")

#-----------------------------------------------#
# Load harmonised filtered data of surface pollen samples ----
#-----------------------------------------------#
surface_pollen_filtered_clim_zone <- 
  read_rds("Inputs/Data/surface_pollen_filtered_030124.rds")

#-----------------------------------------------#
# Load full backbone tree ----
#-----------------------------------------------#
# NOTE: run this only for the first time!
# After we create a regional tree, no need to re-run this code

mytree <-
  ape::read.tree(
    paste(
      "Inputs/Data/Ramirez_Barahona_etal_2020_phylogeny/",
      "Ramirez_Barahona_etal_2020_raxml_processed.tre",
      sep = ""
    ) 
  )

# To prune the backbone tree to keep only keep the taxa from Asian data, we need 
# the list of taxa present in the data.
req_taxa <-
  colnames(data_for_phylodiv_estimation[,-(1:7)]) %>%
        enframe(name = NULL,
                value = "taxa") %>%
  distinct() %>% 
  arrange(taxa) # 105 taxa in the datasets

req_taxa %>% 
  dplyr::filter(!taxa %in% mytree$tip.label )
# Viburnaceae are not in the backbone tree
# Viburnaceae changed to Adoxaceae (synonym)

#-----------------------------------------------#
# Prune the full backbone tree ----
#-----------------------------------------------#
# Taxa in the backbone tree that are absent in the datasets  
drop_list <- 
  mytree$tip.label[!mytree$tip.label %in% req_taxa$taxa] # 338

# Remove all the families that do not occur in our datasets using drop.tip() 
final_tree <- ape::drop.tip(mytree, drop_list) # 148 families     

#-----------------------------------------------#
# Save the regional backbone tree ----
#-----------------------------------------------#
ape::write.tree(
  final_tree, 
  paste(
    "Inputs/Data/Ramirez_Barahona_etal_2020_phylogeny/",
    "pruned_Ramirez_Barahona_etal_2020_raxml_surface_samples.tre",
    sep = ""
  )
)

#-----------------------------------------------#
# Load harmonised filtered data of top 1000 yr fossil pollen samples ----
#-----------------------------------------------#
top_fossil_pollen_1000yr_filtered <- 
  read_rds("Inputs/Data/top_fossil_pollen_1000yr_filtered_080124.rds")

#-----------------------------------------------#
# Load full backbone tree ----
#-----------------------------------------------#
# NOTE: run this only for the first time!
# After we create a regional tree, no need to re-run this code

mytree <-
  ape::read.tree(
    paste(
      "Inputs/Data/Ramirez_Barahona_etal_2020_phylogeny/",
      "Ramirez_Barahona_etal_2020_raxml_processed.tre",
      sep = ""
    ) 
  )

# To prune the backbone tree to keep only keep the taxa from Asian data, we need 
# the list of taxa present in the data.
req_taxa_1000yr <-
  colnames(top_fossil_pollen_1000yr_filtered[,-(1:8)]) %>%
  enframe(name = NULL,
          value = "taxa") %>%
  distinct() %>% 
  arrange(taxa) # 81 taxa in the datasets

req_taxa_1000yr %>% 
  dplyr::filter(!taxa %in% mytree$tip.label )
# Viburnaceae are not in the backbone tree
# Viburnaceae changed to Adoxaceae (synonym)

#-----------------------------------------------#
# Prune the full backbone tree ----
#-----------------------------------------------#
# Taxa in the backbone tree that are absent in the datasets  
drop_list_1000yr <- 
  mytree$tip.label[!mytree$tip.label %in% req_taxa_1000yr$taxa] # 362

# Remove all the families that do not occur in our datasets using drop.tip() 
final_tree_top_fossil_samples_1000yr <- ape::drop.tip(mytree, drop_list_1000yr) # 81 families     
#-----------------------------------------------#
# Save the regional backbone tree ----
#-----------------------------------------------#
ape::write.tree(
  final_tree_top_fossil_samples_1000yr, 
  paste(
    "Inputs/Data/Ramirez_Barahona_etal_2020_phylogeny/",
    "pruned_Ramirez_Barahona_etal_2020_raxml_top_fossil_samples.tre",
    sep = ""
  )
)

#-----------------------------------------------#
# Load harmonised filtered data of modern vegetation occurrences (GBIF) ----
#-----------------------------------------------#
taxon_data_transform <- 
  read_rds("Inputs/Data/data_gbif_for_phylodiversity_160124.rds")

#-----------------------------------------------#
# Load full backbone tree ----
#-----------------------------------------------#
# NOTE: run this only for the first time!
# After we create a regional tree, no need to re-run this code

mytree <-
  ape::read.tree(
    paste(
      "Inputs/Data/Ramirez_Barahona_etal_2020_phylogeny/",
      "Ramirez_Barahona_etal_2020_raxml_processed.tre",
      sep = ""
    ) 
  )

# To prune the backbone tree to keep only keep the taxa from Asian data, we need 
# the list of taxa present in the data.
req_taxa <-
  colnames(taxon_data_transform) %>%
  enframe(name = NULL,
          value = "taxa") %>%
  distinct() %>% 
  arrange(taxa) # 274 taxa in the datasets

req_taxa %>% 
  dplyr::filter(!taxa %in% mytree$tip.label )
# All taxa are in the backbone tree

#-----------------------------------------------#
# Prune the full backbone tree ----
#-----------------------------------------------#
# Taxa in the backbone tree that are absent in the datasets  
drop_list <- 
  mytree$tip.label[!mytree$tip.label %in% req_taxa$taxa] # 169

# Remove all the families that do not occur in our datasets using drop.tip() 
final_tree <- ape::drop.tip(mytree, drop_list) # 274 families     

#-----------------------------------------------#
# Save the regional backbone tree ----
#-----------------------------------------------#
ape::write.tree(
  final_tree, 
  paste(
    "Inputs/Data/Ramirez_Barahona_etal_2020_phylogeny/",
    "pruned_Ramirez_Barahona_etal_2020_raxml_modern_vegetation.tre",
    sep = ""
  )
)

