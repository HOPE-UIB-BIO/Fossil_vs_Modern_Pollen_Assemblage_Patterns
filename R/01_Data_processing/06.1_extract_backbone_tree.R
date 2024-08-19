
#----------------------------------------------------------#
# Fossil pollen data can reconstruct robust spatial patterns of biodiversity 
#                        in the past
#
#                         K. Bhatta 
#
#                           2024
#----------------------------------------------------------#

#----------------------------------------------------------#
#
# Extract time calibrated backbone tree for estimating phylogenetic relatedness ----
#
#----------------------------------------------------------#
#-----------------------------------------------#
# Load configuration ----
#-----------------------------------------------#
source("R/00_Config_file.R")


# A. Surface pollen data ----

#-----------------------------------------------#
# Load harmonised filtered data of surface pollen samples ----
#-----------------------------------------------#
surface_pollen_filtered_clim_zone <- 
  readr::read_rds("Inputs/Data/surface_pollen_filtered_030124.rds")

#-----------------------------------------------#
# Load full backbone tree ----
#-----------------------------------------------#
# NOTE: run this only for the first time!
# After we create a regional tree, no need to re-run this code
# We used maximum likelihood angiosperm phylogeny (NEWICK format) estimated with 
#  RAxML: Data1_RAxMLTree.tre by Ramírez-Barahona et al. 2020 The delayed and 
#  geographically heterogeneous diversification of flowering plant families. 
#  Nature Ecology and Evolution. 4: 1232-1238"

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
  colnames(surface_pollen_filtered_clim_zone[,-(1:7)]) %>%
        enframe(name = NULL,
                value = "taxa") %>%
  distinct() %>% 
  arrange(taxa) # 105 taxa in the datasets

req_taxa %>% 
  dplyr::filter(!taxa %in% mytree$tip.label )
# Viburnaceae are not in the backbone tree
# Viburnaceae changed to Adoxaceae (synonym)
# It was corrected in the harmonisation table

#-----------------------------------------------#
# Prune the full backbone tree ----
#-----------------------------------------------#
# Taxa in the backbone tree that are absent in the datasets  
drop_list <- 
  mytree$tip.label[!mytree$tip.label %in% req_taxa$taxa] # 338

# Remove all the families that do not occur in our datasets using drop.tip() 
final_tree <- ape::drop.tip(mytree, drop_list) # 105 families     

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


# B. Fossil pollen data ----

#-----------------------------------------------#
# Load harmonised filtered data of fossil pollen samples ----
#-----------------------------------------------#
top_fossil_pollen_500yr_filtered_clim_zone <-
  readr::read_rds("Inputs/Data/top_fossil_pollen_500yr_filtered_080124.rds") 

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
req_taxa_500yr <-
  colnames(top_fossil_pollen_500yr_filtered_clim_zone[,-(1:8)]) %>%
  enframe(name = NULL,
          value = "taxa") %>%
  distinct() %>% 
  arrange(taxa) # 79 taxa in the datasets

req_taxa_500yr %>% 
  dplyr::filter(!taxa %in% mytree$tip.label )
# Viburnaceae are not in the backbone tree
# Viburnaceae changed to Adoxaceae (synonym)
# Corrected in the harmonisation table

#-----------------------------------------------#
# Prune the full backbone tree ----
#-----------------------------------------------#
# Taxa in the backbone tree that are absent in the datasets  
drop_list_500yr <- 
  mytree$tip.label[!mytree$tip.label %in% req_taxa_500yr$taxa] # 364

# Remove all the families that do not occur in our datasets using drop.tip() 
final_tree_top_fossil_samples_500yr <- ape::drop.tip(mytree, drop_list_500yr) # 81 families     

#-----------------------------------------------#
# Save the regional backbone tree ----
#-----------------------------------------------#
ape::write.tree(
  final_tree_top_fossil_samples_500yr, 
  paste(
    "Inputs/Data/Ramirez_Barahona_etal_2020_phylogeny/",
    "pruned_Ramirez_Barahona_etal_2020_raxml_top_fossil_samples.tre",
    sep = ""
  )
)

