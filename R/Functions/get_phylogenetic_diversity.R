#-----------------------------------------------#
# Function for estimating phylogenetic diversity ----
#-----------------------------------------------#
get_phylogenetic_diversity <- 
  function(counts,
           backbone_tree,
           type = "mpd",
           null.model = "phylogeny.pool", 
           abundance.weighted = TRUE,
           runs = 9999,
           ...) {
    
  dat <- 
    counts %>%
    as.data.frame() %>% 
    tibble::column_to_rownames("sample_id")
  
  final_tree <- backbone_tree
  # Make a vector of families from the tree that are missing from the counts data
  drop_list <- 
    final_tree$tip.label[!final_tree$tip.label %in% colnames(dat)]
  
  # Remove all the families that do not occur in our sample using drop.tip() 
  #  function in the 'ape' package.
  pruned_tree <- ape::drop.tip(final_tree, drop_list) 
  
  # Arrange taxa in the same order as the pruned tree
  data_ordered <- dat[,c(pruned_tree$tip.label)]
  
  # Calculate MPD and MNTD
  # Create cophenetic distance from the pruned tree.
  phy_dist <- cophenetic(pruned_tree) 
  
  # MPD using ses.mpd() function in the 'picante' package
  
  if (type == "mpd") {
    set.seed(1234)
    mpd_phylogeny <- 
      picante::ses.mpd(
        data_ordered,
        phy_dist,
        null.model = null.model, 
        abundance.weighted = abundance.weighted,
        runs
        ) %>% 
     tibble:: rownames_to_column("sample_id") %>% 
      tibble::as_tibble() %>% 
      dplyr::rename(ses_mpd = mpd.obs.z)
    
    return(mpd_phylogeny) 
    
  } else if (type == "mntd") {
    
    set.seed(1234)
    mntd_phylogeny <- 
      picante::ses.mntd(
        data_ordered,
        phy_dist,
        null.model = null.model, 
        abundance.weighted = abundance.weighted,
        runs = runs
        ) %>% 
      tibble::rownames_to_column("sample_id") %>% 
      tibble::as_tibble() %>% 
      dplyr::rename(ses_mntd = mntd.obs.z)
    
    return(mntd_phylogeny) 
  } else {
    stop("type must be either 'mpd' or 'mntd'")
  }
}
