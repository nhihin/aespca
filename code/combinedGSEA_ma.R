combinedGSEA_ma <- function(exprs, fit, idx, design, contrasts){
  # MODIFIED FOR MICROARRAY ANALYSIS WHICH DOESN'T USE VOOM
  # Performs the following gene set enrichment testing methods: fry, fgsea, camera,
  # and combines the raw p-values obtained with each method to obtain an overall
  # "ensemble" p-value. 
  #
  # Args:
  #   exprs: expression matrix
  #   fit: fit object
  #   idx: Named list of gene sets, with numeric values corresponding to rows in v. 
  #        Can be prepared using limma::ids2indices(). 
  #   design: design matrix prepared during limma workflow. 
  #   contrasts: contrasts matrix prepared during limma workflow. 
  #
  # Returns:
  #   A nested list of the results for each contrast. The indivTest
  #   slot gives the results returned by each gene set enrichment testing method. 
  #   The combTest slot gives the results from combining p-values from the individual
  #   methods, along with FDR and Bonferroni-adjusted p-values. 
  
  
  # fgsea  requires the index to be gene names rather than indices like 
  # fry and camera. 
  fgsea_idx <- lapply(idx, function(x){
    exprs %>% magrittr::extract(x,) %>% rownames() 
  })
  
  
  # Run the gene set enrichment testing methods for each contrast.
  res <- colnames(contrasts) %>%
    lapply(function(x){
      
      # For each contrast, we need to extract out a ranked list of genes for
      # fgsea, ranked by t-statistic for DE in that particular comparison.
      fgsea_ranks <- topTable(fit, coef = x, number = Inf) %>% 
        dplyr::arrange(t) %>%
        set_rownames(rownames(exprs))
      fgsea_ranks <- fgsea_ranks$t %>% set_names(fgsea_ranks%>%rownames)
      
      # The methods used are limma::fry, limma::camera, and fgsea.
      res2 <- list(
        fry = fry(
          exprs, 
          idx, 
          design = design, 
          contrast = contrasts[, x]
        ) %>% 
          rownames_to_column("Geneset") %>%
          dplyr::mutate(pval = PValue.Mixed),
        
        mroast = mroast(
          exprs, 
          idx, 
          design = design, 
          contrast = contrasts[, x]
        ) %>% 
          rownames_to_column("Geneset") %>%
          dplyr::mutate(pval = PValue.Mixed),
        
        camera = camera(
          exprs, 
          idx,
          design = design,  
          contrast = contrasts[, x],
          allow.neg.cor = TRUE,
          inter.gene.cor = NA
        ) %>%
          rownames_to_column("Geneset") %>%
          dplyr::mutate(pval = PValue),
        
        fgsea = fgsea(
          fgsea_idx, 
          fgsea_ranks, 
          nperm = 1e5, 
          nproc = 6
        ) %>%
          dplyr::rename(Geneset = pathway)
      )
      
    }) %>% set_names(colnames(contrasts))
  
  # Extract out only the relevant Geneset and p-value columns so that 
  # the rows of each data.frame can be bound together. 
  res3 <- res %>% lapply(function(x){
    x %<>% lapply(function(y){
      y %<>% dplyr::select(Geneset, pval)
      return(y)
    })
    return(x) %>% .[c("fry","fgsea","camera")] 
  }) %>%
    lapply(function(x){
      x %>% do.call("rbind", .) %>%
        as.data.frame %>%
        dplyr::group_by(Geneset) }) %>%
    lapply(function(x){
      x %>% dplyr::summarise(wilkinsonp = metap::wilkinsonp(pval, r = 1)$p) %>%
        dplyr::mutate(fdr = p.adjust(wilkinsonp, method = "fdr"),
                      bonferroni = p.adjust(wilkinsonp, method = "bonferroni"))%>%
        dplyr::arrange(wilkinsonp)
    })
    #   %>%
    #     dplyr::summarise(wilkinsonp = metap::wilkinsonp(pval, r = 1)$p) %>%
        # dplyr::mutate(fdr = p.adjust(wilkinsonp, method = "fdr"),
        #               bonferroni = p.adjust(wilkinsonp, method = "bonferroni"))%>%
        # dplyr::arrange(wilkinsonp)
    # })
  
  # Perform adjustment for multiple testing
  res4 <- res3 %>% bind_rows(.id = "id") %>%
    mutate(fdr = p.adjust(wilkinsonp, "BH"),
           bonferroni = p.adjust(wilkinsonp, "holm")) %>%
    split(f=.$id)
  
  # Return relevant results
  results <- list(
    indivTest = res,
    combTest = res4
  )
  
  return(results)
}