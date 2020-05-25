combinedGSEA <- function(v, idx, design, contrasts){
  # Performs the following gene set enrichment testing methods: fry, fgsea, camera,
  # and combines the raw p-values obtained with each method to obtain an overall
  # "ensemble" p-value.
  #
  # Args:
  #   v: voom object obtained from using limma::voom() or limma::voomWithQualityWeights().
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


  # Unlike fry and camera, fgsea needs a ranked list of genes, which we will
  # prepare from t-statistics obtained from fitting a linear model.
  fit <- lmFit(v, design) %>%
    contrasts.fit(contrasts) %>%
    eBayes(robust = TRUE)

  # fgsea also requires the index to be gene names rather than indices like
  # fry and camera.
  fgsea_idx <- lapply(idx, function(x){
    v %>% magrittr::extract(x,) %>% rownames()
  })


  # Run the gene set enrichment testing methods for each contrast.
  res <- colnames(contrasts) %>%
    lapply(function(x){

      # For each contrast, we need to extract out a ranked list of genes for
      # fgsea, ranked by t-statistic for DE in that particular comparison.
      fgsea_ranks <- topTable(fit, coef = x, number = Inf) %>%
        #rownames_to_column("ensembl_gene_id") %>%
        dplyr::arrange(t)
      fgsea_ranks <- fgsea_ranks$t %>% set_names(fgsea_ranks$ensembl_gene_id)

      # The methods used are limma::fry, limma::camera, and fgsea.
      res2 <- list(
        fry = fry(
          v,
          idx,
          design = design,
          contrast = contrasts[, x]
          ) %>%
          rownames_to_column("Geneset") %>%
          dplyr::mutate(pval = PValue.Mixed),
        mroast = mroast(
          v,
          idx,
          design = design,
          contrast = contrasts[, x]
        ) %>%
          rownames_to_column("Geneset") %>%
          dplyr::mutate(pval = PValue.Mixed),
        camera = camera(
          v,
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
    return(x)
  }) %>%
    magrittr::extract(. != "mroast") %>%
    lapply(function(x){
      x %>% do.call("rbind", .) %>%
        dplyr::group_by(Geneset) %>%
        dplyr::summarise(wilkinsonp = metap::wilkinsonp(pval, r = 1)$p) %>%
        dplyr::mutate(fdr = p.adjust(wilkinsonp, method = "fdr"),
                      bonferroni = p.adjust(wilkinsonp, method = "bonferroni"))%>%
        dplyr::arrange(wilkinsonp)
    })

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
