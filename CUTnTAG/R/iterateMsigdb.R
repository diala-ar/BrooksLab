#' iterateMsigdb
#' @description Iterates through all levels of msigdb that are relevant
#' for geneset analysis. With the levels loaded, the user can specify
#' the fun() param to pass in their own analysis for these levels.
#' 
#' @param species Either 'Homo sapiens' or 'Mus musculus'
#' @param mlvl a list of all msigdb levels to analyze. For example:
#'    list('H'=list(NULL),
#'         'C2'=list('CP:REACTOME'))
#' @param fun A function passed with a defined msig_ds param passed in
#' from the msigdb levels (msig_ds=msig_ds), meant to analyze the data using
#' these values; for example:
#'    gseaFun <- function(msig_db, lfc_v){
#'      GSEA(sort(na.omit(lfc_v), decreasing = T), TERM2GENE = msig_ds, pvalueCutoff = 1)
#'    }
#' @param ... pass-through params to your function
#'
#' @examples
#' gseaFun <- function(msig_db, lfc_v){
#'   GSEA(sort(na.omit(lfc_v), decreasing = T), TERM2GENE = msig_ds, pvalueCutoff = 1)
#' }
#' iterateMsigdb(species='Homo sapiens', fun=gseaFun, lfc_v = lfc_v)
iterateMsigdb <- function(species, msig_lvls=NULL, fun, ...){
  ## Select the significant genes and filters
  if(is.null(msig_lvls)) {
    warning("Running msigdb levels: H, C2-CP:REACTOM, C5-GO:BP,GO:CC,GO:MF, C8")
    msig_lvls <- list('H'=list(NULL),                       # hallmark gene sets
                      'C2'=list('CP:REACTOME'),             # curated gene sets
                      'C5'=list('GO:BP', 'GO:CC', 'GO:MF'), # ontology gene sets
                      'C8'=list(NULL))                      # cell type signature gene sets
  }
  
  ##  Load MSigDb database
  accepted_species <- c('Mus musculus', 'Homo sapiens')
  stopifnot("species not Homo sapiens or Mus musculus"=species %in% accepted_species)
  if(species == 'Homo sapiens'){
    library(org.Hs.eg.db)
    annotation <- 'org.Hs.eg.db'
    org_code <- 'hsa'
    genome <- org.Hs.eg.db
  } else if(species == 'Mus musculus'){
    library(org.Mm.eg.db)
    annotation <- 'org.Mm.eg.db'
    org_code <- 'mma'
    genome <- org.Mm.eg.db
  }
  
  
  msig_obj <- lapply(names(msig_lvls), function(mlvl){
    if(mlvl == 'custom'){
      sub_obj <- lapply(names(msig_lvls[[mlvl]]), function(sublvl){
        msig_ds <- data.frame("gs_name"=sublvl,
                              "entrez_gene"=msig_lvls[[mlvl]][[sublvl]])
        
        obj <- fun(msig_ds=msig_ds, ...)
        return(obj)
      })
      names(sub_obj) <- names(msig_lvls[[mlvl]])
    } else {
      sub_obj <- lapply(msig_lvls[[mlvl]], function(sublvl){
        print(paste0(">", mlvl, ":", sublvl, "..."))
        msig_ds <- msigdbr(species = species, category = mlvl, subcategory = sublvl) %>%
          dplyr::select(gs_name, entrez_gene) %>%
          as.data.frame()
        
        # overrepresentation analysis
        obj <- fun(msig_ds=msig_ds, ...)
        return(obj)
      })
      ids <- unlist(msig_lvls[[mlvl]])
      ids <- if(is.null(ids)) 'base' else ids
      names(sub_obj) <- ids
    }
    
    return(sub_obj)
  })
  names(msig_obj) <- names(msig_lvls)
  
  return(msig_obj)
}

#' GSEA Wrapper function
#' @description  GSEA wrapper function to insert as the function for 
#' iterateMsigdb()
#' 
#' @param msig_ds parameter passed in via iterateMsigdb()
#' @param lfc_v vector of log-fold change values 
#'
#' @examples
#' lfc_v <- setNames(reslfc$log2FoldChange, reslfc$entrez)
#' gseas <- iterateMsigdb(species=species, fun=gseaFun, 
#'                        lfc_v=lfc_v)
gseaFun <- function(msig_ds, lfc_v){
  gsea <- tryCatch({
    GSEA(sort(na.omit(lfc_v), decreasing = T), 
         TERM2GENE = msig_ds, pvalueCutoff = 1)
  }, error=function(e){NULL})
  return(gsea)
}

#' Gene-set Enrichment Wrapper function
#' @description  GO wrapper function to insert as the function for 
#' iterateMsigdb()
#'
#' @param msig_ds parameter passed in via iterateMsigdb()
#' @param entrez_genes vector of entrez-id genes
#'
#' @examples
#' oras <- iterateMsigdb(species='Mus musculus', fun=oraFun, 
#' entrez_genes=module_genes$entrez)
#' oras <- lapply(oras, summarizeOra, keep_genes=TRUE)
oraFun <- function(msig_ds, entrez_genes){
  require(clusterProfiler)
  # overrepresentation analysis
  sig_ora <- tryCatch({
    enricher(gene = na.omit(entrez_genes), TERM2GENE = msig_ds)@result
  }, error=function(e){NULL})
  return(sig_ora)
}

#' ssGSEA Wrapper function
#' @description  ssGSEA wrapper function to insert as the function for 
#' iterateMsigdb()
#' 
#' @param msig_ds parameter passed in via iterateMsigdb()
#' @param lfc_v matrix of expression with genes in entrez-id form
#' in the rows and samples in the column
#'
#' @examples
#' tcnts <- as.data.frame(assay(counts))
#' rownames(tcnts) <- ens2entrez_ids[rownames(tcnts)]
#' ssgseas <- iterateMsigdb(species=species, fun=ssGseaFun, 
#'                          lfc_v=tcnts,  msig_lvls = list('H'=list(NULL)))
ssGseaFun <- function(msig_ds, lfc_v, ss_method='ssgsea'){
  require(GSVA)
  ssgsea <- tryCatch({
    sig_ens_gs <- split(setNames(msig_ds$entrez_gene, msig_ds$entrez_gene), 
                        f=msig_ds$gs_name)
    gsva(lfc_v, sig_ens_gs, verbose=FALSE, method=ss_method)
  }, error=function(e){NULL})
  return(ssgsea)
}


#' summarizeOra
#' @description Takes a list of ORA analysis from enricher
#' using the msigdb database, adn reduces it to a cleaned up
#' data structure without all the gene names or extra data
#' 
#' @param ora A list of datadrames from enricher ORA analysis
#' @param qcutoff qcutoff to show only genesets below this signficance
summarizeOra <- function(ora, qcutoff=0.05, keep_genes=FALSE){
  keep_idx <- c(1,3,4,5,6,7,9)
  if(keep_genes) keep_idx <- sort(c(keep_idx, 8))
  ora_df <- as.data.frame(do.call(rbind, ora)[,keep_idx])
  splitFrac <- function(i){
    sapply(strsplit(i, split="/"), function(j){
      round(as.integer(j[1]) / as.integer(j[2]), 3)
    })
  }
  #https://en.wikipedia.org/wiki/Cohen%27s_h
  gener <- 2*asin(sqrt(splitFrac(ora_df$GeneRatio)))
  bgr <- 2*asin(sqrt(splitFrac(ora_df$BgRatio)))
  ora_df$cohens_h <- abs(gener - bgr)
  
  rownames(ora_df) <- NULL
  ora_df <- ora_df[ora_df$`p.adjust`<qcutoff, ]
  return(ora_df)
}