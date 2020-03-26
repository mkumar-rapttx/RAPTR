
#' Perform gene set enrichment analysis
#'
#' @param geneSets a list of gene name vectors
#' @param expressionVals a vector of gene scores such as differential expression values
#' @param genes a vector of gene names in the same order as expressionVals
#' @param minSetSize minimum number of genes that must be in gene set after filtering
#'
#' @return
#' @export
#'
#' @examples
geneSetEnrichment = function( geneSets, expressionVals, genes, minSetSize=0 )
{
    #browser()
    geneSets.filtered <- lapply(geneSets, intersect, y = genes)
    geneSets.filtered.length <- sapply(geneSets.filtered, length)
    names(geneSets.filtered) <- names(geneSets)

    # Filter out gene sets that are too small
    smallSets <- which(geneSets.filtered.length < minSetSize)
    if (length(smallSets) > 0) {
        warning(length(smallSets), ' sets are smaller than minSetSize')
        if (all(geneSets.filtered.length < minSetSize)) {
            return(NULL)
        } else {
            geneSets.filtered <- geneSets.filtered[ -smallSets ]
        }
    }

    # Order expression values
    o <- order(expressionVals, decreasing=TRUE)
    expressionVals <- expressionVals[o]
    names(expressionVals) <- genes[o]

    #
    resultList <- lapply(names(geneSets.filtered),
                         function(setName)
                             {
                                geneSet <- geneSets.filtered[[ setName ]]
                                n <- length(expressionVals) # number of genes in expression data
                                nSet <- length(geneSet) # length of gene set

                                inSetScore <- 0.5 / nSet # weight for genes in the gene set
                                outSetScore <- -0.5 / (n - nSet) # weight for genes not in the gene set

                                # Define weights
                                scoreVector <- as.numeric(names(expressionVals) %in% geneSet)
                                scoreVector[ scoreVector == 1 ] <- inSetScore
                                scoreVector[ scoreVector == 0 ] <- outSetScore

                                cumScore <- cumsum(scoreVector) # the running enrichment score


                                # find where expression data cross the origin (separates up and down genes)
                                zeroPoint <- which.min(abs(expressionVals))

                                # generate a normal kernel with max at zero (most upregulated) and min at zeroPoint
                                positive.kernel <- c( dnorm(seq(0,4, length=zeroPoint), mean=0, sd=1),
                                                      rep(0.0, n - zeroPoint) ) # zero after zero point

                                # generate a normal kernel with min at zeroPoint and max at the end (most downregulated)
                                negative.kernel <- c( rep(0.0, zeroPoint), # zero before zero point
                                                      rev(dnorm(seq(0,4, length=n-zeroPoint), mean=0, sd=1)))

                                # the best potential patterns, where all the set genes are maximally upregulated (first in scoreVector)
                                # or maximally downregulated (last in scoreVector)
                                bestPos <- c(rep(inSetScore, nSet), rep(outSetScore, (n-nSet)))
                                bestNeg <- c(rep(outSetScore, n-nSet), rep(inSetScore, nSet))

                                # apply the kernel and cumsum to get the best possible scores
                                maxPos <- sum( cumsum(bestPos) * positive.kernel )
                                maxNeg <- sum( cumsum(bestNeg) * negative.kernel )

                                # the actual scores as a fraction of the best scores
                                posScore <- sum( cumScore * positive.kernel ) / maxPos
                                negScore <- sum( cumScore * negative.kernel ) / maxNeg

                                gseaScores =
                                    tibble(signature = setName,
                                              genePosition = seq_len(n),
                                              enrichmentScore = cumScore,
                                              gene = names(expressionVals),
                                              isMax = FALSE,
                                              posScore = posScore,
                                              negScore = negScore )

                                if (max(gseaScores$enrichmentScore) > -min(gseaScores$enrichmentScore)) {
                                    gseaScores$isMax = gseaScores$enrichmentScore == max(gseaScores$enrichmentScore)
                                } else {
                                    gseaScores$isMax = gseaScores$enrichmentScore == min(gseaScores$enrichmentScore)
                                }

                                return(gseaScores)
                            }
    )

    goodResults <- lapply(resultList, function(x) !is.null(x))
    if (length(goodResults) == 0) return(NULL)
    return(do.call(rbind,resultList))

}


#' Plot the results of GSEA
#'
#' @param enrichmentScores output of geneSetEnrichment
#' @param plot.title optional title for the output plot
#'
#' @return
#' @export
#'
#' @examples
plotGeneSetEnrichment <- function( enrichmentScores, plot.title="Gene Set Enrichment Plot" ) {


    p <- ggplot(enrichmentScores, aes(x=genePosition, y=enrichmentScore, group=signature, color=signature, fill=signature)) +
        geom_hline( yintercept=0, linetype="dashed") +
        geom_line() +
        geom_label_repel(
            data=subset(enrichmentScores, isMax == TRUE), aes(label=signature),
            color="white", segment.color="#000000", segment.size=1,
            fontface = 'bold', size=2.5, box.padding = unit(0.25, "lines"),
            point.padding = unit(0.5, "lines")) +
        ggtitle( plot.title )

    print(p)
}


#' Return tibble with max enrichment score for each gene set
#'
#' @param enrichmentScores output of geneSetEnrichment
#'
#' @return
#' @export
#'
#' @examples
getMaxScores = function(enrichmentScores)
{
    enrichmentScores %>% filter(isMax)
}

# load("~/Documents/science/msigdb.rdata")
#
# # enrichment <- geneSetEnrichment( list(IFN_I=IFN_I, IFN_II=IFN_II), TCGAwICGC.Z[1,], colnames(TCGAwICGC.Z) )
# # plotGeneSetEnrichment(enrichment)
# # enrichment[ , lapply(.SD, function(x) x[1]), .SDcols=c("posScore", "negScore"), by=signature]
# #
#
# hpk1_data <- read_tsv("/Users/gene/Downloads/HPK1_corr_sort_PBMC.txt", col_names=FALSE)
# colnames(hpk1_data) <- c("Gene","Correlation")
#
# pos.list <- hpk1_data %>% arrange(-Correlation) %>% slice(1:10) %>% pull(Gene)
# gse.out <- geneSetEnrichment( list(POS=pos.list), hpk1_data$Correlation, hpk1_data$Gene )
# gse.out[ , lapply(.SD, function(x) x[1]), .SDcols=c("posScore", "negScore"), by=signature]
#
# # create a named vector of gene set lists
# gene.set.lists <- msigdb$Genes
# names(gene.set.lists) <- msigdb$Pathway
# # run analysis in parallel
# enrichment <- mclapply(seq_len(nrow(msigdb)), function(i) {geneSetEnrichment( gene.set.lists[i], hpk1_data$Correlation, hpk1_data$Gene, minSetSize=10)},  mc.cores=12)
# # filter out null results (should just be gene sets that are too small)
# names(enrichment) <- msigdb$Pathway
# enrichment <- enrichment[ !sapply(enrichment, is.null) ]
# # create a summary table of the results
# enrichment.summary <- lapply(enrichment, function(x) x[1, .(signature, posScore, negScore)]) %>% bind_rows
# # sort by highest pos and neg scores
# enrichment.summary[order(-posScore)][1:10]
# enrichment.summary[order(-negScore)][1:10]
