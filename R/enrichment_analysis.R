#' Run GSEA given differential expression analysis and list of gene sets
#'
#' @param resObj ouput of running results on a DESEQ2 dds object
#' @param gene_sets list of gene sets to run enrichment. See read_enrichR
#' @param orthologs table containing mapping information between human and mouse gene symbols
#' @param symbols table containing mapping between mouse ensembl_ids and gene symbols
#'
#' @return
#' @export
#'
#' @examples
run_fgsea = function(resObj, gene_sets, orthologs = human_mouse_orthologs, symbols = gene_symbols) {
    # Format results for gsea
    orderedResults =
        resObj %>%
        as_tibble() %>%
        add_column(ensembl_gene_id = rownames(resObj)) %>%
        left_join(symbols, by = "ensembl_gene_id") %>%
        select(GeneSymbol = mgi_symbol, everything()) %>%
        #left_join(orthologs, by = c("GeneSymbol" = "mouseSymbol"))
        nest_join(human_mouse_orthologs, by = c("GeneSymbol" = "mouseSymbol")) %>%
        mutate(numSym = map(human_mouse_orthologs, function(.x) {nrow(.x)})) %>% # number of human orthologs for each mouse gene
        filter(numSym == 1 ) %>% # use only genes with unique orthologs
        unnest(human_mouse_orthologs) %>%
        select(-numSym) %>%
        filter(!is.na(GeneSymbol) & !is.na(humanSymbol)) %>%
        filter( !is.na(padj) & !is.na(log2FoldChange)) %>%
        mutate(rankStat = -log10(padj) * sign(log2FoldChange) ) %>% # ranking statistic used for GSEA
        filter(is.finite(rankStat)) %>%
        arrange(desc(rankStat) )
    #filter(GeneSymbol %in%  converted_symbols$mouseSymbol) %>% # use only symbols with a unique human ortholog
    #left_join(converted_symbols, by = c("GeneSymbol" = "mouseSymbol")

    rankedGenes = orderedResults %>% pull(rankStat)
    names(rankedGenes) = orderedResults %>% pull(humanSymbol)

    # Run GSEA
    fgseaResults =
        fgsea(pathways = gene_sets,
              stats = rankedGenes,
              #nperm = 10000,
              minSize = 5)


    # Collapse redundant pathways and plot running enrichment scores of top pathways
    collapsedPathways = collapsePathways( fgseaResults %>% arrange(padj) %>% filter(padj < 0.01),
                                          gene_sets, rankedGenes)

    mainPathways <- fgseaResults[pathway %in% collapsedPathways$mainPathways][
        order(-NES), pathway]

    topPathwaysUp <-
        fgseaResults %>% filter(pathway %in% mainPathways, ES > 0) %>% arrange(padj) %>% pull(pathway) %>% head(n=10)
    topPathwaysDown <-
        fgseaResults %>% filter(pathway %in% mainPathways, ES < 0) %>% arrange(padj) %>% pull(pathway) %>% head(n=10)
    topPathways <-
        c(topPathwaysUp, rev(topPathwaysDown))


    plt = plotGseaTable(gene_sets[topPathways], rankedGenes, fgseaResults, gseaParam=0.5, render = FALSE)
    plot(plt)



    return(list("enrichment" =
                    fgseaResults %>%
                    arrange(padj,desc(NES)),
                "input" = rankedGenes,
                "top_pathways" = topPathways))
}


#' Plot expression of leading edge genes given enrichment result
#'
#' @param enrichment_result result of running run_fgsea
#' @param pathway_name pathway to plot leading genes
#' @param dds DESEQ2 dds object used for running fgsea
#' @param orthologs table containing mapping information between human and mouse gene symbols
#' @param symbols table containing mapping between mouse ensembl_ids and gene symbols
#'
#' @return
#' @export
#'
#' @examples
plot_leading_edge = function(enrichment_result, pathway_name, dds,  orthologs = human_mouse_orthologs, symbols = gene_symbols) {

    human_leading_edge = enrichment_result$enrichment %>% filter(pathway == pathway_name) %>% pull(leadingEdge) %>% unlist() %>% as.character()
    mouse_leading_edge = orthologs %>% filter(humanSymbol %in% human_leading_edge ) %>%  pull(mouseSymbol)

    counts(dds) %>%
        as_tibble(rownames = "ensembl_gene_id") %>%
        left_join(symbols, by = "ensembl_gene_id") %>%
        select(-ensembl_gene_id) %>%
        select(GeneSymbol = mgi_symbol, everything()) %>%
        filter(GeneSymbol %in% mouse_leading_edge) %>%
        pivot_longer(cols = -GeneSymbol,names_to = "sample_id",values_to = "counts") %>%
        left_join(tibble(sample_id = dds$sample_id, condition = dds$condition), by = "sample_id") %>%
        ggplot(aes(x = condition, y = counts)) +
        geom_point() +
        geom_boxplot(outlier.shape = NA) + geom_jitter(size = 0.01) + #geom_text(aes(label = sample_id)) +
        facet_wrap(~ GeneSymbol, scales = "free_y") +
        theme_bw(base_size = 16) +
        scale_y_log10() +
        xlab("") + ylab("counts") + ggtitle(pathway_name) +
        theme(plot.title = element_text(hjust = 0.5),
              strip.background = element_blank(),
              panel.grid.minor = element_blank(),
              panel.grid.major.x = element_blank(),
              axis.text.x = element_text(angle = 0))
}

#' Import a collection of gene sets downloaded from enrichR
#'
#' @param file_path path to file containing gene sets
#'
#' @return
#' @export
#'
#' @examples
read_enrichR = function(file_path){
    gene_sets = read_tsv(file_path, col_names = FALSE) %>%
        transpose_tibble(row.names = "X1",col_label = "GeneSymbol") %>%
        select(-GeneSymbol) %>%
        as.list() %>%
        lapply( function(x){ x[!is.na(x)]})
}
