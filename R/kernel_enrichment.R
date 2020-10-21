# geneSetEnrichment: calculate enrichment scores for enriched and depleted gene sets
# geneSets: list of gene name vectors
# expressionVals: vector of gene scores such as differential expression values
# genes: vector of gene names in the same order as expressionVals
# minSetSize:  the smallest gene set to analyze (based on useable genes)
# parallel: use `mclapply` function
kernelEnrichment <- function( geneSets, expressionVals, genes, minSetSize=0, parallel=FALSE ) {
#	require(data.table)
	if (parallel) {
		require(parallel)
		which_lapply <- mclapply
	} else {
		which_lapply <- lapply
	}

	results <- tibble( i=seq_along(geneSets), geneSet = geneSets )
	results <- results %>% mutate(
					signature.genes = map(geneSet, ~intersect(.x, y=genes)),
					signature.length = map_dbl(signature.genes, length),
					signature.name = names(geneSets)
				)
	if (sum(results$signature.length < minSetSize) > 0 ) warning(sum(results$signature.length < minSetSize), ' sets are empty')
	if (sum(results$signature.length >= minSetSize) == 0 ) return(NULL)

	# group by gene set size to improve efficiency
	results <- results %>% dplyr::filter(signature.length >= minSetSize)
	unique.ns <- sort(unique(results$signature.length))

	o <- order(expressionVals, decreasing=TRUE)
	expressionVals <- expressionVals[o]
	genes <- genes[o]
#	names(expressionVals) <- genes[o]
	n <- length(expressionVals)
	

	# calculate scoring for signatures batched by signature length
	scoring <- which_lapply(unique.ns, function(nSet) {
		here <- results %>% dplyr::filter(signature.length == nSet) %>% pull(i)
		nNotSet <- n - nSet
		# where does the expression data cross the origin?  this separates up and down genes
		zeroPoint <- which.min(abs(expressionVals))
		left <- seq_len(n) <= zeroPoint
		right <- seq_len(n) > zeroPoint
		left.n <- sum(left)
		right.n <- sum(right)
		do.right <- right.n > 0
		# set scoring: positive for in set, negative for out of set
		inSetScore <-  1 / nSet # accounts for all hits on the positive or the negative side
		outSetScoreL <- -1 / left.n
		if (do.right) outSetScoreR <- -1 / right.n

		# generate a normal kernel with max at zero (most upregulated) and min at zeroPoint
		positive.kernel <- c( dnorm(seq(0, 4, length=zeroPoint), mean=0, sd=1), rep(0.0, n - zeroPoint) )
		# generate a normal kernel with min at zeroPoint and max at the end (most downregulated)
		if (do.right) negative.kernel <- c( rep(0.0, zeroPoint), rev(dnorm(seq(0,4, length=n-zeroPoint), mean=0, sd=1)))
		# the best potential patterns, where all the set genes are maximally upregulated (first in scoreVector)
		# or maximally downregulated (last in scoreVector)
		bestPos <- c(rep(inSetScore,  nSet),    rep(outSetScoreL, nNotSet))
		if (do.right) bestNeg <- c(rep(outSetScoreR, nNotSet), rep(inSetScore,  nSet))
		# apply the kernel and cumsum to get the best and worst possible scores
		maxPos <- sum( cumsum(bestPos) * positive.kernel ) / left.n
		if (do.right) maxNeg <- sum( rev(cumsum(rev(bestNeg))) * negative.kernel ) / right.n
		minPos <- sum(cumsum(rep(outSetScoreL, n))  * positive.kernel) / left.n
		if (do.right) minNeg <- sum(rev(cumsum(rep(outSetScoreR, n)))  * negative.kernel) / right.n
		
		#########################################
		# calculate empirical score distributions
		scoreVector.zeroed <- c(rep(outSetScoreL, left.n), rep(outSetScoreR, right.n))
		# randomize 1000 times
		randomized <- which_lapply(seq_len(1000), function(j) {
			scoreVector <- copy(scoreVector.zeroed)
			scoreVector[ sample(seq_len(n), size=nSet, replace=FALSE) ] <- inSetScore
			posScore <- sum( cumsum(scoreVector) * positive.kernel ) / left.n
			negScore <- NULL
			if (do.right) negScore <- sum( rev(cumsum(rev(scoreVector))) * negative.kernel ) / right.n 
			return(list(pos=posScore, neg=negScore))
		})
		randomized.pos <- sapply(randomized, function(x) x$pos)
		randomized.pos.mean <- mean(randomized.pos)
		randomized.pos.sd <- sd(randomized.pos)
		rm(randomized.pos)

		if (do.right) {
			randomized.neg <- sapply(randomized, function(x) x$neg)
			randomized.neg.mean <- mean(randomized.neg)
			randomized.neg.sd <- sd(randomized.neg)
			rm(randomized.neg)
		}
		rm(randomized)


		#########################################
		# calcaulte scores for each signature of this length
		results.here <- which_lapply(here, function(j) {
			hits <- genes %in% results$geneSet[[j]]
			scoreVector <- copy(scoreVector.zeroed)
			scoreVector[ hits ] <- inSetScore

			# add up the score from the left and from the right
			posScore <- sum( cumsum(scoreVector) * positive.kernel ) / left.n
			if (do.right) negScore <- sum( rev(cumsum(rev(scoreVector))) * negative.kernel ) / right.n 

			# make tibbles for reporting scores and details
			# calculate p-values based on background distribution
			result_scores <- tibble( i = j, 
									posScore = posScore,  
									posScoreStandardized=(posScore - minPos)/(maxPos-minPos),
									posPval = pnorm(posScore, mean=randomized.pos.mean, sd=randomized.pos.sd, lower.tail=FALSE)
				 				)
			if (do.right) {
				result_scores <- result_score %>% mutate(
						negScore = negScore,
						negScoreStandardized=(negScore - minNeg)/(maxNeg-minNeg),
						negPval = pnorm(negScore, mean=randomized.neg.mean, sd=randomized.neg.sd, lower.tail=FALSE)
				)
			}


			result_detail <- tibble( i = j,
									 X = seq_len(n),
									 hitPosition = hits, 
									 isMax = seq_len(n) == which.max(abs(result_detail$enrichmentScore)),
									 zeroPoint = seq_len(n) == zeroPoint,
									 enrichmentScore =  if_else(do.right,
									 	c((cumsum(scoreVector) * positive.kernel)[left], -1 * (rev(cumsum(rev(scoreVector))) * negative.kernel)[right]),
									 	(cumsum(scoreVector) * positive.kernel)
									 )
									)

			# return the combined results
			return(list(results=result_scores, details=result_detail))
		})
		return(results.here)
	})

	## unnest two levels of nesting
	# get rid of any null results, then double-unlist
	null.scoring <- sapply(scoring, is.null)
	if (all(null.scoring)) return(NULL)
	if (any(null.scoring)) {
		compiled_scoring <- scoring[!null.scoring] %>% unlist(recursive=FALSE) %>% unlist(recursive=FALSE)
	} else {
		compiled_scoring <- scoring %>% unlist(recursive=FALSE) %>% unlist(recursive=FALSE)
	}
	rm(null.scoring)
	# bind rows and join into original `results` tibble	
	are.results <- which(names(compiled_scoring) == "results")
	are.details <- which(names(compiled_scoring) == "details")

	scoring_results <- bind_rows(compiled_scoring[are.results])
	scoring_details <- bind_rows(compiled_scoring[are.details]) %>% nest(details=c(X, hitPosition, isMax, enrichmentScore))
	results <- results %>% left_join(scoring_results, by="i") %>% left_join(scoring_details, by="i")

	return(results)
}


# enrichmentScoreDT is the data.table output of geneSetEnrichment
# plot.title is the optional title for the output plot
kernelEnrichmentPlotter <- function( data, group, plot.title="Gene Set Enrichment Plot" ) {
	require(ggplot2)
	require(ggrepel)

	p <- ggplot(data, aes(x=X, y=enrichmentScore, group={{group}}, color={{group}})) +
		geom_hline( yintercept=0, linetype="dashed") +
		geom_line() +
		geom_label_repel(
			data=subset(data, isMax == TRUE), aes(label={{group}}),
				color="white", segment.color="#000000", segment.size=1,
				fontface = 'bold', size=2.5, box.padding = unit(0.25, "lines"),
				point.padding = unit(0.5, "lines")) +
		ggtitle( plot.title )
	print(p)
}



