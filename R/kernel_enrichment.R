# kernelEnrichment: calculate enrichment scores for enriched and depleted gene sets
# geneSets: list of gene name vectors
# expressionVals: vector of gene scores such as differential expression values
# genes: vector of gene names in the same order as expressionVals
# minSetSize:  the smallest gene set to analyze (based on useable genes)
# parallel: use `mclapply` function
kernelEnrichment <- function( geneSets, expressionVals, genes, minSetSize=0, randomization_n=1000, parallel=FALSE ) {
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
	
	# where does the expression data cross the origin?  this separates up and down genes
	zeroPoint <- which.min(abs(expressionVals))
	left <- 1:zeroPoint
	right <- 1:(n-zeroPoint)
	right_rev <- n:(zeroPoint+1)
	left.n <- length(left)
	right.n <- length(right)
	do.right <- zeroPoint < n
	
	# calculate scoring for signatures batched by signature length
	scoring <- which_lapply(unique.ns, function(nSet) {
		here <- results %>% dplyr::filter(signature.length == nSet) %>% pull(i)
		nNotSet <- n - nSet
		# set scoring: positive for in set, negative for out of set
		inSetScore <-  1 / nSet # accounts for all hits on the positive or the negative side
		outSetScore <- -1 / nNotSet
		# outSetScoreL <- -1 / left.n
		# if (do.right) outSetScoreR <- -1 / right.n
		
		# generate a normal kernel with max at zero (most upregulated) and min at zeroPoint
		positive.kernel <- dnorm(seq(0, 4, length=zeroPoint), mean=0, sd=1)
		# generate a normal kernel with min at zeroPoint and max at the end (most downregulated)
		if (do.right) negative.kernel <- dnorm(seq(0,4, length=n-zeroPoint), mean=0, sd=1)
		# the best potential patterns, where all the set genes are maximally upregulated (first in scoreVector)
		# or maximally downregulated (last in scoreVector)
		optimalScores <- c(rep(inSetScore,  nSet), rep(outSetScore, nNotSet))  #bestPos <- c(rep(inSetScore,  nSet),    rep(outSetScoreL, nNotSet))
		# apply the kernel and cumsum to get the best and worst possible scores
		maxPos <- sum( cumsum(optimalScores)[left] * positive.kernel ) / left.n
		if (do.right) maxNeg <- sum( cumsum(optimalScores[right]) * negative.kernel ) / right.n
		minPos <- sum(cumsum(rep(outSetScore, left.n))  * positive.kernel) / left.n
		if (do.right) minNeg <- sum(cumsum(rep(outSetScore, right.n))  * negative.kernel) / right.n
		
		#########################################
		# calculate empirical score distributions
		scoreVector.zeroed <- rep(outSetScore, n) # scoreVector.zeroed <- c(rep(outSetScoreL, left.n), rep(outSetScoreR, right.n))
		# randomize 1000 times
		randoms <- lapply(seq_len(randomization_n), function(ri) {
			scoreVector <- scoreVector.zeroed
			scoreVector[ sample(n, size=nSet, replace=FALSE) ] <- inSetScore
			scoreVector
		})

		randomized <- which_lapply(seq_len(randomization_n), function(ri) {
			posScore <- sum( cumsum(randoms[[ri]][left]) * positive.kernel ) / left.n
			negScore <- NULL
			if (do.right) negScore <- sum( cumsum(randoms[[ri]][right]) * negative.kernel ) / right.n 
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
		# calculate scores for each signature of this length
		results.here <- which_lapply(here, function(j) {
			hits <- genes %in% unlist(results$geneSet[ results$i == j ])
			scoreVector <- copy(scoreVector.zeroed)
			scoreVector[ hits ] <- inSetScore
			
			# add up the score from the left and from the right
			posScore <- sum( cumsum(scoreVector[left]) * positive.kernel ) / left.n
			if (do.right) negScore <- sum( cumsum(scoreVector[right_rev]) * negative.kernel ) / right.n 
			
			# make tibbles for reporting scores and details
			# calculate p-values based on background distribution
			result_scores <- tibble( i = j, 
									 posScore = posScore,  
									 posScoreStandardized=(posScore - minPos)/(maxPos-minPos),
									 posPval = pnorm(posScore, mean=randomized.pos.mean, sd=randomized.pos.sd, lower.tail=FALSE)
			)
			if (do.right) {
				result_scores <- result_scores %>% mutate(
					negScore = negScore,
					negScoreStandardized=(negScore - minNeg)/(maxNeg-minNeg),
					negPval = pnorm(negScore, mean=randomized.neg.mean, sd=randomized.neg.sd, lower.tail=FALSE)
				)
			}
			
			
			result_detail <- tibble( i = j,
									 X = seq_len(n),
									 hitPosition = hits, 
									 zeroPoint = seq_len(n) == zeroPoint,
									 enrichmentScore =  #ifelse(do.right,  # this doesn't work with `if_else`
									 	#c(cumsum(scoreVector[left]), -1 * rev(cumsum(rev(scoreVector[right])))),
									 	cumsum(scoreVector)
									 #)
			) %>% mutate(isMax = seq_len(n) == which.max(abs(enrichmentScore)))
			
			
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
	scoring_details <- bind_rows(compiled_scoring[are.details]) %>% nest(details=-i)
	results <- results %>% left_join(scoring_results, by="i") %>% left_join(scoring_details, by="i")
	
	return(results)
}


# kernelEnrichment: calculate enrichment scores for enriched and depleted gene sets
# geneSets: list of gene name vectors
# expressionVals: matrix/data.table of gene scores such as differential expression values (rows=genes, cols=samples)
# genes: vector of gene names in the same order as expressionVals
# minSetSize:  the smallest gene set to analyze (based on useable genes)
# parallel: use `mclapply` function
kernelEnrichmentMulti <- function( geneSets, expressionVals, genes, minSetSize=0, randomization_n=1000, parallel=FALSE ) {
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
	
	n <- nrow(expressionVals)
	sample.n <- ncol(expressionVals)
	sns <- seq_len(sample.n)
	sample_names <- colnames(expressionVals)
	if (length(sample_names)==0) sample_names <- sns
	
	# where does the expression data cross the origin?  this separates up and down genes
	expressionOrders <- apply(expressionVals, 2, order, decreasing=TRUE)
	zeroPoints <- sapply(sns, function(j) which.min(abs(expressionVals[,j][expressionOrders[,j]])))
	lefts <- lapply(sns, function(j) which(seq_len(n) <= zeroPoints[j]))
	rights_rev <- lapply(sns, function(j) rev(which(seq_len(n) > zeroPoints[j])))
	left.ns <- sapply(lefts, length)
	right.ns <- sapply(rights_rev, length)
	do.rights <- sapply(right.ns, function(x) x > 0)
	
	# calculate scoring for signatures batched by signature length
	scoring <- which_lapply(unique.ns, function(nSet) {
		here <- results %>% dplyr::filter(signature.length == nSet) %>% pull(i)
		nNotSet <- n - nSet
		# set scoring: positive for in set, negative for out of set
		inSetScore <-  1 / nSet # accounts for all hits on the positive or the negative side
		outSetScore <- -1 / nNotSet
		
		# generate a normal kernel with max at zero (most upregulated) and min at zeroPoint
		positive.kernel.Ls <- lapply(zeroPoints, function(zp) dnorm(seq(0, 4, length=zp), mean=0, sd=1) )
		# generate a normal kernel with min at zeroPoint and max at the end (most downregulated)
		negative.kernel.Rs <- lapply(sns, function(i) dnorm(seq(0,4, length=n-zeroPoints[[i]]), mean=0, sd=1) )

		# the best potential patterns, where all the set genes are maximally upregulated (first in scoreVector)
		# or maximally downregulated (last in scoreVector)
		optimalScores <- c(rep(inSetScore,  nSet), rep(outSetScore, nNotSet))
		# apply the kernel and cumsum to get the best and worst possible scores
		maxPos <- sapply(positive.kernel.Ls, function(pk) sum( cumsum(optimalScores)[1:length(pk)] * pk ) ) / left.ns
		maxNeg <- sapply(negative.kernel.Rs, function(nk) sum( cumsum(optimalScores)[1:length(nk)] * nk )) / right.ns
		minPos <- sapply(positive.kernel.Ls, function(pk) sum( cumsum(rep(outSetScore, length(pk)))  * pk)) / left.ns
		minNeg <- sapply(negative.kernel.Rs, function(nk) sum( cumsum(rep(outSetScore, length(nk))) * nk)) / right.ns
		
		#########################################
		# calculate empirical score distributions
		scoreVector.zeroed <- rep(outSetScore, n) 
		# randomize 1000 times, for each column of expressionVals
		randoms <- lapply(seq_len(randomization_n), function(ri) {
			scoreVector <- scoreVector.zeroed
			scoreVector[ sample(n, size=nSet, replace=FALSE) ] <- inSetScore
			scoreVector
		})

		randomized <- lapply(sns, function(i) {
			left <- lefts[[i]]; right <- rights_rev[[i]]
			left.n <- left.ns[i]; right.n <- right.ns[i]
			positive.kernel.L <- positive.kernel.Ls[[i]]
			negative.kernel.R <- negative.kernel.Rs[[i]]
			randomized <- which_lapply(seq_len(randomization_n), function(ri) {
				sv <- randoms[[ri]]
				posScore <- sum( cumsum(sv[left]) * positive.kernel.L ) / left.n
				negScore <- sum( cumsum(sv[right]) * negative.kernel.R ) / right.n 
				return(list(pos=posScore, neg=negScore))
			})
			randomized.pos <- sapply(randomized, function(x) x$pos)
			randomized.pos.mean <- mean(randomized.pos)
			randomized.pos.sd <- sd(randomized.pos)
			if (do.rights[i]) {
				randomized.neg <- sapply(randomized, function(x) x$neg)
				randomized.neg.mean <- mean(randomized.neg)
				randomized.neg.sd <- sd(randomized.neg)
			} else {
				randomized.neg.mean <- NA
				randomized.neg.sd <- NA
			}
			return(c(pos.m=randomized.pos.mean, pos.sd=randomized.pos.sd, neg.m=randomized.neg.mean, neg.sd=randomized.neg.sd))
		})
		rm(randoms)
		
		#########################################
		# calculate scores for each signature of this length
		results.here <- which_lapply(here, function(j) {
			hits <- lapply(sns, function(i) genes[expressionOrders[,i]] %in% unlist(results$geneSet[ results$i == j ]))
			scoreVectors <- lapply(sns, function(i) {V=scoreVector.zeroed; V[hits[[i]]] <- inSetScore; V})
			
			# add up the score from the left and from the right
			posScores <- sapply(sns, function(i) sum( cumsum(scoreVectors[[i]][lefts[[i]]]) * positive.kernel.Ls[[i]] )) / left.ns
			posScoresStandardized <- (posScores - minPos)/(maxPos-minPos)
			posPvals <- sapply(sns, function(i) pnorm(posScores[i], mean=randomized[[i]][["pos.m"]], sd=randomized[[i]][["pos.sd"]], lower.tail=FALSE))
			
			negScores <- sapply(sns, function(i) sum( cumsum(scoreVectors[[i]][rights_rev[[i]]]) * negative.kernel.Rs[[i]] )) / right.ns
			negScoresStandardized <- (negScores - minNeg)/(maxNeg-minNeg)
			negPvals <- sapply(sns, function(i) if(do.rights[i]) pnorm(negScores[i], mean=randomized[[i]][["neg.m"]], sd=randomized[[i]][["neg.sd"]], lower.tail=FALSE) else NA_real_)
			
			# make tibbles for reporting scores and details
			# calculate p-values based on background distribution
			return(tibble( i = j,
						   sample_name = sample_names,
						   posScore = posScores, posScoreStandardized = posScoresStandardized, posPval = posPvals,
						   negScore = negScores, negScoreStandardized = negScoresStandardized, negPval = negPvals 
			))
		})
		return(results.here)
	})
	
	## unnest lists and combine
	null.scoring <- sapply(scoring, is.null)
	if (all(null.scoring)) return(NULL)
	compiled_scoring <- scoring[!null.scoring] %>% unlist(recursive=FALSE) %>% bind_rows
	results <- results %>% left_join(compiled_scoring, by="i")
	
	return(results)
}


# kernelEnrichmentPlotter: Generate enrichment score plots
# data: a data frame as found in the `details` section of the kernelEnrichment output
# group: the column in data to group (should either be the sample or the gene set column)
# plot.title: an optional title for the output plot
# label.lines: label each line on the plot instead of in the legend
kernelEnrichmentPlotter <- function( data, group, plot.title="Gene Set Enrichment Plot", label.lines=TRUE ) {
	require(ggplot2)
	require(ggrepel)
	
	zeroPoints <- data %>% filter(zeroPoint == TRUE) # point in the gene expression data that goes from positive to negative
	if (label.lines) label_data <- data %>% filter(isMax == TRUE) # positions of labels
	
	p <- ggplot(data, aes(x=X, y=enrichmentScore, group={{group}}, color={{group}})) +
		geom_hline( yintercept=0, linetype="dashed") +
		geom_vline( data=zeroPoints, aes(xintercept=X, color={{group}}), linetype="dashed", alpha=0.5) +
		geom_line( ) +
		ggtitle( plot.title )
	if (label.lines) p <- p +
		geom_label_repel(
			data=label_data, aes(label={{group}}),
			color="black", segment.color="#000000", segment.size=0.5,
			fontface = 'bold', size=2.5, box.padding = unit(0.25, "lines"),
			point.padding = unit(0.5, "lines")) +
		theme(legend.position="none")
	
	print(p)
}



