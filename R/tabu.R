# unified tabu search implementation (both optimized and by spec).

tabu.search = function(x, start, whitelist, blacklist, score, extra.args,
    max.iter, maxp, optimized, tabu, debug = FALSE) {
  # cache nodes' labels.
  nodes = names(x)
  # cache the number of nodes.
  n.nodes = length(nodes)
  # set the iteration counter.
  iter = 1
  # check whether the score is score-equivalent.
  score.equivalence = is.score.equivalent(score, nodes, extra.args)
  # check whether the score is decomposable.
  score.decomposability = is.score.decomposable(score, extra.args)
  # allocate the cache matrix.
  cache = matrix(0, nrow = n.nodes, ncol = n.nodes)
  # nodes to be updated (all of them in the first iteration).
  updated = seq_len(n.nodes) - 1L
  # allocate the tabu list.
  tabu.list = vector("list", tabu)
  # maximum number of iteration the algorithm can go on without
  # improving the best network score.
  max.loss.iter = tabu
  # set the counter for suc iterations.
  loss.iter = 0
  # keep track of the best score value.
  best.score = -Inf
  # set the reference score.
  reference.score = per.node.score(network = start, score = score,
                      targets = nodes, extra.args = extra.args, data = x)
  # convert the blacklist to an adjacency matrix for easy use.
  if (!is.null(blacklist))
    blmat = arcs2amat(blacklist, nodes)
  else
    blmat = matrix(0L, nrow = n.nodes, ncol = n.nodes)
  # convert the whitelist to an adjacency matrix for easy use.
  if (!is.null(whitelist))
    wlmat = arcs2amat(whitelist, nodes)
  else
    wlmat = matrix(0L, nrow = n.nodes, ncol = n.nodes)
  if (debug) {
    cat("----------------------------------------------------------------\n")
    cat("* starting from the following network:\n")
    print(start)
    cat("* current score:", sum(reference.score), "\n")
    cat("* whitelisted arcs are:\n")
    if (!is.null(whitelist)) print(whitelist)
    cat("* blacklisted arcs are:\n")
    if (!is.null(blacklist)) print(blacklist)
    # set the metadata of the network; othewise the debugging output is
    # confusing and not nearly as informative.
    start$learning$algo = "tabu"
    start$learning$ntests = 0
    start$learning$test = score
    start$learning$args = extra.args
    start$learning$optimized = optimized
  }#THEN
total_score <- 0 #AZAR
best_scores_list <- list() #AZAR
    # Create an empty list to store adjacency matrices #AZAR
  adjacency_matrices_list <- list() #AZAR
     best_scores_all <- list()  #AZAR
    best_params_list <- list()  # New list to store the number of parameters for each iteration #AZAR
  
repeat {
  current = as.integer((iter - 1) %% tabu)
  
  # keep the best network seen so far and its score value for the evaluation
  # of the stopping rule; but always create a "best network" in the first
  # iteration using the starting network.
  if ((sum(reference.score) == -Inf) && (best.score == -Inf) && (iter > 1)) {
    old.score = per.node.score(network = best.network, score = score,
                                targets = nodes, extra.args = extra.args, data = x)
    singular.old.nodes = (old.score == -Inf)
    singular.new.nodes = (reference.score == -Inf)
    if (all(singular.old.nodes == singular.new.nodes)) {
      delta = robust.score.difference(
        sum(reference.score[!singular.new.nodes]),
        sum(old.score[!singular.old.nodes]))
    } else if (sum(singular.new.nodes) > sum(singular.old.nodes)) {
      delta = -Inf
    } else if (sum(singular.new.nodes) < sum(singular.old.nodes)) {
      delta = +Inf
    } else {
      delta = robust.score.difference(sum(reference.score), best.score)
    }
    if (delta > 0) {
      best.network = start
      best.score = sum(reference.score)
    }
  } else if ((robust.score.difference(sum(reference.score), best.score) > 0) ||
               (iter == 1)) {
    best.network = start
    best.score = sum(reference.score)
  }
  
  if (debug)
    cat("* iteration", iter, "using element", current, "of the tabu list.\n")
  
  amat = arcs2amat(start$arcs, nodes) #AZAR
  adjacency_matrices_list[[iter]] <- amat #AZAR
  
  nparents = colSums(amat) #AZAR
  
  .Call(call_tabu_hash,
        amat = amat,
        nodes = nodes,
        tabu.list = tabu.list,
        current = current)
  
  .Call(call_score_cache_fill,
        nodes = nodes,
        data = x,
        network = start,
        score = score,
        extra = extra.args,
        reference = reference.score,
        equivalence = score.equivalence && optimized,
        decomposability = score.decomposability,
        updated = (if (optimized) updated else seq(length(nodes)) - 1L),
        amat = amat,
        cache = cache,
        blmat = blmat,
        debug = debug)
  
  to.be.added = arcs.to.be.added(amat = amat, nodes = nodes,
                                  blacklist = blmat, whitelist = NULL, nparents = nparents,
                                  maxp = maxp, arcs = FALSE)
  
  bestop = .Call(call_tabu_step,
                 amat = amat,
                 nodes = nodes,
                 added = to.be.added,
                 cache = cache,
                 reference = reference.score,
                 wlmat = wlmat,
                 blmat = blmat,
                 tabu.list = tabu.list,
                 current = current,
                 baseline = 0,
                 nparents = nparents,
                 maxp = maxp,
                 debug = debug)
  
  cur_score <- bestop$score #AZAR
  cat("* BIC.score in iteration", iter, ":", sum(reference.score), "\n") #AZAR
  
  if (bestop$op == FALSE) {
    if (loss.iter >= max.loss.iter) {
      start = best.network
      if (debug) {
        cat("----------------------------------------------------------------\n")
        cat("* maximum number of iterations without improvements reached, stopping.\n")
        cat("* best network ever seen is:\n")
        print(best.network)
      }
      break
    } else {
      loss.iter = loss.iter + 1
    }
    
    if (debug) {
      cat("----------------------------------------------------------------\n")
      cat("* network score did not increase (for", loss.iter,
          "times), looking for a minimal decrease :\n")
    }
    
    bestop = .Call(call_tabu_step,
                   amat = amat,
                   nodes = nodes,
                   added = to.be.added,
                   cache = cache,
                   reference = reference.score,
                   wlmat = wlmat,
                   blmat = blmat,
                   tabu.list = tabu.list,
                   current = current,
                   baseline = -Inf,
                   nparents = nparents,
                   maxp = maxp,
                   debug = debug)
    
    if (bestop$op == FALSE) {
      if (debug) {
        cat("----------------------------------------------------------------\n")
        cat("* no more possible operations.\n")
        cat("@ stopping at iteration", iter, ".\n")
      }
      if (loss.iter > 0)
        start = best.network
      break
    }
  } else {
    if (robust.score.difference(sum(reference.score), best.score) > 0)
      loss.iter = 0
  }
  
  start = arc.operations(start, from = bestop$from, to = bestop$to,
                         op = bestop$op, check.cycles = FALSE, check.illegal = FALSE,
                         update = TRUE, debug = FALSE)
  
  BIC_score = sum(reference.score) #AZAR
  
  if (iter == 1) {  
    num_parameters = 0  
  } else {
    num_parameters = best_params_list[[length(best_params_list)]]
  } #AZAR
  
  log_likelihood = BIC_score + (num_parameters / 2) * log(n) #AZAR
  total_score <- total_score + log_likelihood  #AZAR
  cat(sprintf("Log-Likelihood in iteration %s: %s\n", iter, log_likelihood)) #AZAR
  best_scores_all[[length(best_scores_all) + 1]] <- sum(reference.score) #AZAR
  params <- nparams.backend(x = start, data = x, debug = debug) #AZAR
  best_params_list[[length(best_params_list) + 1]] <- params #AZAR
   weighted_matrix <- bestop$weights * amat #AZAR
  best_scores_list <- c(best_scores_list, log_likelihood ) #AZAR
    if (bestop$op == "reverse")
    updated = which(nodes %in% c(bestop$from, bestop$to)) - 1L
  else
    updated = which(nodes %in% bestop$to) - 1L
  
  if (debug) {
    start$learning$ntests = test.counter()
    cat("----------------------------------------------------------------\n")
    cat("* best operation was: ")
    if (bestop$op == "set")
      cat("adding", bestop$from, "->", bestop$to, ".\n")
    else if (bestop$op == "drop")
      cat("removing", bestop$from, "->", bestop$to, ".\n")
    else
      cat("reversing", bestop$from, "->", bestop$to, ".\n")
    cat("* current network is :\n") 
    print(start)
    cat(sprintf("* best score up to now: %s (delta: %s)\n", 
                format(best.score),
                format(robust.score.difference(sum(reference.score), best.score)))) #AZAR
  }
  
  if (iter >= max.iter) {
    if (debug)
      cat("@ stopping at iteration", max.iter, ".\n")
    
    if (loss.iter > 0)
      start = best.network
    
    break
  } else iter = iter + 1
}

cat("Total Summation of likelihood Scores after", iter, "iterations:", total_score, "\n") #AZAR
adjusted_scores <- as.numeric(unlist(best_scores_list)) / total_score #AZAR
cat("Total Summation of Scores after", iter , "iterations:", total_score, "\n")    #AZAR

cat("Best Scores List (at every 10 iterations):\n") #AZAR
print(best_scores_all) #AZAR

multiplied_scores <- lapply(1:length(adjusted_scores), function(i) {
  adjusted_scores[[i]] * adjacency_matrices_list[[i]]
}) #AZAR

print(adjacency_matrices_list) #AZAR

cat("List of Multiplied Scores:\n") #AZAR
print(multiplied_scores) #AZAR

final_matrix <- Reduce(`+`, multiplied_scores) #AZAR
cat("Final Matrix (sum of multiplied scores):\n") #AZAR
print(final_matrix) #AZAR

final_symmetric_matrix = final_matrix + t(final_matrix) #AZAR
cat("Final Symmetric Matrix:\n") #AZAR
print(final_symmetric_matrix) #AZAR


        
final_graph <- graphviz.plot(start)  #AZAR

return(list(adjacency_matrices_list = adjacency_matrices_list, 
            best_scores_list = best_scores_list, 
            best_params_list = best_params_list,
            multiplied_scores = multiplied_scores,
            final_matrix = final_matrix,
            final_network = start,
            final_graph = final_graph,
           final_symmetric_matrix = final_symmetric_matrix)) #AZAR
          
best_scores_df <- data.frame(iteration = seq_along(best_scores_list), score = best_scores_list) #AZAR
#write.csv(best_scores_df, file = "C:/Azar_Drive/relationships-between-variables1/01_preprocessing/best_scores.csv", row.names = FALSE) #AZAR
cat("Best scores list saved to 'best_scores.csv'\n") #AZAR
#write.csv(final_symmetric_matrix, file = "C:/Azar_Drive/relationships-between-variables1/01_preprocessing/final_symmetric_matrix.csv", row.names = FALSE) #AZAR
cat("Best scores list saved to 'final_symmetric_matrix.csv'\n") #AZAR

        
  #return(start)
}#TABU.SEARCH

