#' Update one effect from the vector of effects
#'
#' @param effects the current value of the effects vector, length at least 2
#' @param effects_index an integer from 1 to length of current effects vector, indicating which entry to update
#' @param distance parameter for proposal, ie, max distance from current effect value
#' @param trait univariate trait for mapping
#' @param residual_variance current value of residual variance.
#' @param collapsed_genotypes genoprobs matrix, collapsed per current configuration matrix
#' @return a list of length two

update_effect <- function(effects,
                          effect_index, # which element of effects is considered for update
                          distance = 1,
                          trait,
                          residual_variance,
                          collapsed_genotypes){
  proposal <- effects
  proposal[effect_index] <- runif(n = 1,
                                  min = effects[effect_index] - distance,
                                  max = effects[effect_index] + distance
  )
  numerator <- dnorm(trait,
                     mean = collapsed_genotypes %*% proposal,
                     sd = sqrt(residual_variance),
                     log = TRUE
  ) +
    dnorm(proposal[effect_index],
          mean = mean(trait), # check this!!
          sd = sqrt(var(trait)), # and this!!
          log = TRUE
    )
  denominator <- dnorm(trait,
                       mean = collapsed_genotypes %*% effects,
                       sd = sqrt(residual_variance),
                       log = TRUE
  ) +
    dnorm(effects[effect_index],
          mean = mean(trait), # check this!
          sd = sqrt(var(trait)),
          log = TRUE
    )
  output <- get_output(numerator = numerator,
                       denominator = denominator,
                       current = effects,
                       proposal = proposal
  )
  return(output)
}


#' Update residual variance for Jannink & Wu 2003 MCMC
#'
#' @param residual_variance current value of residual variance
#' @param distance max distance between current value and proposed value of residual variance
#' @param trait a univariate trait for QTL mapping
#' @param collapsed_genotypes genoprobs matrix collapsed per current value of configuration matrix
#' @param effects current value of effects vector.
#' @references Jannink & Wu (2003). \url{https://www.ncbi.nlm.nih.gov/pubmed/12872915}
update_residual_variance <- function(residual_variance,
                                     distance = 0.2,
                                     trait,
                                     collapsed_genotypes,
                                     effects){
  proposal <- runif(n = 1,
                    min = max(0, residual_variance - distance),
                    max = residual_variance + distance
  )
  numerator <- dnorm(trait,
                     mean = collapsed_genotypes %*% effects,
                     sd = sqrt(proposal),
                     log = TRUE
  ) + log (1 / min(2 * distance, distance + proposal))
  denominator <- dnorm(trait,
                       mean = collapsed_genotypes %*% effects,
                       sd = sqrt(residual_variance),
                       log = TRUE
  ) + log(1 / min(2 * distance, distance + residual_variance))
  output <- get_output(numerator = numerator,
                       denominator = denominator,
                       current = residual_variance,
                       proposal = proposal
  )
  return(output)
}

#' Update allelic number (and, of course, configuration, too) in Jannink and Wu 2003 MCMC.
#'
#' @param configuration an allelic series matrix, ie, binary matrix of 0s and 1s mapping founders to alleles at a single QTL
#' @param genoprobs genotype probabilities matrix for a single marker, ie, n_subjects by n_founders at 1 marker
#' @param effects vector containing current value of allelic effects, one entry per allele
#' @param trait a univariate trait matrix, n by 1
#' @param residual_variance residual variance
#' @param prior takes values "poisson" or "uniform" to specify which prior to use.
#' @param poisson_prior_mean Lambda, mean value, for the untruncated Poisson distribution
#' @details If a poisson prior is used for the allelic number, it is truncated to 2, 3, ..., number of founders (8 in the case of Diversity Outbred mice).

update_allelic_number <- function(configuration,
                                  genoprobs,
                                  effects,
                                  trait,
                                  residual_variance,
                                  prior,
                                  poisson_prior_mean = 2){
  allelic_number <- ncol(configuration)
  founder_number <- nrow(configuration)
  # get proposal
  if ((allelic_number < founder_number) & (allelic_number > 2)){
    increase_indicator <- as.logical(rbinom(n = 1,
                                            size = 1,
                                            prob = 0.5
    )
    )
  }
  if (allelic_number == 2){
    increase_indicator <- TRUE
  }
  if (allelic_number == founder_number){
    increase_indicator <- FALSE
  }
  if (increase_indicator){
    cp_out <- calc_prob_allelic_number_increase(configuration,
                                                genoprobs,
                                                effects,
                                                trait,
                                                residual_variance,
                                                prior,
                                                poisson_prior_mean
    )
  } else { # propose a decrease in allelic number
    cp_out <- calc_prob_allelic_number_decrease(configuration,
                                                genoprobs,
                                                effects,
                                                trait,
                                                residual_variance,
                                                prior,
                                                poisson_prior_mean
    )
  }
  output <- get_output(numerator = cp_out$numerator,
                       denominator = cp_out$denominator,
                       current = list(effects = effects, configuration = configuration),
                       proposal = list(effects = cp_out$effects_proposal, configuration = cp_out$configuration_proposal)
  )
  return(output)
}



# update configuration for a fixed number of alleles
update_configuration <- function(configuration,
                                 trait,
                                 genoprobs,
                                 effects,
                                 residual_variance){
  if (ncol(configuration) < nrow(configuration)){
    # determine which alleles have multiple founders
    alleles_with_multiple_founders <- which(colSums(configuration) > 1)
    # randomly choose an allele with multiple founders
    if (length(alleles_with_multiple_founders) > 1) {
      allele_to_change <- sample(alleles_with_multiple_founders, size = 1)
    } else {
      allele_to_change <- alleles_with_multiple_founders
    }
    # choose a founder with 'allele_to_change'
    founders_with_allele_to_change <- which(configuration[, allele_to_change] == 1)
    # randomly choose one from 'founders_with_allele_to_change'
    founder_to_change <- sample(founders_with_allele_to_change, size = 1)
    #
    v_possible_alleles <- (1:ncol(configuration))[- allele_to_change]
    if (length(v_possible_alleles) > 1){
      v_allele <- sample(v_possible_alleles, size = 1)
    } else {
      v_allele <- v_possible_alleles
    }
    configuration_proposal <- propose_configuration(configuration,
                                                    founder_to_change = founder_to_change,
                                                    allele_to_change = allele_to_change,
                                                    v_allele = v_allele
    )
    eta_u <- colSums(configuration)[allele_to_change]
    eta_v <- colSums(configuration_proposal)[v_allele]
    l2 <- sum(colSums(configuration) > 1)
    l2_proposal <- sum(colSums(configuration_proposal) > 1)
    collapsed_genotypes <- genoprobs %*% configuration
    collapsed_genotypes_proposal <- genoprobs %*% configuration_proposal
    numerator <- log(eta_u * l2) + dnorm(trait,
                                         mean = collapsed_genotypes_proposal %*% effects,
                                         sd = sqrt(residual_variance),
                                         log = TRUE
    )
    denominator <- log(eta_v * l2_proposal) + dnorm(trait,
                                                    mean = collapsed_genotypes %*% effects,
                                                    sd = sqrt(residual_variance),
                                                    log = TRUE
    )
    output <- get_output(numerator = numerator,
                         denominator = denominator,
                         current = configuration,
                         proposal = configuration_proposal
    )
  } else {
    output <- list(out = configuration, bernoulli = 1, alpha = 1, proposal = configuration)
  }
  return(output)
}

#' Draw samples from a MCMC chain per Jannink & Wu 2003
#'
#' @param initial_values a named list of length 3, where names are 'configuration', 'effects', and 'residual_variance'.
#' @param niter number of iterations
#' @param genoprobs a genoprobs matrix for QTL
#' @param trait a univariate trait affected by the QTL
#' @param allelic_number_prior Takes either of two values, "poisson" or "uniform".
#' @param poisson_prior_mean lamba, ie, mean, for the untruncated poisson.
#' @return a list of length two, where each entry is itself a list of results from each iteration
#' @export

jannink_mcmc <- function(initial_values,
                         niter = 10000,
                         genoprobs,
                         trait,
                         allelic_number_prior = "poisson",
                         poisson_prior_mean = 2){
  outs <- list()
  current <- initial_values
  all_outs <- list()
  for (i in 1:niter){
    # update allelic number
    ua_out <- update_allelic_number(configuration = current$configuration,
                                    genoprobs = genoprobs,
                                    effects = current$effects,
                                    trait = trait,
                                    residual_variance = current$residual_variance,
                                    prior = allelic_number_prior,
                                    poisson_prior_mean = poisson_prior_mean
    )
    current$effects <- ua_out$out[[1]]
    current$configuration <- ua_out$out[[2]]
    # update configuration for fixed allelic number
    uc_out <- update_configuration(configuration = current$configuration,
                                   trait = trait,
                                   genoprobs = genoprobs,
                                   effects = current$effects,
                                   residual_variance = current$residual_variance
    )
    current$configuration <- uc_out$out
    # iterate over effects vector
    ue_outs <- list()
    for (k in seq_along(current$effects)){
      ue_outs[[k]] <- update_effect(effects = current$effects,
                                    effect_index = k,
                                    trait = trait,
                                    residual_variance = current$residual_variance,
                                    collapsed_genotypes = genoprobs %*% current$configuration
      )
      current$effects <- ue_outs[[k]]$out
    }
    ##
    urv_out <- update_residual_variance(residual_variance = current$residual_variance,
                                        trait = trait,
                                        collapsed_genotypes = genoprobs %*% current$configuration,
                                        effects = current$effects
    )
    current$residual_variance <- urv_out$out
    #current$residual_variance <- 1
    outs[[i]] <- current
    all_outs[[i]] <- list(ua_out = ua_out,
                          uc_out = uc_out,
                          ue_outs = ue_outs,
                          urv_out = urv_out
                          )
  }
  return(list(outs = outs, all_outs = all_outs))
}


#' @export
jannink_mcmc_rv_only <- function(initial_values,
                         niter = 10000,
                         genoprobs,
                         trait,
                         allelic_number_prior = "poisson",
                         poisson_prior_mean = 2){
  outs <- list()
  current <- initial_values
  all_outs <- list()
  for (i in 1:niter){
    # update allelic number
    ##
    urv_out <- update_residual_variance(residual_variance = current$residual_variance,
                                        trait = trait,
                                        collapsed_genotypes = genoprobs %*% current$configuration,
                                        effects = current$effects
    )
    current$residual_variance <- urv_out$out
    #current$residual_variance <- 1
    outs[[i]] <- current
    all_outs[[i]] <- list(urv_out = urv_out
    )
  }
  return(list(outs = outs, all_outs = all_outs))
}


#' @export

jannink_mcmc_rv_fx_only <- function(initial_values,
                         niter = 10000,
                         genoprobs,
                         trait,
                         allelic_number_prior = "poisson",
                         poisson_prior_mean = 2){
  outs <- list()
  current <- initial_values
  all_outs <- list()
  for (i in 1:niter){
    # iterate over effects vector
    ue_outs <- list()
    for (k in seq_along(current$effects)){
      ue_outs[[k]] <- update_effect(effects = current$effects,
                                    effect_index = k,
                                    trait = trait,
                                    residual_variance = current$residual_variance,
                                    collapsed_genotypes = genoprobs %*% current$configuration
      )
      current$effects <- ue_outs[[k]]$out
    }
    ##
    urv_out <- update_residual_variance(residual_variance = current$residual_variance,
                                        trait = trait,
                                        collapsed_genotypes = genoprobs %*% current$configuration,
                                        effects = current$effects
    )
    current$residual_variance <- urv_out$out
    #current$residual_variance <- 1
    outs[[i]] <- current
    all_outs[[i]] <- list(ue_outs = ue_outs,
                          urv_out = urv_out
    )
  }
  return(list(outs = outs, all_outs = all_outs))
}
