get_output <- function(numerator, denominator, current, proposal){
  alpha <- min(0, sum(numerator - denominator))
  bernoulli <- rbinom(n = 1, size = 1, prob = exp(alpha))
  if (bernoulli) {out <- proposal} else {out <- current}
  return(list(out = out, bernoulli = bernoulli, alpha = alpha, proposal = proposal, current = current))
}

enumerate_configuration <- function(configuration){
  return(apply(X = configuration, MARGIN = 1, FUN = function(x)which(x == 1)))
}
#' Use configuration and effects vector to get each founder's effect, even when there are fewer alleles than founders
expand_effects <- function(effects, enumerated_configuration){
  return(effects[enumerated_configuration])
}


calc_prob_allelic_number_decrease <- function(configuration,
                                              genoprobs,
                                              effects,
                                              trait,
                                              residual_variance){
  # propose a decrease in allelic number
  allelic_number <- ncol(configuration)
  founder_number <- nrow(configuration)
  # choose two alleles to combine into one
  alleles_to_combine <- sample(1:allelic_number, size = 2) # vector of length 2
  # make configuration proposal matrix
  configuration_proposal <- calc_configuration_combining_alleles(configuration, alleles_to_combine)
  verify_configuration(configuration_proposal)
  # make effects proposal vector
  effects_proposal <- effects[ - max(alleles_to_combine)]
  # define collapsed genotypes
  collapsed_genotypes_proposal <- genoprobs %*% configuration_proposal
  collapsed_genotypes <- genoprobs %*% configuration
  numerator <- log(allelic_number) +
    log(allelic_number - 1) +
    log(count_configurations(founder_number = founder_number,
                             allelic_number = allelic_number
    )) +
    dnorm(trait,
          mean = collapsed_genotypes_proposal %*% effects_proposal,
          sd = sqrt(residual_variance),
          log = TRUE
    )
  # define etau plus etav
  eta_u_plus_v <- colSums(configuration_proposal)[min(alleles_to_combine)]
  l2_proposal <- sum(colSums(configuration_proposal) > 1)
  denominator <-log(2 * (2 ^ (eta_u_plus_v - 1) - 1)) +
    log(l2_proposal) +
    log(count_configurations(founder_number = founder_number,
                             allelic_number = allelic_number - 1
    )) +
    dnorm(trait,
          mean = collapsed_genotypes %*% effects,
          sd = sqrt(residual_variance),
          log = TRUE
    )
  return(list(numerator = numerator,
              denominator = denominator,
              configuration_proposal = configuration_proposal,
              effects_proposal = effects_proposal
  )
  )
}
##
calc_configuration_combining_alleles <- function(configuration,
                                                 alleles_to_combine){
  configuration[ , min(alleles_to_combine)] <- configuration[ , min(alleles_to_combine)] +
    configuration[ , max(alleles_to_combine)]
  out <- configuration[, - max(alleles_to_combine)]
  verify_configuration(out)
  return(out)
}
##
# make new configuration matrix
create_wider_configuration <- function(current_configuration){
  configuration_proposal <- current_configuration
  founder_number <- nrow(current_configuration)
  allelic_number <- ncol(current_configuration)
  binary <- rep(0, founder_number)
  shared_alleles <- which(colSums(current_configuration) > 1) # numeric vector of length up to 4 (for DO)
  if (length(shared_alleles) > 1) {
    allele_to_split <- sample(x = shared_alleles, size = 1)
  } else {
    allele_to_split <- shared_alleles
  }
  founders_to_split <- which(current_configuration[ , allele_to_split] == 1)
  group_indicators <- sample_partition_into_two(founders_to_split) # vector of 0s and 1s
  binary[founders_to_split[as.logical(group_indicators)]] <- 1
  configuration_proposal[ , allele_to_split] <- binary
  ## add a new column
  binary0 <- rep(0, founder_number)
  binary0[founders_to_split[as.logical(1 - group_indicators)]] <- 1
  configuration_proposal <- cbind(configuration_proposal, binary0)
  colnames(configuration_proposal) <- NULL
  verify_configuration(configuration_proposal)
  return(list(configuration_proposal = configuration_proposal, allele_to_split = allele_to_split))
}

sample_partition_into_two <- function(vector){
  all_partitions <- clue::as.cl_partition(vector)$.Data
  partitions_of_size_two <- all_partitions[, apply(FUN = max, X = all_partitions, MARGIN = 2) == 1]
  if (length(partitions_of_size_two) > 1){
    part_out_index <- sample(1:ncol(partitions_of_size_two), size = 1)
  } else {
    part_out_index <- 1
  }
  return(partitions_of_size_two[, part_out_index])
}

calc_prob_allelic_number_increase <- function(configuration,
                                              genoprobs,
                                              effects,
                                              trait,
                                              residual_variance){
  cwc <- create_wider_configuration(configuration)
  configuration_proposal <- cwc$configuration
  allele_to_split <- cwc$allele_to_split
  founder_number <- nrow(configuration)
  allelic_number <- ncol(configuration)
  #draw new effect for proposal
  new_effect <- rnorm(n = 1,
                      mean = mean(trait), # IS THIS RIGHT??? maybe try mean(trait)??
                      sd = sqrt(var(trait))
  ) # check this!!
  effects_proposal <- c(effects, new_effect)
  eta_u <- colSums(configuration)[allele_to_split]
  # how many founders have the allele to 'split'??
  # collapse genoprobs in each of two ways
  collapsed_genotypes <- genoprobs %*% configuration
  collapsed_genotypes_proposal <- genoprobs %*% configuration_proposal

  numerator <- log(2 * (2 ^ (eta_u - 1) - 1)) +
    log(sum(colSums(configuration) > 1)) +
    log(count_configurations(founder_number = founder_number,
                             allelic_number = allelic_number
    )) +
    dnorm(trait,
          mean = collapsed_genotypes_proposal %*% effects_proposal,
          sd = sqrt(residual_variance),
          log = TRUE
    )
  # how many alleles are carried by at least 2 founders?
  denominator <- log(allelic_number + 1) +
    log(allelic_number) +
    log(count_configurations(founder_number = founder_number,
                             allelic_number = 1 + allelic_number
    )) +
    dnorm(trait,
          mean = collapsed_genotypes %*% effects,
          sd = sqrt(residual_variance),
          log = TRUE
    )
  return(list(numerator = numerator,
              denominator = denominator,
              configuration_proposal = configuration_proposal,
              effects_proposal = effects_proposal
  )
  )
}


count_configurations <- function(founder_number = 8, allelic_number){
  # allele labels don't matter, but founder labels matter.
  # we need stirling numbers of the second kind
  return(as.integer(gmp::Stirling2(n = founder_number, k = allelic_number)))
}


# for update_configuration step
propose_configuration <- function(configuration, founder_to_change, allele_to_change, v_allele){
  configuration[founder_to_change, allele_to_change] <- 0
  configuration[founder_to_change, v_allele] <- 1
  verify_configuration(configuration)
  return(configuration)
}


verify_configuration <- function(configuration){
  if ((sum(configuration) != nrow(configuration)) |
      (!identical(rowSums(configuration), rep(1, nrow(configuration)))) |
      (sum(configuration == 0) != (ncol(configuration) - 1) * nrow(configuration)) |
      (prod(colSums(configuration)) == 0)){
    print(configuration)
    stop("invalid configuration")
  }
}
