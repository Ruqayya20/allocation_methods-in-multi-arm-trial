## December/2021
## Updated/Jan 2023
## This code is written by Michael Grayling
## This is simulation codes to evaluate the performance measures for six different allocation methods. 
## We perform our investigations for two different settings based on real clinical trials. 

#install.packages("Rfast")
#install.packages("tibble")
set.seed(123)
evaluate_randomisation_routine <- function(eff_tr = rep(target_effect, k), # The treatment effect for each experimental
           target_effect = 0.006,
           n = ceiling(base_n*c(0.5, 1, 2)), # Total sample size
           base_n = 400,
           method = "simple", # Method for randomising - can be "simple",
           # "block", "stratified", "minimisation"
           prob_cov = rep(0.25, 4), # Probability of a patient having each
           # covariate - the length of this determines
           # the number of covariates
            eff_cov = 0.006*c(1, 0.5, 0.1, 0), # The effect on the outcome data of each
           # covariate
           num_cov_in_rand = 2, # The number of covariates used in the
           # randomisation routine (assumed to be the first
           # num_cov_in_rand covariates) - only used if
           # method is "stratified" or "minimisation"
           B = 6, # Length of the blocks - only used if method is
           # "block" or "stratified" or "block urn"
           p = c(0.7,0.9), # Only used when method is "minimisation"
           burn_in = ceiling(0.1*n), # Only used when method is "minimisation"
           w = 1, # Only used method method is "urn" - signifies number of balls initially in urn of each colour
           a = 1, # Only used method method is "urn"
           b = 2, # Only used method method is "urn"
           sd = 0.011,
           replicates = 10000) { # Number of replicate simulations

  treatment_assignment      <-
    function(data, K, n, method, B, p, J, all_cov_combs, numeric_n, seq_J, Kp1,
             seq_Kw0, seq_js, cov_wt, diag_ratio, num_matching_cov, imbalance,
             seq_nm1, seq_Kp1, assignment_probs, J_rand, seq_J_rand, burn_in,
             w, a, b) {
      if (method == "simple") {
        treatment                          <- sample.int(Kp1, n,
                                                         replace = TRUE) - 1L
      } else if (method == "block") {
        treatment                          <- numeric_n
        each_block                         <- rep(seq_Kw0, B/Kp1)
        for (i in 1:floor(n/B)) {
          treatment[(B*(i - 1) + 1):(B*i)] <- sample(each_block, B)
        }
        if ((n/B)%%1 != 0) {
          treatment[(B*floor(n/B) + 1):n]  <- sample(each_block,
                                                     n - B*floor(n/B))
        }
      } else if (method == "stratified") {
        treatment                          <- numeric_n
        num_strata                         <- 2^J_rand
        each_block                         <- rep(seq_Kw0, B/Kp1)
        blocks                             <- matrix(0, n, num_strata)
        for (i in 1:(n/B)) {
          range_rows                       <- (B*(i - 1) + 1):(B*i)
          for (j in 1:num_strata) {
            blocks[range_rows, j]          <- sample(each_block, B)
          }
        }
        for (i in 1:num_strata) {
          which_strata_i                   <- which(data$strata == i)
          treatment[which_strata_i]        <- blocks[1:length(which_strata_i), i]
        }
      } else if (method == "minimisation") {
        treatment                          <- numeric_n
        cov_mat                            <- as.matrix(data[, seq_J_rand])
        treatment[1:burn_in]               <- sample.int(Kp1, burn_in,
                                                         replace = TRUE) - 1L
        for (j in (burn_in + 1):n) {
          matching_cov                     <-
            (cov_mat[seq_js[[j]], , drop = FALSE] ==
               matrix(cov_mat[j, ], j - 1, J_rand, byrow = TRUE))
          for (k in seq_Kp1) {
            num_matching_cov[, k]          <-
              Rfast::colsums(matching_cov[treatment[seq_js[[j]]] == k - 1, ,
                                          drop = FALSE])
          }
          for (k in seq_Kp1) {
            temp                           <- num_matching_cov
            temp[, k]                      <- temp[, k] + 1
            range_level                    <- apply(temp%*%diag_ratio, 1, range)
            imbalance[k]                   <- sum(cov_wt*(range_level[2, ] -
                                                            range_level[1, ]))
          }
          trt_mini                         <- seq_Kw0[imbalance == min(imbalance)]
          len_trt_mini                     <- length(trt_mini)
          if (len_trt_mini == Kp1) {
            treatment[j]                   <- sample.int(Kp1, 1) - 1L
          } else {
            treatment[j]                   <-
              sample(c(trt_mini, seq_Kw0[-(trt_mini + 1)]), 1,
                     prob = assignment_probs[[len_trt_mini]])
          }
        }
      } else if (method == "urn") {
        treatment                          <- numeric_n
        balls                              <- rep(w, Kp1)
        for (i in 1:n) {
          treatment[i]                     <-
            sample.int(Kp1, 1, prob = balls/sum(balls)) - 1L
          balls[treatment[i] + 1L]         <- balls[treatment[i] + 1L] + a
          balls[-(treatment[i] + 1L)]      <- balls[-(treatment[i] + 1L)] + b
        }
      } else if (method == "block urn") {
        treatment                           <- numeric_n
        balls_active                        <- rep(B/Kp1, Kp1)
        balls_inactive                      <- numeric(Kp1)
        for (i in 1:n) {
          treatment[i]                      <-
            sample.int(Kp1, 1, prob = balls_active/sum(balls_active)) - 1L
          balls_active[treatment[i] + 1L]   <- balls_active[treatment[i] + 1L] - 1L
          balls_inactive[treatment[i] + 1L] <- balls_active[treatment[i] + 1L] + 1L
          if (all(balls_inactive > 0)) {
            balls_active                    <- balls_active + 1L
            balls_inactive                  <- balls_inactive - 1L
          }
        }
      }
      treatment
    }
  
  measurements              <- function(rep, measures, data, J, K, seq_J, seq_K,
                                        seq_n, numeric_n, numeric_Kp1, seq_Kw0,
                                        covariate_names) {
    assignments                             <- list()
    for (k in seq_Kw0) {
      assignments[[k + 1]]                  <- (data$treatment == k)
    }
    treatment0                              <- sum(assignments[[1]])
    for (k in seq_K) {
      measures$measure1[rep, k]             <-
        sum(assignments[[k + 1]]) - treatment0
    }
    measures$measure1[rep, ]                <- abs(measures$measure1[rep, ])
    for (k in seq_K) {
      measures$measure7[rep, k]             <-
        (sd)^2/(sum(assignments[[k + 1]]) + treatment0)-(sum(assignments[[k + 1]]) - treatment0)^2/(sum(assignments[[k + 1]]) + treatment0) 
    }
    measures$measure7[rep, ]                <- abs(measures$measure7[rep, ])
    for (j in seq_J) {
      control_j                             <-
        data[[covariate_names[j]]][assignments[[1]]]
      if (length(control_j) > 0) {
        for (k in seq_K) {
          arm_k                               <- data[[covariate_names[j]]][assignments[[k + 1]]]
          treatment_k                         <- sum(assignments[[k + 1]])
          
          prop_control                        <- sum(control_j==1)/treatment0
          prop_arm                            <- sum(arm_k==1)/treatment_k
          
          if (length(arm_k) > 0) {
            measures$measure2[[j]][rep, k]    <- abs(prop_control - prop_arm)
          }
        }
      }
    }
    guesses                                 <- numeric_n
    running_totals                          <- numeric_Kp1
    for (i in seq_n) {
      which_min                             <-
        which(running_totals == min(running_totals))
      if (length(which_min) == 1) {
        guesses[i]                          <- which_min - 1L
      } else {
        guesses[i]                          <- sample(which_min, 1) - 1L
      }
      running_totals[data$treatment[i] + 1] <-
        running_totals[data$treatment[i] + 1] + 1L
    }
    measures$measure3[rep]                  <- sum(guesses == data$treatment)/n
    model                                   <-
      stats::lm(stats::formula(paste0("outcome ~ as.factor(treatment) + ",
                                      paste0("covariate", seq_J,
                                             collapse = " + "))), data)
    measures$measure4[rep, rownames(summary(model)$coefficients)] <-
      (summary(model)$coefficients[, 4] <= 0.05)
    measures
  }
  
  summary_measurements      <- function(measures, J, K, seq_J, seq_K) {
    performance                         <- list()
    performance$mean                    <- performance$sd <- list()
    performance$mean$measure1           <- colMeans(measures$measure1)
    performance$sd$measure1             <- numeric(K)
    names(performance$mean$measure1)    <- names(performance$sd$measure1) <-
      paste0("measure1_arm", seq_K)
    performance$mean$measure7           <- colMeans(measures$measure7)
    performance$sd$measure7             <- numeric(K)
    names(performance$mean$measure7)    <- names(performance$sd$measure7) <-
      paste0("measure7_arm", seq_K)
    for (k in seq_K) {
      performance$sd$measure1[k]        <- sd(measures$measure1[, k])
      performance$sd$measure7[k]        <- sd(measures$measure7[, k])
    }
    named_matrix                        <- matrix(0, J, K)
    rownames(named_matrix)              <- paste0("covariate", seq_J)
    colnames(named_matrix)              <- paste0("arm", seq_K)
    performance$max$measure2            <- performance$min$measure2        <-
      performance$mean$measure2 <- numeric(J)
    for (j in seq_J) {
      performance$min$measure2[j]       <- mean(apply(measures$measure2[[j]], 1, min))
      performance$max$measure2[j]       <- mean(apply(measures$measure2[[j]], 1, max))
      performance$mean$measure2[j]      <- mean(rowMeans(measures$measure2[[j]]))
    }
    performance$mean$measure3           <- mean(measures$measure3)
    performance$sd$measure3             <- sd(measures$measure3)
    names(performance$mean$measure3)    <- names(performance$sd$measure3) <-
      "measure3"
    performance$mean$measure4           <- colMeans(measures$measure4,
                                                    na.rm = TRUE)
    performance$sd$measure4             <- numeric(K + J + 1)
    names(performance$mean$measure4)    <- names(performance$sd$measure4) <-
      c("measure4_intercept", paste0("measure4_arm", seq_K),
        paste0("measure4_covariate", seq_J))
    for (k in 1:(K + J + 1)) {
      performance$sd$measure4[k]        <- sd(measures$measure4[, k],
                                              na.rm = TRUE)
    }
    performance$mean$measure5           <- numeric(nrow(measures$measure4))
    performance$mean$measure6           <- numeric(nrow(measures$measure4))
    for (i in 1:nrow(measures$measure4)) {
      performance$mean$measure5[i]      <- any(measures$measure4[i, 2:(K + 1)] == TRUE, na.rm = TRUE)
      performance$mean$measure6[i]      <- all(measures$measure4[i, 2:(K + 1)] == TRUE, na.rm = TRUE)
    }
    performance$mean$measure5           <- mean(performance$mean$measure5)
    performance$mean$measure6           <- mean(performance$mean$measure6)
    performance$sd$measure5             <- sd(performance$mean$measure5)
    performance$sd$measure6             <- sd(performance$mean$measure6)
    performance
  }
  
  K                         <- length(eff_tr)
  J                         <- length(prob_cov)
  J_rand                    <- num_cov_in_rand
  seq_J                     <- 1:J
  seq_J_rand                <- 1:J_rand
  seq_K                     <- 1:K
  seq_n                     <- 1:n
  base_data                 <-
    tibble::as_tibble(cbind(matrix(0L, n, J), 1L, 0L), .name_repair = "minimal")
  colnames(base_data)       <- c(paste0("covariate", seq_J), "strata",
                                 "treatment")
  measures                  <- list()
  measures$measure1         <- matrix(NA, replicates, K)
  measures$measure2         <- rep(list(measures$measure1), J)
  measures$measure3         <- numeric(replicates)
  measures$measure4         <- matrix(NA, replicates, K + J + 1)
  measures$measure7         <- matrix(NA, replicates, K)
  colnames(measures$measure4) <- c("(Intercept)",
                                   paste0("as.factor(treatment)", 1:K),
                                   paste0("covariate", 1:J))
  all_cov_combs             <- matrix(0, 2^J_rand, J_rand)
  for (j in seq_J_rand) {
    all_cov_combs[, j]      <- rep(0:1, each = 2^(j - 1),
                                   times = 2^(J_rand - j))
  }
  numeric_n                 <- numeric(n)
  Kp1                       <- K + 1
  numeric_Kp1               <- numeric(Kp1)
  seq_Kw0                   <- 0:K
  covariate_names           <- paste0("covariate", seq_J)
  if (method == "minimisation") {
    seq_nm1                 <- 2:n
    seq_js                  <- assignment_probs <- list()
    for (j in seq_nm1) {
      seq_js[[j]]           <- 1:(j - 1)
    }
    cov_wt                  <- rep(1/J_rand, J_rand)
    diag_ratio              <- diag(1/rep(1, Kp1))
    num_matching_cov        <- matrix(0, J_rand, Kp1)
    imbalance               <- numeric(Kp1)
    seq_Kp1                 <- seq_Kw0 + 1
    for (k in seq_K) {
      assignment_probs[[k]] <- c(rep(p/k, k), rep((1 - p)/(Kp1 - k), Kp1 - k))
    }
  } else {
    seq_js                  <- cov_wt    <- diag_ratio <- num_matching_cov <-
      imbalance <- seq_nm1    <- seq_Kp1          <-
      assignment_probs        <- NULL
  }
  for (rep in 1:replicates) {
    data                    <- base_data
    for (j in seq_J) {
      data[, j]             <- stats::rbinom(n, 1, prob_cov[j])
    }
    for (j in seq_J_rand) {
      data$strata           <- data$strata +
        data[[paste0("covariate", j)]]*as.integer(2^(j - 1))
    }
    data$treatment          <-
      treatment_assignment(data, K, n, method, B, p, J, all_cov_combs,
                           numeric_n, seq_J, Kp1, seq_Kw0, seq_js, cov_wt,
                           diag_ratio, num_matching_cov, imbalance, seq_nm1,
                           seq_Kp1, assignment_probs, J_rand, seq_J_rand,
                           burn_in, w, a, b)
    data$outcome            <-
      Rfast::rowsums(as.matrix(data[, seq_J])*matrix(eff_cov, n, J,
                                                     byrow = TRUE)) +
      c(0, eff_tr)[data$treatment + 1]
    data$outcome            <- stats::rnorm(n, mean = data$outcome, sd = sd)
    measures                <- measurements(rep, measures, data, J, K, seq_J,
                                            seq_K, seq_n, numeric_n,
                                            numeric_Kp1, seq_Kw0,
                                            covariate_names)
    message("done replicate ", rep)
  }
  summary_measurements(measures, J, K, seq_J, seq_K)
}


st_time1     <- Sys.time()
set.seed(123)
simple       <- evaluate_randomisation_routine(eff_tr = rep(1.2, 5), n = 85,
                                               method = "simple",
                                               prob_cov = rep(0.25, 4),
                                               eff_cov = 1.2*c(1, 0.5, 0.1, 0),
                                               num_cov_in_rand = 4,
                                               replicates = 10000)

end_time1    <- Sys.time()
end_time1 - st_time1
st_time2     <- Sys.time()
block        <- evaluate_randomisation_routine(eff_tr = rep(0.006, 6), n = 175,
                                               method = "block",
                                               prob_cov = rep(0.25, 4),
                                               eff_cov = 0.006*c(1, 0.5, 0.1, 0),
                                               num_cov_in_rand = 2,
                                               B = 7,
                                               replicates = 10000)
end_time2    <- Sys.time()
end_time2 - st_time2
st_time3     <- Sys.time()
stratified   <- evaluate_randomisation_routine(eff_tr = rep(0.006, 6), n = 350,
                                               method = "stratified",
                                               prob_cov = rep(0.25, 4),
                                               eff_cov = 0.006*c(1, 0.5, 0.1, 0),
                                               num_cov_in_rand = 2,
                                               B = 21,
                                               replicates = 10000)
end_time3    <- Sys.time()
end_time3 - st_time3
st_time4     <- Sys.time()
minimisation <- evaluate_randomisation_routine(eff_tr = rep(1.2, 5), n = 170,
                                               method = "minimisation",
                                               prob_cov = rep(0.25, 4),
                                               eff_cov = 1.2*c(1, 0.5, 0.1, 0),
                                               num_cov_in_rand = 2,
                                               p = 0.9, burn_in = ceiling(0.1*170),
                                               replicates = 10000)
end_time4    <- Sys.time()
end_time4 - st_time4
st_time5     <- Sys.time()
urn          <- evaluate_randomisation_routine(eff_tr = rep(1.2, 5), n = 43,
                                               method = "urn",
                                               eff_cov = 1.2*c(1, 0.5, 0.1, 0),
                                               prob_cov = rep(0.25, 4),
                                               num_cov_in_rand = 4,
                                               w = 1, a = 1, b = 2,
                                               replicates = 10000)
end_time5    <- Sys.time()
end_time5 - st_time5
st_time6     <- Sys.time()
block_urn    <- evaluate_randomisation_routine(eff_tr = rep(1.2, 5), n = 170,
                                               method = "block urn",
                                               prob_cov = rep(0.25, 4),
                                               eff_cov = 1.2*c(1, 0.5, 0.1, 0),
                                               num_cov_in_rand = 2,
                                               w=1, a=1, b=2,
                                               B = 18,
                                               replicates = 10000)
