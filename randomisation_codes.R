################################################################################
# Filname:       figures_and_tables.R                                          #
# Author(s):     Ruqayya Azher (r.a.o.azher2@newcastle.ac.uk),                 #
#                Michael Grayling (michael.grayling@newcastle.ac.uk)           #
# Last modified: 2023/06/01                                                    #
# Description:   This simulation codes evaluates the performance of several    #
#                randomisation methods within the context of a multi-arm       #
#                trial. It can be used to reproduce the results from "Azher RA #
#                et al (2023) A comparison of randomization methods for        #
#                multi-arm clinical trials"                                    #
################################################################################

##### INSTALL / LOAD REQUIRED PACKAGES #########################################

#install.packages("patchwork")
#install.packages("Rfast")
#install.packages("snowfall")
#install.packages("tidyverse")
library(patchwork)
library(Rfast)
library(snowfall)
library(tidyverse)

##### HELPER FUNCTIONS #########################################################
# These functions are used by evaluate_randomisation_routine()

# Assigns treatments
treatment_assignment <- function(data, K, n, method, B, p, J, J_rand, burn_in,
                                 w, a, b) {
    treatment                               <- numeric(n)
    if (method == "SR") {
      treatment                             <- sample.int(K + 1, n,
                                                          replace = TRUE) - 1L
    } else if (method == "PBR") {
      treatment                             <-
        as.vector(sapply(1:ceiling(n/B),
                         function(i) { sample(rep(0:K, B/(K + 1)), B) }))[1:n]
    } else if (method == "SBR") {
      blocks                                <-
        matrix(as.vector(sapply(1:(2^J_rand*ceiling(n/B)),
                                function(i) { sample(rep(0:K, B/(K + 1)),
                                                     B) })),
               B*ceiling(n/B), 2^J_rand)
      for (i in 1:(2^J_rand)) {
        which_strata_i                      <- which(data$strata == i)
        treatment[which_strata_i]           <-
          blocks[1:length(which_strata_i), i]
      }
    } else if (method == "Mini") {
      covariate_matrix                      <- as.matrix(data[, 1:J_rand])
      num_matching_covariates               <- matrix(0, J_rand, K + 1)
      treatment[1:burn_in]                  <- sample.int(K + 1, burn_in,
                                                          replace = TRUE) - 1L
      for (j in (burn_in + 1):n) {
        matching_covariates                 <-
          (covariate_matrix[1:(j - 1), , drop = FALSE] ==
             matrix(covariate_matrix[j, ], j - 1, J_rand, byrow = TRUE))
        for (k in 1:(K + 1)) {
          num_matching_covariates[, k]      <-
            Rfast::colsums(matching_covariates[treatment[1:(j - 1)] == k - 1, ,
                                        drop = FALSE])
        }
        imbalance                           <- numeric(K + 1)
        for (k in 1:(K + 1)) {
          temp                              <- num_matching_covariates
          temp[, k]                         <- temp[, k] + 1
          range_level                       <-
            apply(temp%*%diag(1/rep(1, K + 1)), 1, range)
          imbalance[k]                      <-
            sum(rep(1/J_rand, J_rand)*(range_level[2, ] - range_level[1, ]))
        }
        minimizers                          <-
          (0:K)[imbalance == min(imbalance)]
        len_minimizers                      <- length(minimizers)
        if (len_minimizers == K + 1) {
          treatment[j]                      <- sample.int(K + 1, 1) - 1L
        } else {
          treatment[j]                      <-
            sample(c(minimizers, (0:K)[-(minimizers + 1)]), 1,
                   prob = c(rep(p/len_minimizers, len_minimizers),
                            rep((1 - p)/(K + 1 - len_minimizers),
                                K + 1 - len_minimizers)))
        }
      }
    } else if (method == "UD") {
      balls                                 <- rep(w, K + 1)
      for (i in 1:n) {
        treatment[i]                        <-
          sample.int(K + 1, 1, prob = balls/sum(balls)) - 1L
        balls[treatment[i] + 1L]            <- balls[treatment[i] + 1L] + a
        balls[-(treatment[i] + 1L)]         <- balls[-(treatment[i] + 1L)] + b
      }
    } else if (method == "BUD") {
      balls_active                          <- rep(B/(K + 1), K + 1)
      balls_inactive                        <- numeric(K + 1)
      for (i in 1:n) {
        treatment[i]                        <-
          sample.int(K + 1, 1, prob = balls_active/sum(balls_active)) - 1L
        balls_active[treatment[i] + 1L]     <-
          balls_active[treatment[i] + 1L] - 1L
        balls_inactive[treatment[i] + 1L]   <-
          balls_inactive[treatment[i] + 1L] + 1L
        if (all(balls_inactive > 0)) {
          balls_active                      <- balls_active + 1L
          balls_inactive                    <- balls_inactive - 1L
        }
      }
    } else if (method == "SBUD") {
      for (i in 1:(2^J_rand)) {
        which_strata_i                      <- which(data$strata == i)
        balls_active                        <- rep(B/(K + 1), (K + 1))
        balls_inactive                      <- numeric(K + 1)
        for (j in which_strata_i) {
          treatment[j]                      <-
            sample.int(K + 1, 1, prob = balls_active/sum(balls_active)) - 1L
          balls_active[treatment[j] + 1L]   <-
            balls_active[treatment[j] + 1L] - 1L
          balls_inactive[treatment[j] + 1L] <-
            balls_inactive[treatment[j] + 1L] + 1L
          if (all(balls_inactive > 0)) {
            balls_active                    <- balls_active + 1L
            balls_inactive                  <- balls_inactive - 1L
          }
        }
      }
    }
  treatment
}

# Computes performance measures for a single simulation replicate
get_measures         <- function(rep, measures, data, n, J, K) {
  ns                                      <- numeric(K + 1)
  for (k in 0:K) {
    ns[k + 1]                             <- sum(data$treatment == k)
  }
  measures$measure1[rep, ]                <- abs(ns[-1] - ns[1])
  for (j in 1:J) {
    covariates_j0                         <-
      data[[paste0("covariate", j)]][data$treatment == 0]
    if (length(covariates_j0) > 0) {
      for (k in 1:K) {
        covariates_jk                     <-
          data[[paste0("covariate", j)]][data$treatment == k]
        if (length(covariates_jk) > 0) {
          measures$measure2[[j]][rep, k]  <-
            abs(sum(covariates_jk == 1)/ns[k + 1] -
                  sum(covariates_j0 == 1)/ns[1])
        }
      }
    }
  }
  guesses                                 <- numeric(n)
  running_totals                          <- numeric(K + 1)
  for (i in 1:n) {
    guesses[i]                            <-
      sample(which(running_totals == min(running_totals)), 1) - 1L
    running_totals[data$treatment[i] + 1] <-
      running_totals[data$treatment[i] + 1] + 1L
  }
  measures$measure3[rep]                  <- sum(guesses == data$treatment)/n
  if (all(ns > 0)) {
    model                                 <-
      stats::lm(stats::formula(paste0("outcome ~ as.factor(treatment) + ",
                                      paste0("covariate", 1:J,
                                             collapse = " + "))), data)
    measures$measure4[rep, rownames(summary(model)$coefficients)] <-
      (summary(model)$coefficients[, 4] <= 0.05)
    measures$measure7[rep, ]                                      <-
      summary(model)$coefficients[paste0("as.factor(treatment)", 1:K), 2]^2
  }
  measures
}

# Summarises performance measures across simulations replicates
summarise_measures   <- function(measures, n, J, K, sd, replicates) {
  performance               <- list()
  performance$measure1      <- mean(apply(measures$measure1, 1, max,
                                          na.rm = TRUE))
  performance$measure2      <- rep(NA, replicates)
  for (i in 1:replicates) {
    max_imbalances          <- rep(NA, J)
    for (j in 1:J) {
      max_imbalances[j]     <- max(measures$measure2[[j]][i, ], na.rm = TRUE)
    }
    performance$measure2[i] <- max(max_imbalances, na.rm = TRUE)
  }
  performance$measure2[is.infinite(performance$measure2)] <- NA
  performance$measure2      <- mean(performance$measure2, na.rm = TRUE)
  performance$measure3      <- 100*mean(measures$measure3, na.rm = TRUE)
  performance$measure4      <- mean(measures$measure4[, 2], na.rm = TRUE)
  performance$measure5      <-
    mean(Rfast::rowsums(measures$measure4[, 2:(K + 1)], na.rm = TRUE) > 0)
  performance$measure6      <-
    mean(Rfast::rowsums(measures$measure4[, 2:(K + 1)], na.rm = TRUE) == K)
  performance$measure7      <-
    100*apply(measures$measure7, 1, max,
               na.rm = TRUE)/(2*sd^2/(n/(K + 1))) - 100
  performance$measure7[is.infinite(performance$measure7)] <- NA
  performance$measure7      <- mean(performance$measure7, na.rm = TRUE)
  unlist(performance)
}

evaluate_randomisation_routine <- function(
    
  target_effect = 0.006,
  
  eff_tr        = rep(target_effect, 6), # The treatment effect for each experimental
  
  n             = 350, # Total sample size
  
  method        = "SR", # Method for randomising - can be "simple", # "block", "stratified", "minimisation"
  
  prob_cov      = rep(0.25, 4), # Probability of a patient having each # covariate - the length of this determines the number of covariates
  
  eff_cov       = target_effect*c(1, 0.5, 0.1, 0), # The effect on the outcome data of each  covariate
  
  J_rand        = 2, # The number of covariates used in the randomisation routine (assumed to be the first num_cov_in_rand covariates) - only used if method is "stratified" or "minimisation"
  
  B             = 6, # Length of the blocks - only used if method is  "block" or "stratified" or "block urn"
  
  p             = 0.7, # Only used when method is "minimisation"
  
  burn_in       = ceiling(0.1*n), # Only used when method is "minimisation"
  
  w             = 1, # Only used method method is "urn" - signifies number of balls initially in urn of each colour
  
  a             = 1, # Only used method method is "urn"
  
  b             = 2, # Only used method method is "urn"
  
  sd            = 0.011,
  
  replicates    = 10000 # Number of replicate simulations
  ) {
  K                           <- length(eff_tr)
  J                           <- length(prob_cov)
  base_data                   <-
    tibble::as_tibble(cbind(matrix(0L, n, J), 1L, 0L), .name_repair = "minimal")
  colnames(base_data)         <- c(paste0("covariate", 1:J), "strata",
                                   "treatment")
  measures                    <- list()
  measures$measure1           <- matrix(NA, replicates, K)
  measures$measure2           <- rep(list(measures$measure1), J)
  measures$measure3           <- numeric(replicates)
  measures$measure4           <- matrix(NA, replicates, K + J + 1)
  measures$measure7           <- matrix(NA, replicates, K)
  colnames(measures$measure4) <- c("(Intercept)",
                                   paste0("as.factor(treatment)", 1:K),
                                   paste0("covariate", 1:J))
  for (rep in 1:replicates) {
    data                      <- base_data
    for (j in 1:J) {
      data[, j]               <- stats::rbinom(n, 1, prob_cov[j])
    }
    for (j in 1:J_rand) {
      data$strata             <- data$strata +
        data[[paste0("covariate", j)]]*as.integer(2^(j - 1))
    }
    data$treatment            <-
      treatment_assignment(data, K, n, method, B, p, J, J_rand, burn_in, w, a,
                           b)
    data$outcome              <-
      Rfast::rowsums(as.matrix(data[, 1:J])*matrix(eff_cov, n, J,
                                                   byrow = TRUE)) +
      c(0, eff_tr)[data$treatment + 1]
    data$outcome              <- stats::rnorm(n, data$outcome, sd)
    measures                  <- get_measures(rep, measures, data, n, J, K)
    if (rep%%100 == 0) message("..done replicate ", rep, "..")
  }
  summarise_measures(measures, n, J, K, sd, replicates)
}

wrapper <- function(i) {
  set.seed(i)
  evaluate_randomisation_routine(
    target_effect = scenarios$target_effect[i],
    eff_tr        = rep(scenarios$target_effect[i], scenarios$K[i]),
    n             = scenarios$n[i],
    method        = scenarios$method[i],
    J_rand        = scenarios$J_rand[i],
    B             = scenarios$B[i],
    p             = scenarios$p[i],
    sd            = scenarios$sd[i],
    replicates    = scenarios$replicates[i])
}

##### EXAMPLES #################################################################

replicates_global <- 10000 

scenarios1        <-
  expand.grid(n             = c(175, 350, 700),
              method        = c("SR", "PBR", "SBR", "Mini", "UD", "BUD",
                                "SBUD"),
              J_rand        = c(2, 4),
              B             = c(7, 21),
              p             = c(0.7, 0.9),
              target_effect = 0.006,
              sd            = 0.011,
              K             = 6,
              replicates    = replicates_global)
scenarios1        <-
  scenarios1[-which(scenarios1$method %in% c("SR", "Mini", "UD") &
                     scenarios1$B > 7), ]
scenarios1        <-
  scenarios1[-which(scenarios1$method %in% c("BUD", "SBUD") &
                      scenarios1$B != 21), ]
scenarios1        <- scenarios1[-which(scenarios1$method != "Mini" &
                                         scenarios1$p == 0.9), ]
scenarios1        <- dplyr::arrange(scenarios1, n, J_rand, B, p, method)
scenarios         <- scenarios1
sfInit(parallel = TRUE, cpus = 18)
sfLibrary(Rfast)
sfLibrary(tibble)
sfExport("treatment_assignment", "get_measures", "summarise_measures",
         "evaluate_randomisation_routine", "scenarios")
results1           <- sfLapply(1:nrow(scenarios1), wrapper)
sfStop()
results_mat1       <- matrix(unlist(results1), nrow(scenarios1), 7,
                             byrow = TRUE)
methods1           <- character(nrow(scenarios1))
for (i in 1:nrow(scenarios1)) {
  if (scenarios1$method[i] == "SR") {
    methods1[i]    <- as.character(scenarios1$method[i])
  } else if (scenarios1$method[i] %in% c("PBR", "SBR")) {
    if (scenarios1$B[i] == 7) {
      methods1[i]  <- paste0(scenarios1$method[i], "(K)")
    } else {
      methods1[i]  <- paste0(scenarios1$method[i], "(3K)")
    }
  } else if (scenarios1$method[i] == "Mini") {
    methods1[i]    <- paste0("Mini(", scenarios1$p[i], ")")
  } else if (scenarios1$method[i] == "UD") {
    methods1[i]    <- "UD(1,1,2)"
  } else if (scenarios1$method[i] %in% c("BUD", "SBUD")) {
    methods1[i]    <- paste0(scenarios1$method[i], "(3)")
  }
}
df1                <- tibble::tibble(
  Setting                                           = "Setting~1",
  `Sample size`                                     = scenarios1$n,
  Method                                            =
    factor(methods1, c("BUD(3)", "Mini(0.7)", "Mini(0.9)", "PBR(K)", "PBR(3K)",
                       "SBR(K)", "SBR(3K)", "SBUD(3)", "SR", "UD(1,1,2)")),
  J_rand                                            =
    paste("italic(J)[rand] ==", scenarios1$J_rand),
  `Maximum group size imbalance`                    = results_mat1[, 1],
  `Maximum covariate imbalance`                     = results_mat1[, 2],
  `Predictability (%)`                              = results_mat1[, 3],
  `Marginal~power`                                  = results_mat1[, 4],
  `Disjunctive~power`                               = results_mat1[, 5],
  `Conjunctive~power`                               = results_mat1[, 6],
  `Maximum treatment effect variance inflation (%)` = results_mat1[, 7])
figure1_1         <-
  ggplot(df1, aes(`Sample size`, `Maximum group size imbalance`,
                  colour = Method)) +
  facet_grid(Setting ~ J_rand, labeller = label_parsed) +
  geom_line() +
  geom_point() +
  theme_bw() +
  theme(legend.position = "bottom",
        legend.title    = element_blank(),
        aspect.ratio    = 1) +
  xlab(expression(paste("Sample size (", italic(N), ")", sep = ""))) +
  guides(colour = guide_legend(nrow = 2, byrow = TRUE))
figure2_1         <-
  ggplot(df1, aes(`Sample size`, `Maximum covariate imbalance`,
                  colour = Method)) +
  facet_grid(Setting ~ J_rand, labeller = label_parsed) +
  geom_line() +
  geom_point() +
  theme_bw() +
  theme(legend.position = "bottom",
        legend.title    = element_blank(),
        aspect.ratio    = 1) +
  xlab(expression(paste("Sample size (", italic(N), ")", sep = ""))) +
  guides(colour = guide_legend(nrow = 2, byrow = TRUE))
figure3_1         <-
  ggplot(df1, aes(`Sample size`,
                  `Maximum treatment effect variance inflation (%)`,
                 colour = Method)) +
  facet_grid(Setting ~ J_rand, labeller = label_parsed) +
  geom_line() +
  geom_point() +
  theme_bw() +
  theme(legend.position = "bottom",
        legend.title    = element_blank(),
        aspect.ratio    = 1,
        axis.title.y    = element_text(size = 8)) +
  xlab(expression(paste("Sample size (", italic(N), ")", sep = ""))) +
  guides(colour = guide_legend(nrow = 2, byrow = TRUE))
figure4_1         <-
  ggplot(df1, aes(`Sample size`, `Predictability (%)`, colour = Method)) +
  facet_grid(Setting ~ J_rand, labeller = label_parsed) +
  geom_line() +
  geom_point() +
  theme_bw() +
  theme(legend.position = "bottom",
        legend.title    = element_blank(),
        aspect.ratio    = 1) +
  xlab(expression(paste("Sample size (", italic(N), ")", sep = ""))) +
  guides(colour = guide_legend(nrow = 2, byrow = TRUE))
df1_gathered      <-
  tidyr::pivot_longer(df1, `Marginal~power`:`Conjunctive~power`)
figure5_1         <-
  ggplot(dplyr::filter(df1_gathered, J_rand == "italic(J)[rand] == 4"),
         aes(as.factor(`Sample size`), 100*value, fill = Method)) +
  geom_bar(position = "dodge", stat = "identity") +
  facet_grid(Setting ~ name, labeller = label_parsed) +
  xlab(expression(paste("Sample size (", italic(N), ")", sep = ""))) +
  ylab("Power (%)") +
  theme_bw() +
  theme(legend.position = "bottom",
        legend.title    = element_blank(),
        aspect.ratio    = 1) +
  guides(fill = guide_legend(nrow = 2, byrow = TRUE))

scenarios2        <-
  expand.grid(n             = c(43, 85, 170),
              method        = c("SR", "PBR", "SBR", "Mini", "UD", "BUD",
                                "SBUD"),
              J_rand        = c(2, 4),
              B             = c(6, 18),
              p             = c(0.7, 0.9),
              target_effect = 1.2,
              sd            = 1,
              K             = 5,
              replicates    = replicates_global)
scenarios2        <-
  scenarios2[-which(scenarios2$method %in% c("SR", "Mini", "UD") &
                      scenarios2$B > 6), ]
scenarios2        <- scenarios2[-which(scenarios2$method %in% c("BUD", "SBUD") &
                                         scenarios2$B != 18), ]
scenarios2        <- scenarios2[-which(scenarios2$method != "Mini" &
                                         scenarios2$p == 0.9), ]
scenarios2        <- dplyr::arrange(scenarios2, n, J_rand, B, p, method)
scenarios         <- scenarios2
sfInit(parallel = TRUE, cpus = 18)
sfLibrary(Rfast)
sfLibrary(tibble)
sfExport("treatment_assignment", "get_measures", "summarise_measures",
         "evaluate_randomisation_routine", "scenarios")
results2          <- sfLapply(1:nrow(scenarios2), wrapper)
sfStop()
results_mat2      <- matrix(unlist(results2), nrow(scenarios2), 7, byrow = TRUE)
methods2          <- character(nrow(scenarios2))
for (i in 1:nrow(scenarios2)) {
  if (scenarios2$method[i] == "SR") {
    methods2[i]   <- as.character(scenarios2$method[i])
  } else if (scenarios2$method[i] %in% c("PBR", "SBR")) {
    if (scenarios2$B[i] == 6) {
      methods2[i] <- paste0(scenarios2$method[i], "(K)")
    } else {
      methods2[i] <- paste0(scenarios2$method[i], "(3K)")
    }
  } else if (scenarios2$method[i] == "Mini") {
    methods2[i]    <- paste0("Mini(", scenarios2$p[i], ")")
  } else if (scenarios2$method[i] == "UD") {
    methods2[i]    <- "UD(1,1,2)"
  } else if (scenarios2$method[i] %in% c("BUD", "SBUD")) {
    methods2[i]    <- paste0(scenarios2$method[i], "(3)")
  }
}
df2                <- tibble::tibble(
  Setting                                           = "Setting~2",
  `Sample size`                                     = scenarios2$n,
  Method                                            =
    factor(methods2, c("BUD(3)", "Mini(0.7)", "Mini(0.9)", "PBR(K)", "PBR(3K)",
                       "SBR(K)", "SBR(3K)", "SBUD(3)", "SR", "UD(1,1,2)")),
  J_rand                                            =
    paste("italic(J)[rand] ==", scenarios$J_rand),
  `Maximum group size imbalance`                    = results_mat2[, 1],
  `Maximum covariate imbalance`                     = results_mat2[, 2],
  `Predictability (%)`                              = results_mat2[, 3],
  `Marginal~power`                                  = results_mat2[, 4],
  `Disjunctive~power`                               = results_mat2[, 5],
  `Conjunctive~power`                               = results_mat2[, 6],
  `Maximum treatment effect variance inflation (%)` = results_mat2[, 7])
figure1_2         <-
  ggplot(df2, aes(`Sample size`, `Maximum group size imbalance`,
                  colour = Method)) +
  facet_grid(Setting ~ J_rand, labeller = label_parsed) +
  geom_line() +
  geom_point() +
  theme_bw() +
  theme(legend.position = "bottom",
        legend.title    = element_blank(),
        aspect.ratio    = 1) +
  xlab(expression(paste("Sample size (", italic(N), ")", sep = ""))) +
  guides(colour = guide_legend(nrow = 2, byrow = TRUE))
figure2_2         <-
  ggplot(df2, aes(`Sample size`, `Maximum covariate imbalance`,
                  colour = Method)) +
  facet_grid(Setting ~ J_rand, labeller = label_parsed) +
  geom_line() +
  geom_point() +
  theme_bw() +
  theme(legend.position = "bottom",
        legend.title    = element_blank(),
        aspect.ratio    = 1) +
  xlab(expression(paste("Sample size (", italic(N), ")", sep = ""))) +
  guides(colour = guide_legend(nrow = 2, byrow = TRUE))
figure3_2         <-
  ggplot(df2, aes(`Sample size`,
                  `Maximum treatment effect variance inflation (%)`,
                  colour = Method)) +
  facet_grid(Setting ~ J_rand, labeller = label_parsed) +
  geom_line() +
  geom_point() +
  theme_bw() +
  theme(legend.position = "bottom",
        legend.title    = element_blank(),
        aspect.ratio    = 1,
        axis.title.y    = element_text(size = 8)) +
  xlab(expression(paste("Sample size (", italic(N), ")", sep = ""))) +
  guides(colour = guide_legend(nrow = 2, byrow = TRUE))
figure4_2         <-
  ggplot(df2, aes(`Sample size`, `Predictability (%)`, colour = Method)) +
  facet_grid(Setting ~ J_rand, labeller = label_parsed) +
  geom_line() +
  geom_point() +
  theme_bw() +
  theme(legend.position = "bottom",
        legend.title    = element_blank(),
        aspect.ratio    = 1) +
  xlab(expression(paste("Sample size (", italic(N), ")", sep = ""))) +
  ylab("Predictability (%)") +
  guides(colour = guide_legend(nrow = 2, byrow = TRUE))
df2_gathered      <-
  tidyr::pivot_longer(df2, `Marginal~power`:`Conjunctive~power`)
figure5_2         <-
  ggplot(dplyr::filter(df2_gathered, J_rand == "italic(J)[rand] == 4"),
         aes(as.factor(`Sample size`), 100*value, fill = Method)) +
  geom_bar(position = "dodge", stat = "identity") +
  facet_grid(Setting ~ name, labeller = label_parsed) +
  xlab(expression(paste("Sample size (", italic(N), ")", sep = ""))) +
  ylab("Power (%)") +
  theme_bw() +
  theme(legend.position = "bottom",
        legend.title    = element_blank(),
        aspect.ratio    = 1) +
  guides(fill = guide_legend(nrow = 2, byrow = TRUE))

figure1           <- figure1_1 / figure1_2 +
  plot_layout(guides = "collect", widths = c(0.5,0.5)) &
  theme(legend.position = "bottom")
figure2           <- figure2_1 / figure2_2 +
  plot_layout(guides = "collect", widths = c(0.5,0.5)) &
  theme(legend.position = "bottom")
figure3           <- figure3_1 / figure3_2 +
  plot_layout(guides = "collect", widths = c(0.5,0.5)) &
  theme(legend.position = "bottom")
figure4           <- figure4_1 / figure4_2 +
  plot_layout(guides = "collect", widths = c(0.5,0.5)) &
  theme(legend.position = "bottom")
figure5           <- figure5_1 / figure5_2 +
  plot_layout(guides = "collect", widths = c(0.5,0.5)) &
  theme(legend.position = "bottom")
ggsave("figure1.pdf", figure1, device = "pdf", width = 6, height = 6,
       units = "in")
ggsave("figure2.pdf", figure2, device = "pdf", width = 6, height = 6,
       units = "in")
ggsave("figure3.pdf", figure3, device = "pdf", width = 6, height = 6,
       units = "in")
ggsave("figure4.pdf", figure4, device = "pdf", width = 6, height = 6,
       units = "in")
ggsave("figure5.pdf", figure5, device = "pdf", width = 6, height = 6,
       units = "in")

scenarios3        <-
  expand.grid(n             = seq(125, 1250, 125),
              method        = "SR",
              J_rand        = 4,
              B             = 7,
              p             = 0.7,
              target_effect = 0.006,
              sd            = 0.011,
              K             = 6,
              replicates    = replicates_global)
scenarios         <- scenarios3
sfInit(parallel = TRUE, cpus = 18)
sfLibrary(Rfast)
sfLibrary(tibble)
sfExport("treatment_assignment", "get_measures", "summarise_measures",
         "evaluate_randomisation_routine", "scenarios")
results3           <- sfLapply(1:nrow(scenarios3), wrapper)
sfStop()
results_mat3       <- matrix(unlist(results3), nrow(scenarios3), 7,
                             byrow = TRUE)
df3                <- tibble::tibble(
  Setting                                           = "Setting~1",
  `Sample size`                                     = scenarios3$n,
  `Marginal power`                                  = results_mat3[, 4],
  `Disjunctive power`                               = results_mat3[, 5],
  `Conjunctive power`                               = results_mat3[, 6])

scenarios4        <-
  expand.grid(n             = seq(24, 24*8, 24),
              method        = "SR",
              J_rand        = 4,
              B             = 7,
              p             = 0.7,
              target_effect = 1.2,
              sd            = 1,
              K             = 5,
              replicates    = replicates_global)
scenarios         <- scenarios4
sfInit(parallel = TRUE, cpus = 18)
sfLibrary(Rfast)
sfLibrary(tibble)
sfExport("treatment_assignment", "get_measures", "summarise_measures",
         "evaluate_randomisation_routine", "scenarios")
results4           <- sfLapply(1:nrow(scenarios4), wrapper)
sfStop()
results_mat4       <- matrix(unlist(results4), nrow(scenarios4), 7,
                             byrow = TRUE)
df4                <- tibble::tibble(
  Setting                                           = "Setting~2",
  `Sample size`                                     = scenarios4$n,
  `Marginal power`                                  = results_mat4[, 4],
  `Disjunctive power`                               = results_mat4[, 5],
  `Conjunctive power`                               = results_mat4[, 6])

df3_gathered      <-
  tidyr::pivot_longer(df3, `Marginal power`:`Conjunctive power`)
figure6_1         <-
  ggplot(df3_gathered,
         aes(`Sample size`, 100*value, colour = name)) +
  geom_point() +
  geom_line() +
  facet_grid(Setting ~ ., labeller = label_parsed) +
  xlab(expression(paste("Sample size (", italic(N), ")", sep = ""))) +
  ylab("Power (%)") +
  theme_bw() +
  theme(legend.position = "bottom",
        legend.title    = element_blank(),
        aspect.ratio    = 1)
df4_gathered      <-
  tidyr::pivot_longer(df4, `Marginal power`:`Conjunctive power`)
figure6_2         <-
  ggplot(df4_gathered,
         aes(`Sample size`, 100*value, colour = name)) +
  geom_point() +
  geom_line() +
  facet_grid(Setting ~ ., labeller = label_parsed) +
  xlab(expression(paste("Sample size (", italic(N), ")", sep = ""))) +
  ylab("Power (%)") +
  theme_bw() +
  theme(legend.position = "bottom",
        legend.title    = element_blank(),
        aspect.ratio    = 1)

figure6           <- figure6_1 / figure6_2 +
  plot_layout(guides = "collect", widths = c(0.5,0.5)) &
  theme(legend.position = "bottom")
ggsave("figure6.pdf", figure6, device = "pdf", width = 4.5, height = 6,
       units = "in")

