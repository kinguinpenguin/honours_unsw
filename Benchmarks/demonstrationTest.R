#  File tests/termTests_rank.R in package ergm.rank, part of the Statnet suite
#  of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free, open
#  source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2008-2025 Statnet Commons
################################################################################
start_time <- Sys.time()
library(ergm.rank)
data("newcomb")
fit_true <- list()
fit_constrained <- list()
for (week in 1:15) {
  cat(sprintf("Week %d\n", week))
  Y_full <- newcomb[[week]]
  rank_matrix <- as.matrix(Y_full, attrname = "descrank")

  n <- nrow(rank_matrix)

  # Create top-5 observed matrix
  newcomb2.top5 <- matrix(1, nrow = n, ncol = n)

  # For each ego, keep only the top 5 highest descrank values
  for (i in 1:nrow(rank_matrix)) {
    valid_indices <- which(!is.na(rank_matrix[i, ]))
    if (length(valid_indices) > 0) {
      ordered <- valid_indices[order(rank_matrix[i, valid_indices], decreasing = TRUE)]
      top5_indices <- head(ordered, 5)
      newcomb2.top5[i, top5_indices] <- rank_matrix[i, top5_indices]
    }
  }

  newcomb2.top5
  nw2 <- simulate(newcomb[[week]] # Start with newcomb time point 2.
                  ~ sum, coef = 0, # A term is required, but its coefficient is set to 0.
                  response="descrank", # Use edge attribute "descrank" as response.
                  reference=~CompleteOrder, # Complete ordering.
                  constraints = ~ adjacent + ranking(newcomb2.top5)
  )


  fit2.true <- ergm(newcomb[[week]] ~ rank.deference+rank.nonconformity("all")+
                      rank.nonconformity("localAND"),
              response="descrank",
              reference=~CompleteOrder,
              constraints = ~ adjacent, # Sample space: make adjacent swap proposals.
              control = snctrl()
  )


  fit2 <- ergm(nw2 ~ rank.deference+rank.nonconformity("all")+
                rank.nonconformity("localAND"),
              response="descrank",
              reference=~CompleteOrder,
              constraints = ~ adjacent, # Sample space: make adjacent swap proposals.
              obs.constraints = ~ ranking(newcomb2.top5), # For the constrained sampler, *also* constrain ranking.
              control = snctrl(init.method = "zeros")
  )

  fit_true[[week]] <- fit2.true
  fit_constrained[[week]] <- fit2
}

# Create an empty data frame to store comparisons
compare_df <- data.frame(
  week = 1:15,
  coef_true = NA,
  coef_constrained = NA,
  se_true = NA,
  se_constrained = NA,
  logLik_true = NA,
  logLik_constrained = NA
)

# Loop through all 15 models
for (i in 1:15) {
  # Check that both fits exist and are valid ergm objects
  if (inherits(fit_true[[i]], "ergm") && inherits(fit_constrained[[i]], "ergm")) {
    # Coefficients
    compare_df$coef_true[i] <- coef(fit_true[[i]])
    compare_df$coef_constrained[i] <- coef(fit_constrained[[i]])

    # Standard errors
    compare_df$se_true[i] <- sqrt(diag(vcov(fit_true[[i]])))
    compare_df$se_constrained[i] <- sqrt(diag(vcov(fit_constrained[[i]])))

    # Log-likelihoods
    compare_df$logLik_true[i] <- as.numeric(logLik(fit_true[[i]]))
    compare_df$logLik_constrained[i] <- as.numeric(logLik(fit_constrained[[i]]))
  } else {
    # If one failed, mark NA
    compare_df[i, 2:7] <- NA
  }
}

# Print to console
print(compare_df)

# Save to file
capture.output(print(compare_df), file = "tests/demonstration_results.txt")

end_time <- Sys.time()
total_time <- as.numeric(difftime(end_time, start_time, units = "secs"))
cat("\nTotal computational time:", round(total_time, 2), "seconds\n",
    file = "tests/demonstration_results.txt", append = TRUE)
