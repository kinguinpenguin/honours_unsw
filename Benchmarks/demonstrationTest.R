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

compare_list <- list()

for (i in 1:15) {
  if (inherits(fit_true[[i]], "ergm") && inherits(fit_constrained[[i]], "ergm")) {
    # Extract coefficients and SEs
    coefs_true <- coef(fit_true[[i]])
    coefs_constrained <- coef(fit_constrained[[i]])

    se_true <- sqrt(diag(vcov(fit_true[[i]])))
    se_constrained <- sqrt(diag(vcov(fit_constrained[[i]])))

    # Log-likelihoods
    logLik_true <- as.numeric(logLik(fit_true[[i]]))
    logLik_constrained <- as.numeric(logLik(fit_constrained[[i]]))

    # Combine into a data frame (one row per parameter × model type)
    df_week <- data.frame(
      week = i,
      term = rep(names(coefs_true), 2),
      model_type = rep(c("Complete", "Top 5"), each = length(coefs_true)),
      estimate = c(coefs_true, coefs_constrained),
      se = c(se_true, se_constrained),
      logLik = rep(c(logLik_true, logLik_constrained), each = length(coefs_true))
    )

    compare_list[[i]] <- df_week
  }
}

compare_df <- do.call(rbind, compare_list)

# Print to console
print(compare_df)

# Save to file
capture.output(print(compare_df), file = "tests/demonstration_results.txt")

end_time <- Sys.time()
total_time <- as.numeric(difftime(end_time, start_time, units = "secs"))
cat("\nTotal computational time:", round(total_time, 2), "seconds\n",
    file = "tests/demonstration_results.txt", append = TRUE)

install.packages("ggplot2")
library(ggplot2)

# Compute confidence interval bounds
compare_df$lower <- compare_df$estimate - compare_df$se
compare_df$upper <- compare_df$estimate + compare_df$se
compare_df$week <- factor(compare_df$week)

# Get unique features (terms)
terms <- unique(compare_df$term)

library(ggplot2)

p <- ggplot(compare_df, aes(x = week, y = estimate,
                            color = model_type, group = model_type)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = model_type),
              alpha = 0.2, color = NA) +
  facet_wrap(~ term, ncol = 1, scales = "free_y") + 
  theme_minimal(base_size = 14) +
  labs(
    title = "Parameter estimates over time",
    x = "Week",
    y = "Estimate",
    color = "Model type",
    fill = "Model type"
  )

print(p)

# Save stacked version
ggsave("tests/coef_trajectories_stacked.png",
       p, width = 10, height = 14, dpi = 300)