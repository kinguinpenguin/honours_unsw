start_time <- Sys.time()
library(ergm)
library(ergm.rank)
data(sampson)
data(samplk)
networks <- list(samplk1, samplk2, samplk3, sampdlk1, sampdlk2, sampdlk3)
complete_ordering <- function(nw, attrname, ref = as.matrix(nw, attrname = attrname)) { 
  diag(ref) <- NA
  r <- apply(ref, 1L, rank, ties.method = "random", na.last = "keep") |>
    t()
  diag(r) <- 0
  nw[,, names.eval=attrname, add.edges=TRUE] <- r
  nw
}
fit_list <- list()
for (i in 1:3) {
  cat(sprintf("Network %d\n", i))
  
  nw <- networks[[i]]
  
  # Complete ranking network used for MCMC proposals
  nw_ranked <- complete_ordering(nw, "score")
  
  # Ranking matrix matches the original observed network (not completed)
  rank_mat <- as.matrix(nw, attrname = "score")

  fit <- ergm(
    nw_ranked ~ rank.deference + rank.nonconformity("all") + rank.nonconformity("localAND"),
    response = "score",
    reference = ~CompleteOrder,
    constraints = ~ adjacent,
    obs.constraints = ~ ranking(rank_mat),
    control = snctrl(init.method = "zeros",
                    parallel = 4)
  )
  
  fit_list[[i]] <- fit
}

library(ggplot2)

compare_list <- list()

### ---- LIKE NETWORKS (1,2,3) ---- ###

for (i in 1:3) {
  if (inherits(fit_list[[i]], "ergm")) {

    coefs <- coef(fit_list[[i]])
    ses   <- sqrt(diag(vcov(fit_list[[i]])))
    logL  <- as.numeric(logLik(fit_list[[i]]))

    df_i <- data.frame(
      week = i,                   
      term = names(coefs),
      estimate = coefs,
      se = ses,
      logLik = logL
    )

    compare_list[[i]] <- df_i
  }
}

### ---- COMBINE ALL ---- ###

compare_df <- do.call(rbind, compare_list)
print(compare_df)
capture.output(print(compare_df), file = "tests/sampson/sampson_simulation.txt")

end_time <- Sys.time()
total_time <- as.numeric(difftime(end_time, start_time, units = "secs"))
cat("\nTotal computational time:", round(total_time, 2), "seconds\n",
    file = "tests/sampson/sampson_simulation.txt", append = TRUE)


### ---- SPLIT INTO LIKE & DISLIKE ---- ###

like_df <- do.call(rbind, compare_list[1:3])
like_df$lower <- like_df$estimate - like_df$se
like_df$upper <- like_df$estimate + like_df$se
like_df$week  <- factor(like_df$week)

plot_estimates <- function(df, title) {
  ggplot(df, aes(x = week, y = estimate)) +
    geom_point(size = 2) +
    geom_errorbar(aes(ymin = lower, ymax = upper),
                  width = 0.15) +
    facet_wrap(~ term, ncol = 1, scales = "free_y") +
    theme_minimal(base_size = 14) +
    labs(
      title = title,
      x = "Time Point",
      y = "Estimate"
    )
}

### ---- GENERATE AND SAVE FIGURES ---- ###

p_like <- plot_estimates(like_df, "Liking Relation Parameter Estimates")
ggsave("tests/sampson/likeness_statistics.png",
       p_like, width = 10, height = 14, dpi = 300)