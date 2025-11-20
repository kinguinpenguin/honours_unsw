start_time <- Sys.time()
library(ergm)
library(ergm.rank)
data(sampson)
data(samplk)
networks <- list(sampdlk1, sampdlk2, sampdlk3)
complete_ordering <- function(nw, attrname, ref = as.matrix(nw, attrname = attrname)) { 
  diag(ref) <- NA
  r <- apply(ref, 1L, rank, ties.method = "random", na.last = "keep") |>
    t()
  diag(r) <- 0
  nw[,, names.eval=attrname, add.edges=TRUE] <- r
  nw
}
fit_list <- list()
for (i in 1:1) {
  cat(sprintf("Network %d\n", i))
  nw <- networks[[i]]
  
  # Complete ranking network
  nw_ranked <- complete_ordering(nw, "score")
  
  # Ranking matrix matches the completed network
  rank_mat <- as.matrix(nw, attrname = "score")
  
  fit <- ergm(
    nw_ranked ~ rank.deference + rank.nonconformity("all") + rank.nonconformity("localAND"),
    response = "score",
    reference = ~CompleteOrder,
    constraints = ~ adjacent,
    obs.constraints = ~ ranking(rank_mat),
    control = snctrl(init.method = "zeros")
  )
  
  fit_list[[i]] <- fit
}
output_file <- "tests/sampson/sampson_simulation_dislike.txt"

# Open the file for writing
sink(output_file)

cat("Sampson Simulation Results\n")
cat("==========================\n\n")

for (i in seq_along(fit_list)) {
  cat(sprintf("Model %d Results\n", i))
  cat("-----------------\n")
  
  # Print summary of the fit
  print(summary(fit_list[[i]]))
  
  cat("\n\n")  # spacing between blocks
}

# Close the sink
sink()

library(ggplot2)

compare_list <- list()


### ---- DISLIKE NETWORKS (4,5,6) ---- ###

for (i in 1:1) {
  if (inherits(fit_list[[i]], "ergm")) {

    coefs <- coef(fit_list[[i]])
    ses   <- sqrt(diag(vcov(fit_list[[i]])))
    logL  <- as.numeric(logLik(fit_list[[i]]))

    df_i <- data.frame(
      week = i,     # 1,2,3 for dislike
      term = names(coefs),
      estimate = coefs,
      se = ses,
      logLik = logL
    )

    compare_list[[i]] <- df_i
  }
}
for (i in 2:3) {
  if (inherits(fit_list[[1]], "ergm")) {

    coefs <- coef(fit_list[[1]])
    ses   <- sqrt(diag(vcov(fit_list[[1]])))
    logL  <- as.numeric(logLik(fit_list[[1]]))

    df_i <- data.frame(
      week = i,     # 1,2,3 for dislike
      term = names(coefs),
      estimate = 0,
      se = 0,
      logLik = 0
    )

    compare_list[[i]] <- df_i
  }
}


### ---- COMBINE ALL ---- ###

compare_df <- do.call(rbind, compare_list)
print(compare_df)

dislike_df <- do.call(rbind, compare_list[4:6])
dislike_df$lower <- dislike_df$estimate - dislike_df$se
dislike_df$upper <- dislike_df$estimate + dislike_df$se
dislike_df$week  <- factor(dislike_df$week)


### ---- FIXED PLOTTING FUNCTION (NO model_type) ---- ###

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
p_dislike <- plot_estimates(dislike_df, "Disliking Relation Parameter Estimates")
ggsave("tests/sampson/dislikeness_statistics.png",
       p_dislike, width = 10, height = 14, dpi = 300)