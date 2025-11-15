start_time <- Sys.time()
library(ergm.rank)
data("newcomb")
fit_true <- list()
fit_constrained <- list()
for (week in 1:15) {
  fit2.true <- ergm(newcomb[[week]] ~ rank.deference+rank.nonconformity("all")+
                      rank.nonconformity("localAND"),
              response="descrank",
              reference=~CompleteOrder, # Sample space: make adjacent swap proposals.
              constraints = ~ adjacent,
              control = snctrl()
  )
  fit_true[[week]] <- fit2.true
}

end_time <- Sys.time()
total_time <- as.numeric(difftime(end_time, start_time, units = "secs"))
cat(sprintf("Total time: %.2f seconds\n", total_time))