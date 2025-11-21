library(ergm.rank)
library(coda)
library(statnet.common)
data(newcomb)

# Create the fit to simulate from if it does not already exist.
if(!file.exists("newcomb.2.fit.rds")) {
  newcomb.2.fit <- ergm(newcomb[[2]]~rank.deference+rank.nonconformity("all")+rank.nonconformity("localAND"), response="descrank", reference=~CompleteOrder)
  saveRDS(newcomb.2.fit, "newcomb.2.fit.rds")
} else newcomb.2.fit <- readRDS("newcomb.2.fit.rds")

# Prints time taken, effective sample size for each statistic, and a
# matrix of effective sample size per second (higher is better) for
# each combination.
print_results <- function(time, sim) {
  cat("Time taken (lower = better): \n")
  print(time)
  
  cat("\nAutocorrelation (lag 1, lower = better):\n")
  print(diag(drop(autocorr(sim[[1]], lags = 1))))

  cat("\nEffective sample sizes (higher = better):\n")
  print(effectiveSize(sim))

  cat("\nEfficiency (effective sample size / second of CPU time, higher = better):\n")  
  print(effectiveSize(sim) / time[1])
}

# Here, if a particular hint is not supported, it'll just give up and move on.
time_con <- function(constraints = ~.) {
  try({
    suppressWarnings(rm(time, sim)) # Don't accidentally use the previous run.
    time <- system.time(sim <- simulate(newcomb.2.fit, output = "stats", nsim = 1000, simplify = FALSE, constraints = constraints))
    cat("\n======================================\nHints: ", deparse1(constraints), "\n======================================\n\n")
    print_results(time, sim)
})
}

time_con() # default

time_con(~adjacent) # adjacent
