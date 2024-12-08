## optimal K and C for 1. stock ADS
set.seed(2024)

## packages
install.packages("foreach")
install.packages("doParallel")
library(foreach)
library(doParallel)

## please adjust wd
setwd("~/Max Thesis")

#### load data, please adjust the paths ####
historical_returns_data <- read.csv("C:\\Users\\Max\\Documents\\Max Thesis\\DAX40_historical_returns.csv")
dirft_data <- read.csv("C:\\Users\\Max\\Documents\\Max Thesis\\DAX40_drift.csv")
volatility_data <- read.csv("C:\\Users\\Max\\Documents\\Max Thesis\\DAX40_volatility.csv")
DAX40_data <- read.csv("C:\\Users\\Max\\Documents\\Max Thesis\\DAX40_Data.csv")

## extract drift, volatility and historcial returns for the stock
stock_name <- "ADS"
historical_returns <- historical_returns_data[[stock_name]]
drift_value <- dirft_data[dirft_data[, 1] == stock_name, 2]
volatility_value <- volatility_data[volatility_data[, 1] == stock_name, 2] 
S0_value <- DAX40_data[1, stock_name]

#### set the bins ####
## max and min
range_min <- drift_value - 3 * volatility_value
range_max <- drift_value + 3 * volatility_value

##number of bins
num_bins <- 21

## generating the bins: (+1 bc of the edges, one more edge then there is bins)
bins <- seq(range_min, range_max, length.out = num_bins + 1)



#### extract frequencies and bin creation####
## check how many point are outisde the defined range
outside_points <- historical_returns[historical_returns < range_min | historical_returns > range_max]
## exclude that one
historical_returns <- historical_returns[historical_returns >= range_min & historical_returns <= range_max]
## histogram
observed_hist <- hist(historical_returns, breaks = bins, plot = FALSE)
observed_frequency <- observed_hist$counts

## now lets check if there are bins filled with less than 1% of all observations
min_frequency <- 0.01 * length(historical_returns)
low_frequency_bins <- which(observed_frequency< min_frequency)

## Adjust bins with low frequencies
if (length(low_frequency_bins) > 0) {
  while (length(low_frequency_bins) > 0) {
    low_bin <- low_frequency_bins[1]
    #if its last bin merge with previous
    if (low_bin == length(observed_frequency)) {
      observed_frequency[low_bin -1] <- observed_frequency[low_bin - 1] + observed_frequency[low_bin]
    }
    # if not merge with next bin
    else {
      observed_frequency[low_bin + 1] <- observed_frequency[low_bin] + observed_frequency[low_bin + 1]
    }
    # remove low bin
    observed_frequency <- observed_frequency[-low_bin]
    bins <- bins[-(low_bin + 1)]
    low_frequency_bins <- which(observed_frequency < min_frequency)
  }
}
# double check
if (any(observed_frequency == 0)) {
  print("There are empty bins!!! - check")
} else {
  print("No empty bins - all good!")
}



#### create IFBM function ####
IFBM_loop <- function(nsim, t, mu, sigma, S0, k, c, dt) {
  # storage
  IFBM <- matrix(ncol = nsim, nrow = t)
  for (simu in 1:nsim) {
    IFBM[1, simu] <- S0
    for (day in 2:t) {
      epsilon <- rnorm(1) 
      # calculate the different parts
      f_epsilon <- (2*exp(-c * (epsilon^2 / 2)) - 1) * atan(epsilon)
      drift <- (mu - (sigma^2 / 2)) * dt
      volatility <- sigma * epsilon * sqrt(dt)
      modification <- (mu - (sigma^2 / 2)) * k * f_epsilon * dt
      # simulation
      IFBM[day, simu] <- IFBM[(day - 1), simu] * exp(drift + volatility + modification)
    }
  }
  return(IFBM)
}

#### run IFBM simulations ####
set.seed(2024)

# set parameters
k_values <- seq(-2, 2, by = 0.1)
c_values <- seq(-2, 2, by = 0.1)
nsim <- 1000
t <- 132
dt <- 1
S0 <- S0_value
mu <- drift_value
sigma <- volatility_value

# clusters
number_cores <- detectCores() - 1
clusters <- makeCluster(number_cores)
# enables foreach package to distributes simulations:
registerDoParallel(clusters)

# this is just for the clusters to stop if there is an error - pc got overloaded on some runs
on.exit({
  stopCluster(clusters)
  registerDoSEQ() 
})

# multiple simulation at once of avge prices
avg_results <- foreach(k = k_values, .combine = 'list', .packages= "base") %:%
  foreach(c = c_values, .combine ='list', .packages = "base") %dopar% {
    # unique seed for each k and c combination - had to be added here otherwise values differed
    set.seed(2024 + as.integer(k * 100 + c * 10))  
    IFBM_simulations <- IFBM_loop(nsim, t, mu, sigma, S0, k, c, dt)
    avg_prices <- rowMeans(IFBM_simulations)
    result <- list()
    result[[paste0("k=", k, "_c=", c)]] <- avg_prices
    return(result)
  }


# make a single list due to the nested loops above
avg_results <- do.call(c, avg_results)
# still had more layers, function to flatten:
flatten_list <- function(x) {
  if (is.list(x) && all(sapply(x, is.list))) {
    do.call(c, lapply(x, flatten_list))
  } else {
    x
  }
}

# apply
avg_results <- flatten_list(avg_results)


#### calculate returns and bin the frequencies####
## returns
# list
avg_returns <- list()
#loop
for (x in names(avg_results)) {
  avg_prices <- avg_results[[x]]
  returns <- diff(avg_prices)/head(avg_prices, -1)
  avg_returns[[x]] <-returns
}


##frequencies (here are the parts were i have some questions)
#list
simulated_frequencies_list <- list()
# loop
for (x in names(avg_returns)) {
  sim_returns <- avg_returns[[x]]
  # filter out returns which were outsider the bin
  sim_returns <- sim_returns[sim_returns >= min(bins) & sim_returns <= max(bins)]
  # sort returns which are in the range to bins
  sim_hist <- hist(sim_returns, breaks = bins, plot = FALSE)
  simulated_frequencies_list[[x]] <- sim_hist$counts
}

# here i filter out simulations were not at leasdt 90% of the returns fall into the predefined bins
# to ensure that there are no low chi squared values only bc of low data coverage
# or too few observations in the bins which could give wrong results

# 90%
threshold <- 0.9

# list
valid_combinations <- list()

# filter out
for (x in names(simulated_frequencies_list)) {
  sim_frequencies <- simulated_frequencies_list[[x]]
  total_simulated <- length(avg_returns[[x]]) 
  in_bins <- sum(sim_frequencies) 
  ratio <- in_bins/total_simulated 
  if (ratio >=threshold) {
    valid_combinations[[x]] <- sim_frequencies
  }
}

#### chi squared test ####
#df
chi_squared_results <- data.frame(
  Combinations = character(0),
  TotalSimulated = numeric(0),
  ChiSquared = numeric(0),
  PValue = numeric(0)
)

# loop chi squared test
for (x in names(valid_combinations)) {
  sim_frequencies <- valid_combinations[[x]]
  # not sure if scaling should be done
  scaled_observed <- observed_frequency / sum(observed_frequency) * sum(sim_frequencies)
  chi_test <- chisq.test(sim_frequencies, p = scaled_observed / sum(scaled_observed))
  # store
  chi_squared_results <- rbind(chi_squared_results, data.frame(
    Combination = x,
    TotalSimulated = sum(sim_frequencies),
    ChiSquared = chi_test$statistic,
    PValue = chi_test$p.value
  ))
}


chi_squared_results <- chi_squared_results[order(chi_squared_results$ChiSquared), ]
print(head(chi_squared_results))


