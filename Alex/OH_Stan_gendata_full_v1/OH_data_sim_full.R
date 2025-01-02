# Set constants
N_id <- 2        # Two individuals
N_board <- 10     # Number of boards in the sampling phase
N_food <- 2       # Two food item types
N_test <- 8 + 15  # Total test observations from both individuals

# Initialize y_sampling array
y_sampling <- array(0, dim = c(N_id, N_board, N_food))

# ID 1: Half the boards have 10 high-value, 0 low-value
#       Half the boards have 0 high-value, 10 low-value
# Let's say the first 5 boards are high-value only, and the next 5 are low-value only.

# First 5 boards: high-value = 10, low-value = 0
y_sampling[1, 1:5, 1] <- 10
y_sampling[1, 1:5, 2] <- 0

# Next 5 boards: high-value = 0, low-value = 10
y_sampling[1, 6:10, 1] <- 0
y_sampling[1, 6:10, 2] <- 10

# ID 2: Each board has 5 high-value and 5 low-value items
y_sampling[2, , 1] <- 5  # high-value items on each board for ID 1
y_sampling[2, , 2] <- 5  # low-value items on each board for ID 1

# Check the structure
y_sampling[1,,]
y_sampling[2,,]

# TEST PHASE DATA
# ID 1: 8 test observations
# ID 2: 15 test observations
# Combine them in a single vector:
id <- c(rep(1, 8), rep(2, 15))

# All observations are low-quality items (coded as 1, per the model definition)
y_test <- rep(1, N_test)

# Trial numbers should be per individual:
# For ID 1: trials 1 through 8
# For ID 2: trials 1 through 15
trial_num <- c(1:8, 1:15)

# Switch behavior:
# ID 1: no switch (1) for first 7, switch (2) on 8th
# ID 2: no switch (1) for first 14, switch (2) on 15th
Switch <- c(rep(1, 7), 2, rep(1, 14), 2)

# Put it all together into a list for Stan
dataList <- list(
  N_id = N_id,
  N_board = N_board,
  N_food = N_food,
  y_sampling = y_sampling,
  N_test = N_test,
  y_test = y_test,
  id = id,
  trial_num = trial_num,
  Switch = Switch
)

str(dataList)
