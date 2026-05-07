library(terra)
library(rsample)
library(tidyverse)
library(sf)
library(caret)

# data splitting examples

## rsample method

# 1. Sample SpatVector (e.g., points)
pts <- vect(matrix(runif(200), 100, 2), crs="EPSG:4326")
pts$class <- rep(c("A", "B"), 50) # Example classification outcome

# Convert to sf for tidy manipulation
pts_sf <- st_as_sf(pts)

# 2. Stratified 3-way split (Train/Validation/Test)
set.seed(123)
data_split <- initial_validation_split(pts_sf, strata = class, prop = c(0.7, 0.15))

# Extract subsets
train_sf <- training(data_split)
val_sf   <- validation(data_split)
test_sf  <- testing(data_split)

# 3. Convert back to SpatVector
train_vect <- vect(train_sf)
val_vect   <- vect(val_sf)
test_vect  <- vect(test_sf)


#### Caret Method

set.seed(123)
# 1. First split: 70% training, 30% for others (validation + test)
train_idx <- createDataPartition(pts_sf$class, p = 0.7, list = FALSE)
train_sf  <- pts_sf[train_idx, ]
others_sf <- pts_sf[-train_idx, ]
# 2. Second split: Split the remaining 30% into two equal parts (15% each)
val_idx <- createDataPartition(others_sf$class, p = 0.5, list = FALSE)
val_sf  <- others_sf[val_idx, ]
test_sf <- others_sf[-val_idx, ]

# 3. Convert back to SpatVector
train_vect <- vect(train_sf)
val_vect   <- vect(val_sf)
test_vect  <- vect(test_sf)

