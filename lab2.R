#Lab 2

library(gbm)
library(tidyr)
library(ggplot2)
library(iml)
library(xgboost)

data<-readxl::read_xlsx("ML_data.xlsx")
names(data)

exposures<-data[,6:61]

scaled_exposures <- as.data.frame(scale(exposures))
colnames(scaled_exposures) <- paste0(colnames(exposures), "_scaled")
# Combine with original data if needed
data <- cbind(data, scaled_exposures)

############
## Gradient boosting
set.seed(123)


### Prepare data for GMB
data.gbm<-data
#data.gbm <- na.omit(data.gbm)

form1 <- as.formula(paste("bw ~", paste(names(data.gbm)[62:117],collapse="+")))

## Step 1: hyperparameters tuning
# Consider tuning parameters one by one and trial and error procedures ( see example in 12.3.3. https://bradleyboehmke.github.io/HOML/gbm.html)
# The following does the tuning for all parameters for simplicity

#
hyper_grid <- expand.grid(
  shrinkage = c(.001, .05, .1),
  interaction.depth = c(1,2,3),
  optimal_trees = 0,               # a place to dump results
  min_RMSE = 0                     # a place to dump results
)

# Generate a training data. Note, might not be necessary if doing 2 splits instead of 3 splits
random_index <- sample(1:nrow(data.gbm), nrow(data.gbm)*0.7)
random_ames_train <- data.gbm[random_index, ]

# Hyperparameters tuning for 3 splits -> more recommended for prediction

for(i in 1:nrow(hyper_grid)) {

  set.seed(123)
  gbm.tune <-gbm(formula = form1,
                 data=random_ames_train,
                 distribution = "gaussian",
                 n.trees = 1000,
                 interaction.depth = hyper_grid$interaction.depth[i],
                 shrinkage = hyper_grid$shrinkage[i],
                 n.cores = NULL,
                 train.fraction = .75, # refer to figure in slides, this is how much of overall training is training for first layer
                 verbose = FALSE)

  # add min training error and trees to grid
  hyper_grid$optimal_trees[i] <- which.min(gbm.tune$valid.error)
  hyper_grid$min_RMSE[i] <- sqrt(min(gbm.tune$valid.error))
}

hyper_grid %>%
  dplyr::arrange(min_RMSE) %>%
  head(10)

# Hyperparameters tuning for 2 splits -> more recommended for mechanisms description

random_index <- sample(1:nrow(data.gbm), nrow(data.gbm))
random_ames_train <- data.gbm[random_index, ]

for(i in 1:nrow(hyper_grid)) {

  set.seed(123)
  gbm.tune <-gbm(formula = form1,
                 data=random_ames_train,
                 distribution = "gaussian",
                 n.trees = 1000,
                 interaction.depth = hyper_grid$interaction.depth[i],
                 shrinkage = hyper_grid$shrinkage[i],
                 n.cores = NULL,
                 train.fraction = .75, # refer to figure in slides, this is how much of overall training is training for first layer
                 verbose = FALSE)

  # add min training error and trees to grid
  hyper_grid$optimal_trees[i] <- which.min(gbm.tune$valid.error)
  hyper_grid$min_RMSE[i] <- sqrt(min(gbm.tune$valid.error))
}

hyper_grid %>%
  dplyr::arrange(min_RMSE) %>%
  head(10)


# fine to add trees as long as shrinking is small.
# overfitting is mostly controlled by eta. n trees will make figures smoother

# For interpretability goals, also fine to select here depths=2 instead of depth=1


# Other hyperparameters
# n.minobsinnode = 10 by default
# bag.fraction = 0.50 by default

##  FINAL GBM model -> Depending on 2 vs 3 splits still use full training or full data. If doing prediction, use training and then predict in test
gbm1 <- gbm(
  formula = form1,
  data=data.gbm,
  distribution = "gaussian",
  n.trees = 1000,
  interaction.depth = 2,
  shrinkage = 0.05,
  train.fraction = 1,
  n.cores = NULL,
  verbose = FALSE)

## Interpretable ML
#1. VIP

### PLOT top 30
par(mar = c(5, 8, 1, 1))
summary(
  gbm1,
  cBars = 30,
  method = relative.influence, # alternative is permutation.test.gbm
  las = 2
)

vip::vip(gbm1)

# 2. PDP and ice. Example with p3

pdp_p3 <- pdp::partial(
  object = gbm1,
  pred.var = "p3_scaled",
  n.trees = gbm1$n.trees,
  grid.resolution = 100,
  train = data.gbm)

# Plot it
autoplot(pdp_p3, rug = TRUE, train=data.gbm)+
  geom_smooth(se = FALSE, method = "loess", span = 0.5)


ice2 <- gbm1 %>%
  pdp::partial(
    pred.var = "p3_scaled",
    n.trees = gbm1$n.trees,
    grid.resolution = 100,
    ice = TRUE
  ) %>%
  autoplot(rug = TRUE, train = random_ames_train, alpha = .1, center=TRUE) +
  ggtitle("Centered")
ice2

#Note, something with small importance, like p1 , will be more wiggly and less robust as present from fewer trees
ice2 <- gbm1 %>%
  pdp::partial(
    pred.var = "p1_scaled",
    n.trees = gbm1$n.trees,
    grid.resolution = 100,
    ice = TRUE
  ) %>%
  autoplot(rug = TRUE, train = random_ames_train, alpha = .1, center=TRUE) +
  ggtitle("Centered")
ice2

# 3. Interactions
# iml package implements model-agnostic procedures so requires prediction
predictor_gbm <- Predictor$new(
  model = gbm1,
  data = data.gbm[, all.vars(form1)],  # only use variables in the model
  y = data.gbm[[as.character(form1[[2]])]]  # extract response variable
)

# Compute interaction strengths
## These will take a few minutes to run

interact <- Interaction$new(predictor_gbm)

# View top interacting variables
interact$results %>%
  arrange(desc(.interaction)) %>%
  head()

plot(interact)

# These even longer
interact_2way <- Interaction$new(predictor_gbm, feature = "ph1_scaled")
interact_2way$results %>%
  arrange(desc(.interaction)) %>%
  top_n(10)

# Two-way PDP using iml
interaction_pdp <- Partial$new(
  predictor_gbm,
  c("ph1_scaled", "m5_scaled"),
  ice = FALSE,
  grid.size = 20
)
plot(interaction_pdp)

##########
#########
##########
#######
#####

##  XGBOOST

data.xgboost<-data.gbm[,c(2,62:117)]
val_ind <- sample.int(nrow(data.xgboost), 0.1 * nrow(data.xgboost))

x_train <- as.matrix(data.xgboost[-val_ind, !names(data.xgboost) %in% c("bw")])
x_label <- as.matrix(data.xgboost[-val_ind,"bw"])

x_val <- xgb.DMatrix(a<-as.matrix(data.xgboost[val_ind, !names(data.xgboost) %in% c("bw")]),
                     label = as.matrix(data.xgboost[val_ind,"bw"]))


# create hyperparameter grid
hyper_grid <- expand.grid(
  eta = c(.01, .05, .1),
  max_depth = c(1, 3),
  optimal_trees = 0,               # a place to dump results
  min_RMSE = 0                     # a place to dump results
)
for(i in 1:nrow(hyper_grid)) {

  # create parameter list
  params <- list(
    eta = hyper_grid$eta[i],
    max_depth = hyper_grid$max_depth[i]
  )

  # reproducibility
  set.seed(123)


  # train model
  xgb.tune  <- xgb.cv(
    params = params,
    data = x_train,
    label = x_label,
    nrounds = 1000,
    nfold = 5,
    objective = "reg:linear",  # for regression models
    verbose = 0,               # silent,
    early_stopping_rounds = 10 # stop if no improvement for 10 consecutive trees
  )

  # add min training error and trees to grid
  hyper_grid$optimal_trees[i] <- which.min(xgb.tune$evaluation_log$test_rmse_mean)
  hyper_grid$min_RMSE[i] <- min(xgb.tune$evaluation_log$test_rmse_mean)
}

hyper_grid %>%
  dplyr::arrange(min_RMSE) %>%
  head(10)

# untuned to make it quicker
#min_child_weight = 5,
#subsample = 0.65,
#colsample_bytree = 1

# parameter list
params <- list(
  eta = 0.01,
  max_depth = 3,
  min_child_weight = 5,
  subsample = 0.65,
  colsample_bytree = 1
)

# train final model
xgb.fit.final <- xgboost(
  params = params,
  data = x_train,
  label = x_label,
  nrounds = 500,
  objective = "reg:linear",
  verbose = 0
)


# predictors ranking

vip(xgb.fit.final)



### FINAL regression

final<-glm(bw~p3_scaled+ph1_scaled+p5_scaled+m11_scaled+p10_scaled+
             ph11_scaled+m9_scaled, data=data)
summary(final)$coefficients

# plus interaction and potential non-linearities

final_int<-glm(bw~p3_scaled+ph1_scaled+p5_scaled+m11_scaled+p10_scaled+
             ph11_scaled+m9_scaled+p5_scaled*ph1_scaled, data=data)
summary(final_int)$coefficients
