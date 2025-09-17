
## LAB 1. ML in epi / September 2025

# 1.	Prepare the dataset

# Packages
#devtools::install_github("laresbernardo/lares")
# Required packages and open data
Packages <- c("readxl", "qgraph", "corrplot", "factoextra","dplyr","naniar",
              "ggplot2","gridExtra","lares","ggrepel","rms","glmnet","caret",
              "glmnetUtils","gbm","xgboost","vip","pdp","iml")
lapply(Packages, library, character.only = TRUE)

data<-readxl::read_xlsx("ML_data.xlsx")
names(data)

### 1 binary outcome, 1 continuous, 16 pollutants, 22 phthalates, 18 metals

# Pre-processing

exposures<-data[,6:61]

### Missingness

gg_miss_var(data)


# 2.	Investigate the data
### Distributions


#### Outcomes
summary(data$bw)
hist(data$bw)
data$lbw<-ifelse(data$bw<2500,1,0)
table(data$lbw)

#### Exposures


p1_plot <- ggplot(data, aes(x = p1)) +
  geom_histogram( fill = "steelblue", color = "black") +
  ggtitle("p1")

ph1_plot <- ggplot(data, aes(x = ph1)) +
  geom_histogram(fill = "darkgreen", color = "black") +
  ggtitle("ph1")

m1_plot <- ggplot(data, aes(x = m1)) +
  geom_histogram(fill = "purple", color = "black") +
  ggtitle("m1")

ph5_plot <- ggplot(data, aes(x = ph5)) +
  geom_histogram( fill = "orange", color = "black") +
  ggtitle("ph5")

grid.arrange(p1_plot, ph1_plot, m1_plot, ph5_plot, ncol = 2)

# will need scaling

scaled_exposures <- as.data.frame(scale(exposures))
colnames(scaled_exposures) <- paste0(colnames(exposures), "_scaled")
# Combine with original data if needed
data <- cbind(data, scaled_exposures)


p1_plot <- ggplot(data, aes(x = p1_scaled)) +
  geom_histogram( fill = "steelblue", color = "black") +
  ggtitle("p1")

ph1_plot <- ggplot(data, aes(x = ph1_scaled)) +
  geom_histogram(fill = "darkgreen", color = "black") +
  ggtitle("ph1")

m1_plot <- ggplot(data, aes(x = m1_scaled)) +
  geom_histogram(fill = "purple", color = "black") +
  ggtitle("m1")

ph5_plot <- ggplot(data, aes(x = ph5_scaled)) +
  geom_histogram( fill = "orange", color = "black") +
  ggtitle("ph5")

grid.arrange(p1_plot, ph1_plot, m1_plot, ph5_plot, ncol = 2)

### Correlations

cor.matrix <- cor(exposures, method="spearman")

corrplot(cor.matrix,
         method="circle",
         order = "hclust",
         addrect =10,
         tl.pos = "l",
         tl.col = "black", tl.cex=0.95, col=COL2("RdYlBu"))

### Network

qgraph(cor.matrix, graph = "pcor", layout = "spring")

a<-qgraph(cor.matrix, layout="spring", borders=FALSE, minimum = 0.5,
          edge.labels=FALSE, labels=names(cor.matrix),
          label.scale=FALSE,label.cex=.5,edge.width=.5)

# highest
a<-qgraph(cor.matrix, layout="spring", borders=FALSE, minimum = 0.8,
          edge.labels=FALSE, labels=names(cor.matrix),
          label.scale=FALSE,label.cex=.5,edge.width=.5)

corr<-as.vector(a$Edgelist$weight)

## how many
length(corr[corr>0.9])
## [1] 69
length(corr[corr>0.8])
## [1] 256
## higher correlations
corr_cross(data[,4:59],top = 20)

## Extend to checkas within categories

##########
##########
# 3.	Check univariate associations

# Regression / univariate
## function for plotting volcano plot Bonferroni
volcano <- function(data) {
  data$Legend <- "FDR>=0.05"
  data$Legend[data$fdr < 0.05] <- "FDR<0.05 "

  data$delabel <- NA
  data$delabel[data$Legend != "FDR>=0.05"] <- data$biomarker[data$Legend != "FDR>=0.05"]

  ggplot(data=data, aes(x=or, y=-log10(pvalue), col=Legend, label=delabel)) +
    geom_point() +
    theme_minimal() + scale_y_continuous(limits=c(0,10))+
    geom_text_repel()+scale_x_continuous(trans='log')+labs(x="OR (per 1 SD)")+
    theme(axis.text.x = element_blank())
}

volcano_linear <- function(data) {
  data$Legend <- "FDR>=0.05"
  data$Legend[data$fdr < 0.05] <- "FDR<0.05"

  data$delabel <- NA
  data$delabel[data$Legend != "FDR>=0.05"] <- data$biomarker[data$Legend != "FDR>=0.05"]

  ggplot(data = data, aes(x = beta, y = -log10(pvalue), col = Legend, label = delabel)) +
    geom_point() +
    theme_minimal() +
    scale_y_continuous(limits = c(0, 10)) +
    geom_text_repel() +
    labs(x = "Beta (per 1 SD)") +
    theme(axis.text.x = element_blank())
}


## Function for unadjusted model [logistic regression here: this can be modified if the endpoint is time-to-event ]

multiple<- function(data,cc,start,end){
  data <- as.data.frame(data)
  exp_nvar=end-start+1
  biomarker=rep(NA, exp_nvar)
  or=rep(NA, exp_nvar)
  pvalue=rep(NA, exp_nvar)
  bonferroni=rep(NA, exp_nvar)
  controls=rep(NA, exp_nvar)
  cases=rep(NA, exp_nvar)
  number=1
  data$cc <- data[ , cc]
  for (j in start:end){
    exposure = colnames(data)[j]
    model <- glm(cc~get(exposure),
                 data=data,family="binomial")
    or[number] = as.numeric(round(exp(summary(model)$coefficients[2,1]),3))
    pvalue[number] = as.numeric(summary(model)$coefficients[2,4])
    biomarker[number] = exposure
    controls[number] = table(model$model$cc)[1]
    cases[number] = table(model$model$cc)[2]
    number = number + 1
  }
  results <- data.frame(biomarker, cases, controls, or, pvalue)
  results<<-results

}


multiple_linear <- function(data, outcome_col, start, end) {
  data <- as.data.frame(data)
  exp_nvar <- end - start + 1

  biomarker <- rep(NA, exp_nvar)
  beta <- rep(NA, exp_nvar)
  pvalue <- rep(NA, exp_nvar)
  bonferroni <- rep(NA, exp_nvar)
  n <- rep(NA, exp_nvar)

  data$outcome <- data[[outcome_col]]
  number <- 1

  for (j in start:end) {
    exposure <- colnames(data)[j]
    model <- lm(outcome ~ get(exposure), data = data)

    beta[number] <- round(summary(model)$coefficients[2, 1], 3)
    pvalue[number] <- summary(model)$coefficients[2, 4]
    biomarker[number] <- exposure
    n[number] <- nrow(model$model)

    number <- number + 1
  }

  results <- data.frame(biomarker, n, beta, pvalue)
  results$fdr <- p.adjust(results$pvalue, method = "fdr")
  results<<-results
}

# continuous outcome
multiple_linear(data=data,outcome_col="bw",start=63,end=118)
results
results$bonferroni= p.adjust(results$pvalue, method = "bonferroni")
results1 <-results[order(results$fdr),]
results1
volcano_linear(results)

# binary outcome

multiple(data=data,cc="lbw",start=63,end=118)
results$bonferroni= p.adjust(results$pvalue, method = "bonferroni")
results$fdr= p.adjust(results$pvalue, method = "fdr")
results1 <-results[order(results$bonferroni),]
results1
volcano(results)


corr_var(exposures, p5, top = 20)

# 4.	Non-Linearity
## Non-linearity checks

nonlinear_cont <- function(data, outcome_col, start, end) {
  data <- as.data.frame(data)
  exp_nvar <- end - start + 1

  biomarker <- rep(NA, exp_nvar)
  pvalue <- rep(NA, exp_nvar)
  number <- 1

  data$outcome <- data[[outcome_col]]

  for (j in start:end) {
    exposure <- colnames(data)[j]

    # Fit linear model with restricted cubic spline
    model <- lm(outcome ~ rcs(get(exposure), 3), data = data)

    # Extract p-value for non-linear term (3rd coefficient)
    pvalue[number] <- summary(model)$coefficients[3, 4]
    biomarker[number] <- exposure

    number <- number + 1
  }

  results <- data.frame(biomarker, pvalue)
  results <- results[order(results$pvalue), ]

  # Select significant biomarkers
  significant <- as.vector(results$biomarker[results$pvalue < 0.05])

  return(list(results = results, significant = significant))
}


nonlinear_cont(data=data,outcome_col="bw",start=63,end=118)


## Splines, one here for providing example

p9rcs <- glm(bw ~ rcs(p9_scaled, 4),
             data = data)
new_data <- data.frame(p9_scaled = seq(min(data$p9_scaled), 5, length.out = 500))

# Predict the fitted response and standard error on the new data
predictions <- predict(p9rcs, newdata = new_data, se.fit = TRUE)

# Add predictions and confidence intervals to the new data frame
new_data$fit <- predictions$fit
new_data$se.fit <- predictions$se.fit
new_data$lwr <- new_data$fit - 1.96 * new_data$se.fit
new_data$upr <- new_data$fit + 1.96 * new_data$se.fit

ggplot(data = data, aes(x = p9_scaled, y = bw)) +
  # Add the fitted spline curve
  geom_line(data = new_data, aes(y = fit), color = "blue", linewidth = 1) +
  # Add the 95% confidence interval as a shaded ribbon
  geom_ribbon(data = new_data, aes(y = fit, ymin = lwr, ymax = upr), fill = "blue", alpha = 0.2) +
  # Customize plot aesthetics
  labs(
    title = "p9_scaled",
    x = "p9_scaled",
    y = "BW"
  ) +
  theme_minimal()

# briefly discuss functions and paper for binary and time to event
source("splines_functions.R")


# 5.	Multivariable regression
##  multiple regression
form1 <- as.formula(paste("bw ~", paste(names(data)[63:118],collapse="+")))

a<-glm(form1, data=data)
a$coefficients
vif(a)
###

# add multiple within categories

###
# 6. Penalized regression
set.seed(123)

X<-as.matrix(data[,6:61])
Y<-data$bw
lambda_grid <- 10^seq(-3, 4, length.out = 100)

?glmnet
# https://cran.r-project.org/web/packages/glmnet/vignettes/glmnet.pdf
# Ridge

# Note that direction of plots has been changed in the latest version
?plot.glmnet

res <- glmnet(X, Y, alpha = 0, lambda = lambda_grid, standardize = TRUE)

res$lambda

plot(res, xvar = "lambda", xlab="-log(lambda)")
legend("bottomright", lwd = 1, legend = colnames(X), cex = .7)

# note type.measure = "mse" by default
ridge_cv <- cv.glmnet(X, Y, alpha = 0,
                      lambda = lambda_grid,
                      standardize = TRUE,
                      nfolds = 1000)


plot(ridge_cv, xlab="-log(lambda)")

# “lambda.min”: the λ at which the smallest MSE is achieved.
# “lambda.1se”: the largest λ at which the MSE is within one standard error of the smallest MSE (default).

# lowest lambda
lambda_cv_min <- ridge_cv$lambda.min
# Best cross-validated lambda
lambda_cv <- ridge_cv$lambda.1se


model_cv <- glmnet(X, Y, alpha = 0,
                   lambda = lambda_cv, standardize = TRUE)


res<-as.matrix(round(model_cv$beta,3))
colnames(res)<-"Estimate"
res
rownames(res)[res[, "Estimate"] != 0]

# LASSO

res_lasso <- glmnet(X, Y, alpha = 1, lambda = lambda_grid, standardize = TRUE)
plot(res_lasso, xvar = "lambda", xlab="-log(lambda)")
legend("bottomright", lwd = 1, legend = colnames(X), cex = .7)

lasso_cv <- cv.glmnet(X, Y, alpha = 1,
                      lambda = lambda_grid,
                      standardize = TRUE,
                      nfolds = 1000)
plot(lasso_cv, xlab="-log(lambda)")

# lowest lambda
lambda_cv_min_lasso <- lasso_cv$lambda.min
# Best cross-validated lambda
lambda_cv_lasso <- lasso_cv$lambda.1se

model_cv_lasso <- glmnet(X, Y, alpha = 1,
                         lambda = lambda_cv_lasso,
                         standardize = TRUE)


res<-as.matrix(round(model_cv_lasso$beta,3))
colnames(res)<-"Estimate"
res
rownames(res)[res[, "Estimate"] != 0]

model_cv_lasso2 <- glmnet(X, Y, alpha = 1,
                          lambda = lambda_cv_min_lasso,
                          standardize = TRUE)


res2<-as.matrix(round(model_cv_lasso2$beta,3))
colnames(res)<-"Estimate"
res2
rownames(res2)[res2[, "Estimate"] != 0]


# Elastic Net

# example with random alpha
enet_cv <- cv.glmnet(X, Y, alpha = 0.3,
                     lambda = lambda_grid,
                     standardize = TRUE,
                     nfolds = 1000)

# Best cross-validated lambda
lambda_cv_enet <- enet_cv$lambda.1se

model_cv_enet <- glmnet(X, Y, alpha = 0.1,
                        lambda = lambda_cv_enet,
                        standardize = TRUE)


res<-as.matrix(round(model_cv_enet$beta,3))
colnames(res)<-"Estimate"
res
rownames(res)[res[, "Estimate"] != 0]



# glmnet does not support tuning together alpha and lambda. Use glmnetUtils

# https://cran.r-project.org/web/packages/glmnetUtils/vignettes/intro.html

# Perform cross-validation for both alpha and lambda
# cva.glmnet automatically handles a range of alpha and lambda values
cv_elastic_net <- cva.glmnet(X, Y, nfolds = 10)
plot(cv_elastic_net)


#min MSE at varying alphas are similar. slightly better with smaller one, but we can decide how much we want to
# Be conservative. The goal here is not prediction but explaining



lambda_cv_enet2.min<-cv_elastic_net$modlist[[which.min(abs(cv_elastic_net$alpha - 0.1))]]$lambda.min
lambda_cv_enet2<-cv_elastic_net$modlist[[which.min(abs(cv_elastic_net$alpha - 0.1))]]$lambda.1se

model_cv_enet2 <- glmnet(X, Y, alpha = 0.1,
                         lambda = lambda_cv_enet2,
                         standardize = TRUE)


res2<-as.matrix(round(model_cv_enet2$beta,3))
colnames(res2)<-"Estimate"
res2
rownames(res2)[res2[, "Estimate"] != 0]

model_cv_enet2 <- glmnet(X, Y, alpha = 0.1,
                         lambda = lambda_cv_enet2.min,
                         standardize = TRUE)


res2<-as.matrix(round(model_cv_enet2$beta,3))
colnames(res2)<-"Estimate"
res2
rownames(res2)[res2[, "Estimate"] != 0]

# Example: how to adjust for covariates

X_conf<-as.matrix(data[,c(3:5,63:118)])

enet_cv_adj <- cv.glmnet(X_conf, Y, alpha = 0.1,
                         lambda = lambda_grid,
                         standardize = FALSE, nfolds = 1000,
                         penalty.factor=c(0,0,0,rep(1,ncol(X_conf)-3)))



# lowest lambda
lambda_cv_min_enet_adj <- enet_cv_adj$lambda.min
# Best cross-validated lambda
lambda_cv_enet_adj <- enet_cv_adj$lambda.1se

model_cv_enet_adj <- glmnet(X_conf, Y, alpha = 0.1,
                            lambda = lambda_cv_enet_adj,
                            standardize = TRUE,
                            penalty.factor=c(0,0,0,rep(1,ncol(X_conf)-3)))



res<-as.matrix(round(model_cv_enet_adj$beta,3))
colnames(res)<-"Estimate"
res
rownames(res)[res[, "Estimate"] != 0]


# family = "cox", poisson, multinomial, binomial. if cox, Y is Surv


### Extend to checks within families of ph, PM, metals
