# Script associated with the paper: "Impact-based early warning of mass movements – A dynamic spatial modelling approach for the Alpine region"
# doi: 

# The R code uses data from the repository to fit and analyze the models and to visualize data and results.
# It also contains code snippets on how certain anlayses (e.g. basin-based cross validation, variable importance) were performed
#
## Structure of the script
# Step 1: Load data and basic visualizations
# Step 2: Fit the models
# Step 3: Calculate fitting performance and plot ROC curves (Fig. 4a in publication)
# Step 4: Cross validation: Snippet code for basin-based cross validation and visualizing final performance (Fig. 4a-d in publication)
# Step 5: Variable importance plots (incl. snippet code for calculating feature importance) (Fig. 5)
# Step 6: Visualizing partial effects from fitted models (Fig. 6, Fig. 7, Fig. 8)

rm(list=ls()); gc() # clean workspace and run garbage collection to free memory
# Load packages ---------------------------------------------------------------------------------------
library(mgcv)      # fitting GAMs (bam/gam)
library(sf)        # spatial vector data (simple features)
library(pROC)      # ROC / AUC calculations
library(gratia)    # diagnostic & plotting helpers for GAMs (draw(), smooths(), etc.)
library(sperrorest) # spatial cross-validation helper (partition_factor_cv)
library(ggplot2)   # plotting
library(tmap)      # thematic mapping for sf objects
library(dplyr)     # data manipulation (select, mutate, %>%)
library(patchwork) # combine ggplot objects
library(vip)       # variable importance (vi_permute)
#
##
###
##########
########## Step 1: Load data ----------------------------------------------------------------------------------
##########
###
##
#
folder <- "/export/home/klidl/ssteger/scratch/xriskcc/CodeDataSharing" # Set Path to folder where data is stored

# Load each object (readRDS used for reproducible storage of R objects) --------------------------------
basins     <- readRDS(file.path(folder, "basins_compressed.rds"))          # Basins with static data (sf object) -- polygons with attributes
AS         <- readRDS(file.path(folder, "Alpine_Space_Boundary_4326.rds")) # Alpine Space outline (sf)
df_slides  <- readRDS(file.path(folder, "df_slides.rds"))                  # Training data for slide-type (SL) model
df_flows   <- readRDS(file.path(folder, "df_flows.rds"))                   # Training data for flow-type (DF) model
df_falls   <- readRDS(file.path(folder, "df_falls.rds"))                   # Training data for fall-type (RF) model
SLcv       <- readRDS(file.path(folder, "mySLcv_5f_10r.rds"))              # Cross validation results for slide type (SL) model
DFcv       <- readRDS(file.path(folder, "myDFcv_5f_10r.rds"))              # Cross validation results for flow type (DF) model
RFcv       <- readRDS(file.path(folder, "myRFcv_5f_10r.rds"))              # Cross validation results for fall type (RF) model
varimpo_raw    <- readRDS(file.path(folder, "varimpo_raw.rds"))            # Results (raw) from variable importance assessment (raw permute outputs)
varimpoSL <- readRDS(file.path(folder, "varimpoSL.rds"))                   # Plot-ready variable importance data for slide-type (SL)
varimpoDF <- readRDS(file.path(folder, "varimpoDF.rds"))                   # Plot-ready variable importance data for flow-type (DF)
varimpoRF <- readRDS(file.path(folder, "varimpoRF.rds"))                   # Plot-ready variable importance data for fall-type (RF)

# Colors for plots -----------------------------------------------------------------------------------
cols = "#56b899" # color slide-type models (SL)
cold = "#9052b6" # color flows-type  models (DF)
colr = "#e06c3e" # color fall-type models (RF)

# Example: Alpine Space outline map --------------------------------------------
tm_shape(AS) + tm_fill(col = "lightgrey") + tm_borders(col = "black", lwd = 2) + tm_title("Alpine Space Boundary")

# Example plot for basins: Mean slope angles in the respective PPAs -----------------------------------
# "Missing" (NA) means no PPA pixel in the entire basin; facets show SL, DF, RF slope maps
# Basin polygons are filled without outlines; use tm_polygons() for visualizing borders.
tm_shape(basins) +
  tm_fill(c("SL_Slope", "DF_Slope", "RF_Slope"), # Slides (SL), Flows (DF), Falls (RF)  -- plot three attributes
    fill.scale = tm_scale(values = c("#d7f4be", "#f0d076", "#e73030")),
    fill.legend = tm_legend(title = "Mean slope in PPA"), fill.free = FALSE) +
  tm_facets(ncol = 1)
##
###
##########
########## Step 2: Fit the models ----------------------------------------------------------------------------------
##########
###
##
#
maxk <- 4  # k-factor to restrict maximum flexibility of smooth terms (keeps smooths conservative/generalized)

# Define formula for slide-type (SL) ------------------------------------------------------------------
fo_slides = 
SL01 ~                              # Binary response variable for slide-type (SL) -- dependent variable in training data
s(tp_2,                   k=maxk) + # Short-term precipitation (smoothed)
s(tp_30,                  k=maxk) + # Antecedent precipitation: 30-day (smoothed)
s(tp_su_mean_annual_mean, k=maxk) + # Mean annual precipitation (smoothed)
s(SL_Slope ,              k=maxk) + # Slope angle in slide-type PPA (smoothed)
SL_dominant_landcover             + # Dominant land cover class in slide-type PPA (parametric categorical)
SL_dominant_lithology             + # Dominant lithology class in slide-type PPA (parametric categorical)
s(log_B_SL_n ,            k=maxk) + # Log Number of buildings in slide-type PPA (smoothed)
s(SL_roads_12 ,           k=maxk) + # Portion of “road” cells in slide PPA (smoothed)
s(doy , bs="cc",          k=maxk) + # Day-of-Year as circular variable (cyclic spline)
s(temp_mean,              k=maxk) + # Mean daily temperature (smoothed)
temp_sub0                         + # Binary: Day below 0 °C (yes/no) (parametric)
s(cat,  bs="re")                  + # Random effect: Sampling location (Basin-ID) using bs="re"
s(year, bs="re")                    # Random effect: Sampling Year using bs="re"

# Fit the slide type model ----------------------------------------------------------------------------
fit_slides = mgcv::bam(fo_slides, data=df_slides, family=binomial, discrete = 100, method="fREML", select=T) # select=T enables automatic smooth selection
summary(fit_slides) # check model summary
plot(fit_slides, pages = 1) # quick partial effect plots (exploratory; use gratia::draw for publication-quality plots --> see Step: 6)

# Define formula for flow types (DF) -------------------------------------------------------------------
fo_flows = 
DF01 ~                              # Response for flow-type (DF) (binary) 
s(tp_2, k=maxk) +                   # Short term precipitation (smoothed)
s(tp_21, k=maxk) +                  # Antecedent precipitation (21-day) (smoothed)
s(tp_su_mean_annual_mean, k=maxk) + # Mean annual precipitation (smoothed)
s(DF_Slope ,   k=maxk) +            # Slope angle in flow-type PPA (smoothed)
s(DF_CI_mean , k=maxk) +            # Convergence index in flow-type PPA (smoothed)
DF_dominant_landcover +             # Dominant land cover class in flow-type PPA (parametric categorical)
DF_dominant_lithology +             # Dominant lithology class in flow-type PPA (parametric categorical)
s(log_B_DF_n , k = maxk)+           # Log Number of buildings in flow-type PPA (smoothed)
s(doy , bs="cc", k=maxk) +          # Day-of-Year as circular variable (cyclic spline)
s(temp_mean, k=maxk) +              # Mean daily temperature (smoothed)
temp_sub0 +                         # Binary: Day below 0 °C (yes/no) (parametric)
s(cat,  bs="re") +                  # Random effect: Sampling location (Basin-ID) using bs="re"
s(year, bs="re")                    # Random effect: Sampling Year using bs="re"

# Fit the flow type model ---------------------------------------------------------------------------
fit_flows = mgcv::bam(fo_flows, data=df_flows, family=binomial, discrete = 100, method="fREML", select=T) # select=T enables automatic smooth selection
summary(fit_flows) # check model summary
plot(fit_flows, pages = 1) # quick partial effect plots (exploratory; use gratia::draw for publication-quality plots --> see Step: 6)

# Define formula for fall types (RF) -------------------------------------------------------------------
fo_falls = 
RF01 ~                               # Response for fall-type (RF) (binary) 
s(tp_2, k=maxk) +                    # Short term precipitation (smoothed)
s(tp_14 , k=maxk) +                  # Antecedent precipitation (14-day) (smoothed)
s(tp_su_mean_annual_mean, k=maxk) +  # Mean annual precipitation (smoothed)
s(RF_Slope ,   k=3) +                # Slope angle in PPA; k=3 chosen for less wiggliness
s(RF_CI_mean , k=maxk) +             # Convergence index in fall-type PPA (smoothed)
s(RF_bare.other , k=maxk) +          # Portion of bare surface in fall-type PPA (smoothed)
RF_dominant_landcover +              # Dominant land cover class in fall-type PPA (parametric categorical)
RF_dominant_lithology +              # Dominant lithology class in fall-type PPA (parametric categorical)
s(log_B_RF_n , k = maxk)+            # Log Number of buildings in fall-type PPA (smoothed)
s(RF_roads_12 , k=maxk) +            # Portion of “road” cells in fall PPA (smoothed)
temp_mixed +                         # Binary: Day crossing 0 °C (yes/no) (parametric)
s(cat,  bs="re") +                   # Random effect: Sampling location (Basin-ID) using bs="re"
s(year, bs="re")                     # Random effect: Sampling Year using bs="re"

# Fit the fall type model ---------------------------------------------------------------------------
fit_falls = mgcv::bam(fo_falls, data=df_falls, family=binomial, discrete = 100, method="fREML", select=T) # select=T enables automatic smooth selection
summary(fit_falls) # check model summary
plot(fit_falls, pages = 1) # quick partial effect plots (exploratory; use gratia::draw for publication-quality plots --> see Step: 6)
#
##
###
##########
########## Step 3:  Calculate fitting performance and plot ROC curves (Fig. 4a in publication) ----------------------------------------------------------------------------------
##########
###
##
#
myexclude = c("s(cat)", "s(year)") # Argument for averaging out the random effects (bs="re") during prediction (excluding effects)

# Fitting performance slides ----------------------------------------------------------------------------
df_slides$probSL = predict(fit_slides, type = "response", newdata = df_slides, exclude = myexclude) # predict probabilities for training data w/ random effects excluded
myrocs <- roc(response=df_slides$SL01, predictor=df_slides$probSL, auc=T, ci=T); myrocs$auc # roc object and AUC (with CI)
bestSL  = pROC::coords(myrocs, x="best", input="threshold",  best.method="c") # identify "best" threshold (c method = closest to (0,1))

# Fitting performance flows ----------------------------------------------------------------------------
df_flows$probDF = predict(fit_flows, type = "response", newdata = df_flows, exclude = myexclude) # predict for flows
myrocd <- roc(response=df_flows$DF01, predictor=df_flows$probDF, auc=T, ci=T); myrocd$auc # roc object and AUC (with CI)
bestDF  = pROC::coords(myrocd, x="best", input="threshold",  best.method="c") # identify "best" threshold (c method = closest to (0,1))

# Fitting performance falls ----------------------------------------------------------------------------
df_falls$probRF = predict(fit_falls, type = "response", newdata = df_falls, exclude = myexclude) # predict for falls
myrocr <- roc(response=df_falls$RF01, predictor=df_falls$probRF, auc=T, ci=T); myrocr$auc # roc object and AUC (with CI)
bestRF  = pROC::coords(myrocr, x="best", input="threshold",  best.method="c") # identify "best" threshold (c method = closest to (0,1))

# Fitting performance plot ROC ----------------------------------------------------------------------------
# === Labels for AUCs ===
auc_labels <- data.frame(
  Process = c("Slide-type", "Flow-type", "Fall-type"), # label/order
  AUC = c(round(auc(myrocs), 3), round(auc(myrocd), 3), round(auc(myrocr), 3)), # AUC rounded for reporting
  x = 0.4,y = c(0.15, 0.11, 0.07)) # coordinates for annotation
auc_vals <- setNames(paste0(auc_labels$Process, " (AUC = ", round(auc_labels$AUC, 2), ")"), auc_labels$Process) # label string used in legend
auc_vals # print to console for quick check

# === Convert ROC objects to data frames ===
roc_slide <- data.frame(
  FPR = 1 - myrocs$specificities, # false positive rate
  TPR = myrocs$sensitivities,     # true positive rate
  Process = "Slide-type")          # label

roc_flow <- data.frame(
  FPR = 1 - myrocd$specificities,
  TPR = myrocd$sensitivities,
  Process = "Flow-type")

roc_fall <- data.frame(
  FPR = 1 - myrocr$specificities,
  TPR = myrocr$sensitivities,
  Process = "Fall-type")

roc_all <- bind_rows(roc_slide, roc_flow, roc_fall) # combine ROC curves into one data frame for ggplot
roc_all$Process     <- factor(roc_all$Process, levels = c("Slide-type", "Flow-type", "Fall-type")) # set factor levels for plotting order

# === Cutpoints (TPR, FPR) ===
cutpoints <- data.frame(
  FPR = c(1 - bestSL$specificity, 1 - bestDF$specificity, 1 - bestRF$specificity), # FPR at best thresholds
  TPR = c(bestSL$sensitivity, bestDF$sensitivity, bestRF$sensitivity), # TPR at best thresholds
  Process = c("Slide-type", "Flow-type", "Fall-type"))
cutpoints$Process  <- factor(cutpoints$Process, levels = c("Slide-type", "Flow-type", "Fall-type")) # ensure order
cutpoints # print cutpoints
# Define a common theme to unify font sizes, etc. -----------------------------------------------------
common_theme <- theme_minimal(base_size = 12) + # base theme for ggplot
  theme(
    plot.title = element_text(size = 12, face = "bold"),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 12),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 12),
    panel.border = element_rect(color = "black", fill = NA, size = 0.5)) # subtle border for clarity

# Build ROC plot ---------------------------------------------------------------------------------------
myroc <- ggplot(roc_all, aes(x = FPR, y = TPR, color = Process)) +
  geom_line(size = 1) + # ROC curves
  geom_abline(linetype = "dashed", color = "gray60") + # no-skill line
  geom_point(data = cutpoints, aes(x = FPR, y = TPR, color = Process), size = 2) + # best-threshold points
  scale_color_manual(
    values = c("Slide-type" = cols, "Flow-type" = cold, "Fall-type" = colr),
    labels = auc_vals,
    breaks = c("Slide-type", "Flow-type", "Fall-type")) +
  labs(title = "Fitting performance",
       x = "False Positive Rate (FPR)",
       y = "True Positive Rate (TPR)",
       color = NULL) +  # removes legend title (we use labels instead)
  coord_equal() + # equal aspect ratio to avoid distortion
  common_theme +  # apply common theme
theme(legend.position = c(0.95, 0.05),  # place legend inside plot at bottom-right
  legend.justification = c("right", "bottom"),
  legend.background = element_rect(fill = alpha("white", 0.7), color = NA),  # translucent bg for legend
  legend.key.size = unit(1.5, "lines"), legend.text = element_text(size = 12))
# plot roc
myroc # render ROC plot (Fig 4a)
#
##
###
##########
########## Step 4: Cross validation: Snippet code for basin-based cross validation and visualizing final performance (Fig. 4a-d in publication) ----------------------------------------------------------------------------------
##########
###
##
#
# This section describes how basin-based cross-validation was performed and provides data and plots
#
# Description:
# - The function "partition_factor_cv" (package: sperrorest) was used to create train/test splits based on basin IDs ("cat").
# - Five folds and ten repetitions were applied.
# - The process was repeated for each model (slide, flow, fall) and for multiple antecedent precipitation time windows (tp_null, tp_7, tp_14, tp_21, tp_30).
#
# Example (commented) usage:
#>nreps <- 10                                     # Number of repetitions (10x repeated CV)                        
#>nfolds <- 5                                     # Number of folds (5-fold CV)
#> parti <- partition_factor_cv(                  # Function to generate basin-based cross-validation splits
#>   df_slides,                                   # Input data frame (here slide-type data)
#>   nfold = nfolds,                              # Number of folds (5-fold CV)
#>   repetition = nreps,                          # Number of repetitions (10x repeated CV)
#>   seed1 = 666,                                 # Random seed for reproducibility
#>   fac = "cat"                                  # Factor variable ("cat") defining basin IDs 
#> )  # Creates a list of partition objects for repeated 5-fold CV
# Example iteration (commented snippet):
#>results_df_slides <- data.frame(Formula = character(),Repetition = integer(),Fold = integer(),AUROC = numeric(),stringsAsFactors = FALSE) # Create an empty df to store results
#> for (formula_name in names(formula_list)) {    # Iterate over each formula (formula_list contain gam-formulas (cf. Step 2) with varying antec. precip. time windows)
#>   fo <- formula_list[[formula_name]]           # Extract current formula by name
#>   print(fo)                                    # Print formula to console for tracking progress
# Loop through repetitions and folds
#>   for (j in 1:nreps) {                         # Outer loop over number of repetitions
#>     partiloop <- parti[[j]]                    # Access the j-th partition (one repetition)
#>     for (i in 1:nfolds) {                      # Inner loop over folds within repetition
#>       first <- partiloop[[i]][[2]]             # Extract test indices for the i-th fold
#>       test <- df_slides[first, ]; ntesti <- nrow(test)  # Create test set and count its rows
#>       train <- df_slides[-first, ]             # Define training data by excluding test indices
#>       # Fit model and calculate AUROC
#>       myfit <- mgcv::bam(                     # Fit GAM model using mgcv::bam (fast REML method)
#>         fo,                                   # Model formula
#>         data = train,                         # Training data subset
#>         family = binomial,                    # Binary family for presence/absence
#>         method = "fREML",                     # Fitting method: fast REML
#>         discrete = 100)                       # Use discrete approximation (speed optimization)
#>       test$prob <- predict.gam(               # Predict probabilities on the test set
#>         myfit,                                # Fitted GAM model
#>         type = "response",                    # Output predicted probabilities (0–1)
#>         newdata = test,                       # Data for prediction
#>         exclude = myexclude)                  # Exclude random effects (as defined above)
#>       myroc <- roc(                           # Compute ROC using pROC::roc
#>         response = test$SL01,                 # Binary response
#>         predictor = test$prob,                # Predicted probabilities
#>         auc = TRUE)                            # Return AUC value
#>       auroc <- round(myroc$auc, 5)            # Round AUC
#>       print(auroc)                            # Print current AUC to console
#>       results_df_slides <- rbind(             # Append results to a data frame
#>         results_df_slides,                    # Existing results table
#>         data.frame(                           # Add new row with results for current fold
#>           Formula = formula_name,             # Formula identifier
#>           Repetition = j,                     # Current repetition number
#>           Fold = i,                           # Current fold number
#>           AUROC = auroc                       # Calculated AUC value
#>         ))}}}
#> print(results_df_slides)                      # Display full table of AUROC results
#> mySLcv_5f_10r <- results_df_slides            # Save cross-validation results object for slides (cf. data loaded in Step 1)
# Repeat the same procedure for flows (DF) and falls (RF)
# --- End of cross-validation example. Below: summarizing final results for all process types.

# extract stats from cross-validation outputs (slide / flow / fall) ------------------------------------
SLcv %>% group_by(Formula) %>% summarise(Mean_AUROC = mean(AUROC), Median_AUROC = median(AUROC), Min_AUROC = min(AUROC), Max_AUROC = max(AUROC), IQR_AUROC = IQR(AUROC)) # slides summary
DFcv %>% group_by(Formula) %>% summarise(Mean_AUROC = mean(AUROC), Median_AUROC = median(AUROC), Min_AUROC = min(AUROC), Max_AUROC = max(AUROC), IQR_AUROC = IQR(AUROC)) # flows summary
RFcv %>% group_by(Formula) %>% summarise(Mean_AUROC = mean(AUROC), Median_AUROC = median(AUROC), Min_AUROC = min(AUROC), Max_AUROC = max(AUROC), IQR_AUROC = IQR(AUROC)) # falls summary

# custom ordering / labels for antecedent precipitation windows ---------------------------------------
custom_order <- c("fo_null", "fo_7", "fo_14", "fo_21", "fo_30") # Replace with order for plotting 
custom_labels <- c("Null", "7 days", "14 days", "21 days", "30 days") # Labels shown on plots

# Convert Formula to factor with custom labels for consistent x-axis labels ----------------------------
SLcv$Formula <- factor(SLcv$Formula, levels = custom_order, labels = custom_labels)
DFcv$Formula <- factor(DFcv$Formula, levels = custom_order, labels = custom_labels)
RFcv$Formula <- factor(RFcv$Formula, levels = custom_order, labels = custom_labels)

# Compute median AUROC per group ----------------------------------------------------------------------------
SLcv_medians <- SLcv %>% group_by(Formula) %>% summarise(median_AUROC = median(AUROC, na.rm = TRUE))
max_medianSL <- max(SLcv_medians$median_AUROC) # used for title annotation
DFcv_medians <- DFcv %>% group_by(Formula) %>% summarise(median_AUROC = median(AUROC, na.rm = TRUE))
max_medianDF <- max(DFcv_medians$median_AUROC)
RFcv_medians <- RFcv %>% group_by(Formula) %>% summarise(median_AUROC = median(AUROC, na.rm = TRUE))
max_medianRF <- max(RFcv_medians$median_AUROC)

# IQR metrics used in text or figure annotations -------------------------------------------------------
(iqr_SL_30 <- IQR(SLcv$AUROC[SLcv$Formula == "30 days"], na.rm = TRUE)) # Slide-type IQR for 30d formula
(iqr_DF_21 <- IQR(DFcv$AUROC[DFcv$Formula == "21 days"], na.rm = TRUE)) # Flow-type IQR for 21d
(iqr_RF_14 <- IQR(RFcv$AUROC[RFcv$Formula == "14 days"], na.rm = TRUE)) # Fall-type IQR for 14d

# Boxplots for each process type with unified theme ----------------------------------------------------
p1 <- ggplot(SLcv, aes(x = Formula, y = AUROC)) +
  geom_boxplot(fill = cols) + # slide boxplot using predefined color
 labs(title = paste0("Slide-type: Predictive performance\nBest: 30 days (median: ", round(max_medianSL, 2), ")"),
       x = "Antecedent precipitation time window",
       y = "AUROC") +
  coord_cartesian(ylim = c(0.725, 0.925)) + # y-limits tuned for visibility across plots
  common_theme

p2 <- ggplot(DFcv, aes(x = Formula, y = AUROC)) +
  geom_boxplot(fill = cold) + # flow boxplot
 labs(title = paste0("Flow-type: Predictive performance\nBest: 21 days (median: ", round(max_medianDF, 2), ")"),
       x = "Antecedent precipitation time window",
       y = "AUROC") +
  coord_cartesian(ylim = c(0.725, 0.925)) +
  common_theme

p3 <- ggplot(RFcv, aes(x = Formula, y = AUROC)) +
  geom_boxplot(fill = colr) + # fall boxplot
 labs(title = paste0("Fall-type: Predictive performance\nBest: 14 days (median: ", round(max_medianRF, 2), ")"),
       x = "Antecedent precipitation time window",
       y = "AUROC") +
  coord_cartesian(ylim = c(0.725, 0.925)) +
  common_theme

## plot for fitting and predictive performance ---------------------------------------------------------
mymargin <- theme(plot.margin = unit(c(0.8, 0.8, 0.8, 0.8), "cm"))  # top, right, bottom, left margins for consistent spacing
myroc <- myroc + mymargin # apply margin to ROC plot
p1 <- p1 + mymargin       # apply margin to CV plots
p2 <- p2 + mymargin
p3 <- p3 + mymargin
myroc_tagged <- myroc + labs(tag = "a)") # tag letters for figure panels
p1_tagged <- p1 + labs(tag = "b)")
p2_tagged <- p2 + labs(tag = "c)")
p3_tagged <- p3 + labs(tag = "d)")

# Combine plots and keep tags using patchwork ----------------------------------------------------------
combined_plot <- (myroc_tagged + p1_tagged) / (p2_tagged + p3_tagged) +  # 2x2 layout: ROC + p1 top row, p2 + p3 bottom
  plot_layout(tag_level = 'keep') &  # keep tags on all panels
  theme(plot.tag = element_text(size = 16, face = "bold", hjust = 0, vjust = 1)) # style tags
print(combined_plot) # render combined figure (Figure 4)
#
##
###
##########
########## Step 5: Variable importance plots (incl. snippet code for calculating feature importance) (Fig. 5) ----------------------------------------------------------------------------------
##########
###
##
#
# Snippet code for slide-type models (SL) -- permutation importance via vip::vi_permute (for fast demo nsim small)
mynsim = 5 # set the number of permutations per variable (set to 100 in original publication -> slow)

set.seed(1) # reproducible permutations
# Select the relevant columns from fitted model
cn_s <- colnames(fit_slides$model)
train_sel_s <- df_slides %>%
  dplyr::select(any_of(cn_s)) %>% # select columns used in model
  dplyr::select(-SL01)            # drop response from training set used for vi_permute
# Define the target variable
target_s <- as.factor(df_slides$SL01) # binary response

# Calculate variable importance using permutation (vi_permute returns Importance & StDev)
result_s <- vi_permute(
  nsim = mynsim,                      # number of permutations (small here for speed)
  object = fit_slides,                # fitted mgcv::bam object
  train = train_sel_s,                # training predictors
  target = target_s,                  # target vector
  event_level = "second",             # which factor level is event (depending on factor coding)
  metric = "roc_auc",                 # metric to evaluate (AUROC)
  pred_wrapper = predict.bam          # use predict.bam as wrapper for consistency
) %>%
  dplyr::arrange(-Importance) %>%     # sort by decreasing importance
  dplyr::mutate(process = "Slide-type") # annotate process type

# Optionally drop random effect variables from importance table
result_s <- result_s %>%
  dplyr::filter(!(Variable %in% c("cat", "year"))) # remove random effect vars since permuting them is not meaningful
print(result_s) # print slide VI results (use in SI or checks)
# repeat for flow-types and fall-types and create "varimpo" data (increase mynsim beforehand)

# Print raw and processed importance outputs ----------------------------------------------------------------
print(varimpo_raw, n = 33) # raw results 
print(varimpoSL) # sorted and renamed importance results for slide-types (SL)
print(varimpoDF) # sorted and renamed importance results for flow-types (DF)
print(varimpoRF) # sorted and renamed importance results for fall-types (RF)

# create plots for variable importance (ggplot2) -------------------------------------------------------
p1 <- ggplot(varimpoSL, aes(x = Importance, y = Variable)) +
  geom_point(color = cols, size = 3) +  # point for importance
  geom_errorbarh(aes(xmin = Importance - StDev, xmax = Importance + StDev), 
                  height = 0.3, color = cols, size = 1) +  # horizontal error bars represent permute SD
  theme_minimal(base_size = 14) +  # readable baseline font size
  labs(title = "Slide-type", x = "Importance", y = "Variable") +
  theme(axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        plot.title = element_text(size = 15, face = "bold"))

p2 <- ggplot(varimpoDF, aes(x = Importance, y = Variable)) +
  geom_point(color = cold, size = 3) +  
  geom_errorbarh(aes(xmin = Importance - StDev, xmax = Importance + StDev), 
                  height = 0.3, color = cold, size = 1) +  
  theme_minimal(base_size = 14) +  
  labs(title = "Flow-type", x = "Importance", y = "Variable") +
  theme(axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        plot.title = element_text(size = 15, face = "bold"))

p3 <- ggplot(varimpoRF, aes(x = Importance, y = Variable)) +
  geom_point(color = colr, size = 3) +  
  geom_errorbarh(aes(xmin = Importance - StDev, xmax = Importance + StDev), 
                  height = 0.3, color = colr, size = 1) +  
  theme_minimal(base_size = 14) +  
  labs(title = "Fall-type", x = "Importance", y = "Variable") +
  theme(axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        plot.title = element_text(size = 15, face = "bold"))

# --- Tagging the variable importance plots ---
p1_tagged <- p1 + labs(tag = "a)") # tag for panel referencing in manuscript
p2_tagged <- p2 + labs(tag = "b)")
p3_tagged <- p3 + labs(tag = "c)")

# Arrange plots vertically, keeping tags (patchwork) --------------------------------------------------
combined_plot <- (p1_tagged) / (p2_tagged) / (p3_tagged) +
  plot_layout(tag_level = 'keep') & # keep patchwork tags across panels
  theme(plot.tag = element_text(size = 16, face = "bold", hjust = 0, vjust = 1))
print(combined_plot) # render VI plots (Fig. 5)
#
##
###
##########
########## Step 6: Visualizing partial effects from fitted models (Fig. 6, Fig. 7, Fig. 8) ----------------------------------------------------------------------------------
##########
###
##
#
# Partial effects created separately for SL, DF, RF using gratia::draw 
# Partial effect plot for slide-type from fitted model (fit_slides)
mys = smooths(fit_slides) # list of smooth terms available in fit_slides
myselect = c(mys[1], mys[2], mys[3], mys[4], mys[5], mys[6], mys[7],mys[8]) # indices of smooths to extract (custom selection)
myselect # print selection to console for verification

# Define custom labels for each plot
custom_titles <- c(
  "Short-term precipitation (mm)", "Antecedent precipitation (mm) — 30 days", "Mean annual precipitation (mm)", "Slope angle (°)", "Number of buildings (log scale)", "Transport infrastructure (coverage %)",
  "Day-of-Year", "Mean daily temperature (°C)", "Land cover", "Lithology", "Day below 0 °C")

# Extract individual plots instead of combining into patchwork directly (draw() from gratia)
plots <- draw(fit_slides, data = df_slides, select = myselect, parametric = TRUE, ci_level = 0.95, ci_col = cols, smooth_col = "#070000",
              return = "list")  # return list of ggplot objects for fine-grained customization

# Post-process plots: apply consistent styling and handle discrete variables ------------------------------------------------
plots <- lapply(seq_along(plots), function(i) {
  p <- plots[[i]] +
    labs(title = custom_titles[i], x = custom_titles[i]) + # set title & x-label
    theme_bw(base_size = 10) + # classical background for publication
    theme(
      axis.text = element_text(size = 8.5),
      axis.title = element_text(size = 9),
      plot.title = element_blank(), # titles handled separately in arrangement
      legend.text = element_text(size = 10),
      legend.title = element_text(size = 100))
  # Manually relabel the ticks for the 11th plot (discrete variable) if present
  if (i == 11) {
    p <- p + scale_x_discrete(labels = c("No", "Yes")) # nicer labels for binary factor}
  return(p)})
length(plots) # check number of plots produced
# Tag letters for each plot ---------------------------------------------------------------------------------
tags <- letters[1:length(plots)] # generates sequence "a", "b", "c", ...
plots <- lapply(seq_along(plots), function(i) {
  p <- plots[[i]] +
    labs(title = custom_titles[i], x = custom_titles[i], tag = paste0(tags[i], ")")) + # tag inserted for panel lettering
    theme_bw(base_size = 10) +
    theme(
      axis.text = element_text(size = 8.5),
      axis.title = element_text(size = 9),
      plot.title = element_blank(),
      plot.tag = element_text(size = 14, face = "bold", hjust = 0, vjust = 1), # style tag (letter)
      legend.text = element_text(size = 10),
      legend.title = element_text(size = 100))  
  if (i == 11) { # discrete variable adjustment
    p <- p + scale_x_discrete(labels = c("No", "Yes"))}
  return(p)})
# Combine plots with tags preserved (wrap_plots from patchwork)
p_SL <- wrap_plots(plots, ncol = 4) + plot_layout(tag_level = 'keep') # arrange into grid, 4 columns
print(p_SL) # render slide-type partial effects (Fig. 6)

# Partial effect plot for flow-type from fitted model (fit_flows)
mys = smooths(fit_flows); mys # show available smooths for flow model
summary(fit_flows) # model diagnostics reminder
myselect = c(mys[1], mys[2], mys[3], mys[4], mys[5], mys[6], mys[7],mys[8]) # selection indices for flow model smooths
myselect # print selection

# Define custom labels for each plot for flows
custom_titles <- c(
  "Short-term precipitation (mm)", "Antecedent precipitation (mm) — 21 days", "Mean annual precipitation (mm)", "Slope angle (°)", "Convergence index", "Number of buildings (log scale)",
  "Day-of-Year", "Mean daily temperature (°C)", "Land cover", "Lithology", "Day below 0 °C")

# Extract individual plots for flow model using gratia::draw()
plots <- draw(fit_flows, data = df_flows, select = myselect, parametric = TRUE, ci_level = 0.95, ci_col = cold, smooth_col = "#070000",
              return = "list") 

# Post-process flow plots: styling and discrete label handling
plots <- lapply(seq_along(plots), function(i) {
  p <- plots[[i]] +
    labs(title = custom_titles[i], x = custom_titles[i]) +
    theme_bw(base_size = 10) +
    theme(
      axis.text = element_text(size = 8.5),
      axis.title = element_text(size = 9),
      plot.title = element_blank(),
      legend.text = element_text(size = 10),
      legend.title = element_text(size = 100))
  # discrete variable relabeling if necessary
  if (i == 11) {
    p <- p + scale_x_discrete(labels = c("No", "Yes"))}
  return(p)})
length(plots) # number of flow plots produced

# Tag letters for Flow-type plots and styling
tags <- letters[1:length(plots)]
plots <- lapply(seq_along(plots), function(i) {
  p <- plots[[i]] +
    labs(title = custom_titles[i], x = custom_titles[i], tag = paste0(tags[i], ")")) +
    theme_bw(base_size = 10) +
    theme(
      axis.text = element_text(size = 8.5),
      axis.title = element_text(size = 9),
      plot.title = element_blank(),
      plot.tag = element_text(size = 14, face = "bold", hjust = 0, vjust = 1),
      legend.text = element_text(size = 10),
      legend.title = element_text(size = 100) )
  if (i == 11) { # discrete variable adjustment
    p <- p + scale_x_discrete(labels = c("No", "Yes")) }
  return(p)})

# Combine and keep tags for flow plots
p_DF <- wrap_plots(plots, ncol = 4) + plot_layout(tag_level = 'keep') # 4 columns layout
print(p_DF) # render flow partials (Fig. 7)

# Partial effect plot for fall-type from fitted model (fit_falls)
mys = smooths(fit_falls); mys # show available smooths for falls
summary(fit_falls) # model summary check
myselect = c(mys[1], mys[2], mys[3], mys[4], mys[5], mys[6], mys[7],mys[8]) # selection indices for fall model
myselect # print selected smooths

# Define custom labels for fall plots
custom_titles <- c(
  "Short-term precipitation (mm)", "Antecedent precipitation (mm) — 14 days", "Mean annual precipitation (mm)", "Slope angle (°)", "Convergence index", "Bare surface (coverage %)", "Number of buildings (log scale)",
  "Transport infrastructure (coverage %)", "Land cover", "Lithology", "Day crossing 0 °C")

# Extract partial plots for fall-type model
plots <- draw(fit_falls, data = df_falls, select = myselect, parametric = TRUE, ci_level = 0.95, ci_col = colr, smooth_col = "#070000",
              return = "list") 

# Style and discrete handling for fall plots
plots <- lapply(seq_along(plots), function(i) {
  p <- plots[[i]] +
    labs(title = custom_titles[i], x = custom_titles[i]) +
    theme_bw(base_size = 10) +
    theme(
      axis.text = element_text(size = 8.5),
      axis.title = element_text(size = 9),
      plot.title = element_blank(),
      legend.text = element_text(size = 10),
      legend.title = element_text(size = 100))
  # Manually relabel the ticks for the 11th plot (discrete variable)
  if (i == 11) {
    p <- p + scale_x_discrete(labels = c("No", "Yes"))}
  return(p)})

length(plots) # number of fall plots
# Tag letters for Fall-type plots
tags <- letters[1:length(plots)]
plots <- lapply(seq_along(plots), function(i) {
  p <- plots[[i]] +
    labs(title = custom_titles[i], x = custom_titles[i], tag = paste0(tags[i], ")")) +
    theme_bw(base_size = 10) +
    theme(
      axis.text = element_text(size = 8.5),
      axis.title = element_text(size = 9),
      plot.title = element_blank(),
      plot.tag = element_text(size = 14, face = "bold", hjust = 0, vjust = 1),
      legend.text = element_text(size = 10),
      legend.title = element_text(size = 100))
  if (i == 11) { # discrete variable adjustment
    p <- p + scale_x_discrete(labels = c("No", "Yes"))}
  return(p)})

# Combine and keep tags for fall-type partials
p_RF <- wrap_plots(plots, ncol = 4) + plot_layout(tag_level = 'keep') # 4 columns grid
print(p_RF) # render fall partials (Fig. 8)

# End of script ---------------------------------------------------------------------------------------


