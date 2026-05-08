# =============================================================================
# Mark Twain EcoClass Modeling — Ranger + Terra (weights-only)
# David White NRCS, NM
# Date: 2026-05-08
# =============================================================================

# ------------------------------ CONFIG ---------------------------------------
# Paths
PATH_COV_VRT   <- "~/cov/cov.vrt"
PATH_TRAIN     <- "~/0508/train.RDS"
DIR_TILES      <- "~/cov/tiles"
DIR_RESULTS    <- "~/0508/results/ranger"
DIR_TEMP       <- "~/0508/temp"
covKeep        <- readRDS("~/0508/kept_covariates.rds")

# Modeling/tuning
SEED_GLOBAL    <- 20260423
MIN_CLASS_N    <- 10            # minimum observations per class retained
RFE_STEP       <- 5            # feature subset step size (1..p by this step)
CV_NUMBER      <- 10           # folds
CV_REPEATS     <- 3            # repeats
USE_CLASS_WEIGHTS <- TRUE      # per-fold inverse-frequency weights (capped)

# Tuning grid breadth (decent size search)
MTRY_FRACS     <- seq(0.05, 0.60, by = 0.05)    # as fractions of p
SPLITRULES     <- c("gini", "extratrees")
MIN_NODE_SIZES <- c(2, 4, 8)
NUM_TREES_SET  <- c(200, 400, 600)

# Prediction/outputs
WRITE_PROB_RASTERS <- TRUE     # write per-class probability rasters (multi-band)
WRITE_CLASS_RASTER <- TRUE     # write hard-class categorical raster
GDAL_COMPRESS      <- "LZW"    # compression for outputs
# INT1U supports 0..255 class codes; switch to INT2U if >255 classes
RASTER_DATATYPE    <- "INT1U"

# Report
REPORT_MD          <- file.path(DIR_RESULTS, "report.md")

# -----------------------------------------------------------------------------
# Packages
required.packages <- c(
  "sf","terra","caret","ranger","randomForest","doParallel","parallel","foreach",
  "aqp","foreign","MLmetrics"
)
new.packages <- required.packages[!(required.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
invisible(lapply(required.packages, require, character.only = TRUE))

# Reproducibility & IO
set.seed(SEED_GLOBAL)
dir.create(DIR_RESULTS, showWarnings = FALSE, recursive = TRUE)
dir.create(DIR_TEMP,    showWarnings = FALSE, recursive = TRUE)
terra::terraOptions(memfrac = 0.8, tempdir = DIR_TEMP, progress = 1)

# ----------------------------- DATA PREP -------------------------------------

# Covariates (VRT stack)
rStack <- terra::rast(PATH_COV_VRT)

# subset raster stack
rStack <- rStack[[covKeep]]

# Training points (primary + misc)
pts  <- readRDS(PATH_TRAIN)

# Extract raster covariates at training points
pts.sv <- terra::extract(rStack, pts, xy = FALSE, bind = TRUE, na.rm = TRUE)
pts.sv <- terra::unique(pts.sv)

# Convert to data.frame
all.pts <- sf::st_drop_geometry(sf::st_as_sf(pts.sv))

# Standardize: use abbreviations as modeling class labels
all.pts$class <- as.factor(all.pts$class)

# Drop observations with any missing covariates
all.pts <- all.pts[complete.cases(all.pts), ]
# Drop columns with NAs
all.pts <- all.pts[, colSums(is.na(all.pts)) == 0]

# Filter classes with sufficient support (MIN_CLASS_N)
tbl_class <- table(all.pts$class)
keep      <- names(which(tbl_class >= MIN_CLASS_N))
comp      <- droplevels(all.pts[all.pts$class %in% keep, ])

write.csv(as.data.frame(summary(comp$class)),
          file = file.path(DIR_RESULTS, "ModeledClassBalance.csv"),
          row.names = FALSE)

# Align raster stack to retained covariates
names(comp)
nm <- names(comp)[-c(1)]
rStack <- subset(rStack, nm)

# ------------------------- RECURSIVE FEATURE ELIM ----------------------------

subsets <- seq(1, ncol(comp) - 1, by = RFE_STEP)
number  <- CV_NUMBER
repeats <- CV_REPEATS

# Seeds for RFE
set.seed(12)
num_resamples <- number * repeats
num_tune      <- length(subsets)
seeds_per     <- num_tune + 1
seeds         <- vector("list", length = num_resamples + 1)
for (i in 1:num_resamples) seeds[[i]] <- sample.int(1000, seeds_per)
seeds[[num_resamples + 1]] <- sample.int(1000, 1)

# Parallel for RFE
cl <- parallel::makeCluster(parallel::detectCores() - 2)
doParallel::registerDoParallel(cl)

ctrl.RFE <- caret::rfeControl(
  functions = caret::rfFuncs,
  method    = "repeatedcv",
  number    = number,
  repeats   = repeats,
  seeds     = seeds,
  verbose   = FALSE,
  allowParallel = TRUE
)

set.seed(9)
rfe_results <- caret::rfe(
  x = comp[, -c(1,2)],
  y = comp$class,
  sizes = subsets,
  rfeControl = ctrl.RFE,
  allowParallel = TRUE
)

stopCluster(cl) 
registerDoSEQ() 
gc()

# Selected predictors
a <- caret::predictors(rfe_results)
saveRDS(a, file = file.path(DIR_RESULTS, "preds.rds"))

# Subset comp to selected features
comp.sub <- comp[, c("class", a)]
comp.sub$class <- factor(comp.sub$class)              # sanity
stopifnot(nlevels(comp.sub$class) > 1)

# Predictor hygiene
X <- comp.sub[, -1]
y <- comp.sub$class

# Character -> factor
is_char <- vapply(X, is.character, logical(1))
if (any(is_char)) X[is_char] <- lapply(X[is_char], factor)


# ----------------- Custom caret model: ranger per-fold weights ---------------
# NOTE: No SMOTE or sampling; weights-only balancing per fold.

ranger_perfold_weights <- list(
  label = "Ranger (per-fold class weights)",
  library = "ranger",
  type   = c("Classification"),
  parameters = data.frame(
    parameter = c("mtry", "splitrule", "min.node.size", "num.trees"),
    class     = c("numeric", "character", "numeric", "numeric"),
    label     = c("mtry", "splitrule", "min.node.size", "num.trees")
  ),
  
  grid = function(x, y, len = NULL, search = "grid") {
    p <- ncol(x)
    # expansive grid: mtry as fractions of p, combinatorial across rules, node sizes, trees
    mtrys <- unique(pmax(1, round(MTRY_FRACS * p)))
    expand.grid(
      mtry          = mtrys,
      splitrule     = SPLITRULES,
      min.node.size = MIN_NODE_SIZES,
      num.trees     = NUM_TREES_SET
    )
  },
  
  fit = function(x, y, wts, param, lev, last, classProbs, ...) {
    y <- droplevels(y)
    
    # Per-fold inverse-frequency weights (capped)
    cw <- NULL
    if (isTRUE(USE_CLASS_WEIGHTS)) {
      freq <- table(y)
      cw   <- as.numeric(median(freq) / freq)
      names(cw) <- names(freq)
      cw[cw > 10] <- 10
    }
    
    dat <- data.frame(.y = y, x)
    
    rf <- ranger::ranger(
      dependent.variable.name = ".y",
      data         = dat,
      mtry         = param$mtry,
      splitrule    = as.character(param$splitrule),
      min.node.size= param$min.node.size,
      num.trees    = if (!is.null(param$num.trees)) param$num.trees else 500,
      class.weights= cw,
      importance   = "impurity",
      probability  = isTRUE(classProbs),
      num.threads  = max(1, parallel::detectCores() - 2),
      ...
    )
    rf$obsLevels <- levels(y)
    rf
  },
  
  predict = function(modelFit, newdata, submodels = NULL) {
    pr <- predict(modelFit, data = newdata)
    p  <- pr$predictions
    
    # If probability = TRUE, ranger returns a matrix
    if (is.matrix(p)) {
      # Get the index of the highest probability
      idx <- max.col(p, ties.method = "first")
      # Map index back to class names
      cls <- colnames(p)[idx]
      return(factor(cls, levels = modelFit$obsLevels))
    }
    
    # Otherwise return the vector/factor directly
    return(p)
  },
  
  prob = function(modelFit, newdata, submodels = NULL) {
    pr <- predict(modelFit, data = newdata)
    p  <- pr$predictions
    
    # Ensure it's a data frame so caret can index it by column name
    if (is.matrix(p)) {
      return(as.data.frame(p))
    } else {
      # If probability wasn't TRUE, this will fail multiClassSummary
      return(NULL)
    }
  },
  
  
  # Ensure prob function returns a data frame with proper class names
  prob = function(modelFit, newdata, submodels = NULL) {
    pr <- predict(modelFit, data = newdata)
    as.data.frame(pr$predictions)
  },
  
  levels = function(x) x$obsLevels,
  tags = c("Random Forest", "ranger", "Multiclass", "Class Weights"),
  sort = function(x) x[order(x$Accuracy, decreasing = TRUE), ]
)

# ---------------------------- TRAINING (CV) ----------------------------------

# Parallel for training
cl <- parallel::makeCluster(parallel::detectCores() - 2)
doParallel::registerDoParallel(cl)
parallel::clusterExport(cl, c("MTRY_FRACS", "SPLITRULES", "MIN_NODE_SIZES", 
                              "NUM_TREES_SET", "USE_CLASS_WEIGHTS"))

fitControl <- caret::trainControl(
  method          = "repeatedcv",
  number          = CV_NUMBER,
  repeats         = CV_REPEATS,
  classProbs      = WRITE_PROB_RASTERS,          # if TRUE, caret asks for probs
  savePredictions = "final",
  verboseIter     = TRUE,
  allowParallel   = TRUE,
  summaryFunction = if (WRITE_PROB_RASTERS) caret::multiClassSummary else caret::defaultSummary
)

set.seed(SEED_GLOBAL)
rfm <- caret::train(
  x         = X,
  y         = y,
  method    = ranger_perfold_weights,
  trControl = fitControl,
  metric    = "Accuracy"
)

stopCluster(cl)
registerDoSEQ()
gc()

# Inspect tuning
print(rfm)
print(rfm$bestTune)
top <- rfm$results[order(rfm$results$Accuracy, decreasing = TRUE), ][1:10, ]
write.csv(top, file = file.path(DIR_RESULTS, "top_tuning_rows.csv"), row.names = FALSE)

# Variable importance (caret-level)
vim <- caret::varImp(rfm)$importance
write.csv(vim, file = file.path(DIR_RESULTS, "varImp_caret.csv"))

# ----------------------- FINAL MODEL (RANGER) --------------------------------
# Fit a single final model with bestTune; probability TRUE if writing probs.

levs <- levels(comp.sub$class)
freq <- table(comp.sub$class)
cw   <- NULL
if (isTRUE(USE_CLASS_WEIGHTS)) {
  cw <- as.numeric(median(freq) / freq)
  names(cw) <- names(freq)
  cw[cw > 10] <- 10
}

rf_final_prob <- ranger::ranger(
  dependent.variable.name = "class",
  data         = comp.sub,
  mtry         = rfm$bestTune$mtry,
  splitrule    = rfm$bestTune$splitrule,
  min.node.size= rfm$bestTune$min.node.size,
  num.trees    = rfm$bestTune$num.trees,
  importance   = "impurity",
  class.weights= cw,
  probability  = WRITE_PROB_RASTERS,
  oob.error    = TRUE,
  num.threads  = max(1, parallel::detectCores() - 2)
)

# A hard-label final model (for OOB confusion matrices)
rf_final_hard <- ranger::ranger(
  dependent.variable.name = "class",
  data         = comp.sub,
  mtry         = rfm$bestTune$mtry,
  splitrule    = rfm$bestTune$splitrule,
  min.node.size= rfm$bestTune$min.node.size,
  num.trees    = rfm$bestTune$num.trees,
  importance   = "impurity",
  class.weights= cw,
  probability  = FALSE,
  oob.error    = TRUE,
  num.threads  = max(1, parallel::detectCores() - 2)
)

# Confusion matrices
# CV-based confusion from caret resamples
cm_cv <- caret::confusionMatrix(rfm$pred$pred, rfm$pred$obs)
write.csv(as.table(cm_cv),                    file = file.path(DIR_RESULTS, "ConMat_CV.csv"))
write.csv(as.matrix(cm_cv, what = "overall"), file = file.path(DIR_RESULTS, "Overall_CV.csv"))
write.csv(as.matrix(cm_cv, what = "classes"), file = file.path(DIR_RESULTS, "ClassAccuracy_CV.csv"))

# OOB confusion from hard-label final model
cm_oob <- caret::confusionMatrix(rf_final_hard$predictions, comp.sub$class)
write.csv(as.table(cm_oob),                    file = file.path(DIR_RESULTS, "ConMat_OOB.csv"))
write.csv(as.matrix(cm_oob, what = "overall"), file = file.path(DIR_RESULTS, "Overall_OOB.csv"))
write.csv(as.matrix(cm_oob, what = "classes"), file = file.path(DIR_RESULTS, "ClassAccuracy_OOB.csv"))

# Top confused pairs (from CV confusion)
tab <- as.matrix(cm_cv$table)  # rows = Pred, cols = Ref
col_sums <- colSums(tab); col_sums[col_sums == 0] <- 1
P_ref <- sweep(tab, 2, col_sums, "/")          # P(pred | true)
classes <- colnames(P_ref)
pairs   <- subset(expand.grid(A = classes, B = classes, stringsAsFactors = FALSE), A != B)
pairs$A_to_B <- mapply(function(a,b) P_ref[b, a], pairs$A, pairs$B)
pairs$B_to_A <- mapply(function(a,b) P_ref[a, b], pairs$A, pairs$B)
pairs$Symmetric <- (pairs$A_to_B + pairs$B_to_A)/2
top_pairs <- pairs[order(-pairs$Symmetric), ]
write.csv(head(top_pairs, 25), file = file.path(DIR_RESULTS, "TopConfusedPairs.csv"), row.names = FALSE)

# ------------------------- TILE PREDICTIONS ----------------------------------
# Write tiles to disk first, then mosaic

# List tiles
tile_files <- list.files(DIR_TILES, pattern = "\\.tif$", full.names = TRUE)

# Create output tile directories
DIR_TILE_OUT_CLASS <- file.path(DIR_RESULTS, "tiles_class")
DIR_TILE_OUT_PROBS <- file.path(DIR_RESULTS, "tiles_probs")
dir.create(DIR_TILE_OUT_CLASS, showWarnings = FALSE, recursive = TRUE)
if (WRITE_PROB_RASTERS) dir.create(DIR_TILE_OUT_PROBS, showWarnings = FALSE, recursive = TRUE)

# Parallel tile prediction
cl <- parallel::makeCluster(parallel::detectCores() - 2)
doParallel::registerDoParallel(cl)

# ---- Hard-class per-tile prediction -> files ----
class_tile_paths <- NULL
if (WRITE_CLASS_RASTER) {
  class_tile_paths <- foreach(i = seq_along(tile_files), .packages = c("terra","ranger")) %dopar% {
    tfile <- tile_files[i]
    rti   <- terra::rast(tfile)
    rsel  <- terra::subset(rti, a)               # select RFE predictors
    
    # FIX: Changed type = "class" to type = "response"
    rpred <- terra::predict(rsel, rf_final_hard, type = "response", na.rm = TRUE, steps = 4)
    
    out_i <- file.path(DIR_TILE_OUT_CLASS, paste0("class_", tools::file_path_sans_ext(basename(tfile)), ".tif"))
    terra::writeRaster(
      rpred, filename = out_i, overwrite = TRUE,
      gdal = c(paste0("COMPRESS=", GDAL_COMPRESS)),
      datatype = RASTER_DATATYPE
    )
    out_i
  }
}


# ---- Probability per-tile prediction -> files ----
prob_tile_paths <- NULL
if (WRITE_PROB_RASTERS) {
  prob_tile_paths <- foreach(i = seq_along(tile_files), .packages = c("terra","ranger")) %dopar% {
    tfile <- tile_files[i]
    rti   <- terra::rast(tfile)
    rsel  <- terra::subset(rti, a)
    
    rpr <- terra::predict(rsel, rf_final_prob, fun = function(mod, data) predict(mod, data)$predictions, na.rm = TRUE)
    names(rpr) <- paste0("prob_", levs)          # name layers per class
    out_i <- file.path(DIR_TILE_OUT_PROBS, paste0("prob_", tools::file_path_sans_ext(basename(tfile)), ".tif"))
    terra::writeRaster(
      rpr, filename = out_i, overwrite = TRUE,
      gdal = c(paste0("COMPRESS=", GDAL_COMPRESS)),
      datatype = "FLT4S"
    )
    out_i
  }
}

stopCluster(cl) 
registerDoSEQ()
gc()
# ---- Mosaic tiles from disk (class & probs) ----

# Class mosaic
out_class <- file.path(DIR_RESULTS, "ecoclass.tif")
if (WRITE_CLASS_RASTER) {
  # Read tile rasters and merge
  class_rasts <- lapply(class_tile_paths, terra::rast)
  pred_class  <- do.call(terra::merge, class_rasts)
  
  # Map predicted factor codes to 1..K if needed & build RAT with original long names
  # Frequencies for VAT
  fq <- terra::freq(pred_class, digits = 0)
  # Factor levels: abbreviations (levs) -> original long names
  RAT <- data.frame(
    VALUE = seq_along(levs),
    ABBR  = levs,
    stringsAsFactors = FALSE
  )
  RAT <- merge(RAT, abr_to_long, by.x = "ABBR", by.y = "ABBR", all.x = TRUE)
  names(RAT)[names(RAT) == "CLASS_LONG"] <- "CLASS"
  
  # Ensure counts joined
  if (!is.null(fq) && nrow(fq) > 0) {
    RAT <- merge(RAT, fq[, c("value","count")], by.x = "VALUE", by.y = "value", all.x = TRUE)
    names(RAT)[names(RAT) == "count"] <- "COUNT"
    RAT$COUNT[is.na(RAT$COUNT)] <- 0
  } else {
    RAT$COUNT <- 0L
  }
  
  # Attach RAT
  levels(pred_class) <- list(RAT)
  
  # Datatype guard: if >255 classes, use INT2U
  dt_out <- if (length(levs) <= 255) RASTER_DATATYPE else "INT2U"
  
  terra::writeRaster(
    pred_class, filename = out_class, overwrite = TRUE,
    gdal = c("TFW=YES", paste0("COMPRESS=", GDAL_COMPRESS)),
    datatype = dt_out
  )
  # ESRI VAT .dbf
  foreign::write.dbf(RAT, file = file.path(DIR_RESULTS, "ecoclass.tif.vat.dbf"))
}

# Probability mosaic
out_probs <- file.path(DIR_RESULTS, "ecoclass_probs.tif")
if (WRITE_PROB_RASTERS) {
  prob_rasts <- lapply(prob_tile_paths, terra::rast)
  probs      <- do.call(terra::merge, prob_rasts)
  names(probs) <- paste0("prob_", levs)
  terra::writeRaster(
    probs, filename = out_probs, overwrite = TRUE,
    gdal = c("TFW=YES", paste0("COMPRESS=", GDAL_COMPRESS)),
    datatype = "FLT4S"
  )
}

# --------------------------- SAVE SESSION ------------------------------------
save.image(file = file.path(DIR_RESULTS, "ecoclass-042226.RData"))
writeLines(paste("seed:", SEED_GLOBAL), file.path(DIR_RESULTS, "seed.txt"))

# ----------------------------- REPORT ----------------------------------------
# Build a Markdown report summarizing the run
vim_df <- as.data.frame(vim)

# If 'Overall' doesn't exist, compute the row maximum across all importance columns
if (!"Overall" %in% colnames(vim_df)) {
  # Drop any non-numeric columns first if they exist
  num_cols <- sapply(vim_df, is.numeric)
  
  if (sum(num_cols) > 0) {
    # Take the maximum importance across all classes for each feature
    vim_df$Overall <- apply(vim_df[, num_cols, drop = FALSE], 1, max)
  } else {
    # Fallback: if no numeric columns, just make it the first column
    vim_df$Overall <- vim_df[[1]]
  }
}


md_lines <- c(
  "# EcoClass Modeling Report (Mark Twain)",
  "",
  sprintf("**Run Date:** %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  sprintf("**Seed:** %d", SEED_GLOBAL),
  "",
  "## Paths",
  paste0("- Covariates (VRT): `", PATH_COV_VRT, "`"),
  paste0("- Training SHP: `", PATH_TRAIN, "`"),
  #paste0("- Misc SHP: `", PATH_MISC_SHP, "`"),
  paste0("- Tiles dir: `", DIR_TILES, "`"),
  paste0("- Results dir: `", DIR_RESULTS, "`"),
  "",
  "## Class Balance (after filtering)",
  "",
  "```\n",
  paste(capture.output(print(summary(comp.sub$class))), collapse = "\n"),
  "\n```",
  "",
  "## Feature Selection (RFE)",
  sprintf("- Candidate subset sizes: %s", paste(subsets, collapse = ", ")),
  sprintf("- Selected predictors count: %d", length(a)),
  "",
  "```\n",
  paste(capture.output(print(a)), collapse = "\n"),
  "\n```",
  "",
  "## Tuning Grid & Best Model",
  sprintf("- Grid size: mtry(%d) × splitrule(%d) × min.node.size(%d) × num.trees(%d) = **%d** combos",
          length(unique(pmax(1, round(MTRY_FRACS * ncol(X))))),
          length(SPLITRULES), length(MIN_NODE_SIZES), length(NUM_TREES_SET),
          length(unique(pmax(1, round(MTRY_FRACS * ncol(X))))) *
            length(SPLITRULES) * length(MIN_NODE_SIZES) * length(NUM_TREES_SET)),
  sprintf("- BestTune: mtry=%d, splitrule=%s, min.node.size=%d, num.trees=%d",
          rfm$bestTune$mtry, rfm$bestTune$splitrule, rfm$bestTune$min.node.size, rfm$bestTune$num.trees),
  "",
  "### Top 10 tuning rows by Accuracy",
  "```\n",
  paste(capture.output(print(head(rfm$results[order(rfm$results$Accuracy, decreasing = TRUE), ], 10))), collapse = "\n"),
  "\n```",
  "",
  "## Performance",
  "### CV Metrics (caret resamples)",
  "```\n",
  paste(capture.output(print(cm_cv$overall)), collapse = "\n"),
  "\n```",
  "### OOB Metrics (final hard-label model)",
  "```\n",
  paste(capture.output(print(cm_oob$overall)), collapse = "\n"),
  "\n```",
  "",
  "## Confusion Diagnostics (CV)",
  "### Top 15 mutually confused class pairs",
  "```\n",
  paste(capture.output(print(head(top_pairs[order(top_pairs$Symmetric, decreasing = TRUE), ], 15))), collapse = "\n"),
  "\n```",
  "",
  "## Variable Importance (caret)",
  "```\n",
  paste(capture.output(print(head(vim_df[order(vim_df$Overall, decreasing = TRUE), , drop = FALSE], 20))), collapse = "\n"),
  "\n```",
  "",
  "## Raster Outputs",
  sprintf("- Class mosaic: `%s`", if (WRITE_CLASS_RASTER) out_class else "(not generated)"),
  sprintf("- Probability mosaic: `%s`", if (WRITE_PROB_RASTERS) out_probs else "(not generated)"),
  "",
  "### RAT (Classification) — Original Long Names",
  "Columns: VALUE (code), ABBR (abbreviation), CLASS (long name), COUNT.",
  "Saved as `ecoclass.tif.vat.dbf` and attached to the GeoTIFF.",
  "",
  "## Tile Processing",
  sprintf("- Tiles processed: %d", length(tile_files)),
  sprintf("- Class tiles dir: `%s`", if (WRITE_CLASS_RASTER) DIR_TILE_OUT_CLASS else "(not generated)"),
  sprintf("- Prob tiles dir: `%s`", if (WRITE_PROB_RASTERS) DIR_TILE_OUT_PROBS else "(not generated)"),
  "",
  "## Notes",
  "- Balancing uses **per-fold inverse-frequency class weights only** (no SMOTE).",
  "- Tuning grid breadth covers mtry, splitrule, min.node.size, and num.trees.",
  "- Predictions were **written per tile to disk first** and then mosaicked to reduce peak RAM.",
  "- Classification raster RAT uses **original long class names** (and ABBR).",
  ""
)

# Define your output path (replace with your desired location and filename)
output_html_path <- file.path(DIR_RESULTS, "EcoClass_Modeling_Report.html")

# 1. Write the character vector to a temporary markdown file
temp_md <- tempfile(fileext = ".md")
writeLines(md_lines, temp_md)

# 2. Render the markdown file to a standalone HTML file
# We use rmarkdown::render which automatically applies clean default styling
library(rmarkdown)
rmarkdown::render(
  input = temp_md,
  output_format = html_document(theme = "cosmo", highlight = "textmate"),
  output_file = output_html_path,
  quiet = TRUE
)

# Clean up the temporary markdown file
unlink(temp_md)

cat("HTML report successfully generated at:", output_html_path, "\n")


#writeLines(md_lines, REPORT_MD)

gc()
# =============================================================================
