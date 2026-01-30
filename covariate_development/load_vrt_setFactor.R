
# loard vrt, assign factor levels

# Apply categorical levels (RAT) from an RDS to the matching bands in a VRT.
# Assumes the RDS is a list of data.frames; each df has 'ID' plus one label column
# named for its band (e.g., 'gmrph_ms_30', 'morpfeat_8', ...).

suppressPackageStartupMessages({
  required.packages <- c("terra")
  new.packages <- required.packages[!(required.packages %in% installed.packages()[, "Package"])]
  if (length(new.packages)) install.packages(new.packages)
  lapply(required.packages, require, character.only = TRUE)
})

## ---- 1) Paths ---------------------------------------------------------------
vrt_path <- "~/spodic/tiles/cov.vrt"
levels_rds_path <- "~/spodic/tiles/categorical_levels.rds"
out_tif <- "~/spodic/tiles/cov_with_levels.tif"  # optional output

## ---- 2) Load VRT ------------------------------------------------------------
rStack <- rast(vrt_path)
message(sprintf("Loaded VRT: %s", vrt_path))
message(sprintf("Bands in VRT: %d", nlyr(rStack)))
vrt_names <- names(rStack)
if (is.null(vrt_names)) {
  warning("VRT bands do not have names. Mapping by name will not work. ",
          "You will need to provide band indices explicitly.")
}
message("Preview VRT band names (first 30): ",
        paste(head(vrt_names, 30), collapse = ", "))

## ---- 3) Load levels RDS -----------------------------------------------------
levels_rds <- readRDS(levels_rds_path)
if (!is.list(levels_rds)) {
  stop("categorical_levels.rds must contain a list of data.frames (levels per categorical band).")
}
message(sprintf("RDS contains %d categorical level tables.", length(levels_rds)))

## ---- 4) Extract band keys from the RDS --------------------------------------
# Each data.frame has 'ID' and one band-specific label column (e.g., 'gmrph_ms_30').
get_band_key <- function(df) {
  # The band key is the non-ID column name; 'value' might also appear as ID in some cases.
  cols <- names(df)
  # exclude common ID column names
  label_candidates <- setdiff(cols, c("ID", "value"))
  if (length(label_candidates) < 1) {
    stop("Could not find a label column (non-ID) in one of the RDS data.frames.")
  }
  label_candidates[1]
}
band_keys <- vapply(levels_rds, get_band_key, FUN.VALUE = character(1))
message("Band keys found in RDS: ", paste(band_keys, collapse = ", "))

## ---- 5) Map band keys to VRT band indices ----------------------------------
# Preferred: exact name match in VRT names
map_keys_to_indices <- function(keys, vrt_names) {
  if (is.null(vrt_names)) {
    # no names present; cannot map by name
    return(rep(NA_integer_, length(keys)))
  }
  idx <- match(keys, vrt_names)
  # Fallback: if exact match fails, try regex (contains) to locate a unique match
  for (i in seq_along(keys)) {
    if (is.na(idx[i])) {
      hits <- grep(keys[i], vrt_names, fixed = TRUE)
      if (length(hits) == 1) {
        idx[i] <- hits
      } else if (length(hits) > 1) {
        warning(sprintf("Multiple VRT bands match key '%s' by substring: indices %s. ",
                        keys[i], paste(hits, collapse = ", ")))
      } else {
        warning(sprintf("No VRT band name matches key '%s'.", keys[i]))
      }
    }
  }
  idx
}
target_idx <- map_keys_to_indices(band_keys, vrt_names)
message("Mapped indices in VRT for keys: ",
        paste(sprintf("%s->%s", band_keys,
                      ifelse(is.na(target_idx), "NA", as.character(target_idx))),
              collapse = "; "))

if (any(is.na(target_idx))) {
  stop("Some keys could not be mapped to VRT bands. ",
       "Please ensure VRT band names contain: ",
       paste(band_keys[is.na(target_idx)], collapse = ", "),
       "\nIf VRT has no names, provide explicit band indices for these keys and re-run.")
}

## ---- 6) Normalize and apply levels to matched bands -------------------------
normalize_levels_df <- function(df, expect_label_name) {
  if (is.null(df)) return(NULL)
  # Rename ID/value -> 'ID'
  if (!("ID" %in% names(df))) {
    if ("value" %in% names(df)) {
      names(df)[names(df) == "value"] <- "ID"
    } else {
      stop("Levels table lacks 'ID' or 'value' column.")
    }
  }
  # Rename the specific label column to 'label'
  if (!(expect_label_name %in% names(df))) {
    stop(sprintf("Expected label column '%s' not present in levels table.", expect_label_name))
  }
  names(df)[names(df) == expect_label_name] <- "label"
  
  # Keep only two columns in order and coerce types
  out <- df[, c("ID", "label"), drop = FALSE]
  out$ID <- as.integer(out$ID)
  out$label <- as.character(out$label)
  # Sort for reproducibility
  out[order(out$ID), ]
}

# Prepare levels list for assigning to VRT
l_vrt <- levels(rStack)
if (length(l_vrt) != nlyr(rStack)) {
  l_vrt <- vector("list", nlyr(rStack))
}

for (i in seq_along(levels_rds)) {
  bname <- band_keys[i]
  bidx  <- target_idx[i]
  lev_df <- normalize_levels_df(levels_rds[[i]], expect_label_name = bname)
  l_vrt[[bidx]] <- lev_df
}

# Apply in one shot
levels(rStack) <- l_vrt

## ---- 7) Verify only the affected bands --------------------------------------
flags <- is.factor(rStack)
message("Factor flags for applied bands: ",
        paste(sprintf("%s(band %d)=%s", band_keys, target_idx, flags[target_idx]),
              collapse = "; "))


message("Done.")

rm(list = setdiff(ls(), "rStack"))

