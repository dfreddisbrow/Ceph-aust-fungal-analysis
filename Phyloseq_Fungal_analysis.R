################################################################################
# PHYLOSEQ WORKFLOW FOR FUNGAL ITS COMMUNITY ANALYSIS 
################################################################################
# Google Gemini Flash was used to help develop and troubleshoot some of the R code: Google. (2025). Gemini 2.5 Flash [Large language model]. https://gemini.google.com
# Set Working Directory
setwd("C:/Users/withe/Desktop/C_aust_fungal_ITS")

# Ensure Rtools compiler PATH is set for Windows build tools
if (Sys.which("make") == "") {
  Sys.setenv(PATH = paste("C:\\rtools44\\usr\\bin", Sys.getenv("PATH"), sep = ";"))
}

# Clear any stale installation file locks
unlink("C:/Users/withe/AppData/Local/R/win-library/4.6/00LOCK*", recursive = TRUE)

# --- 0. PACKAGE INSTALLATION, ENVIRONMENT & OUTPUT DIRECTORY SETUP ---
options(pkgType = "both", INSTALL_opts = "--no-lock")

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

# Install Bioconductor core packages
BiocManager::install(c("phyloseq", "decontam", "edgeR"), update = FALSE)

# Install CRAN dependencies
cran_pkgs <- c("dplyr", "tidyr", "stringr", "ggplot2", "ggsci", "plotly", 
               "scatterpie", "sf", "rnaturalearth", "rnaturalearthdata", 
               "scales", "geosphere", "ggrepel", "colorspace", "remotes", 
               "rlang", "Polychrome", "tibble")
new_pkgs <- cran_pkgs[!(cran_pkgs %in% installed.packages()[,"Package"])]
if(length(new_pkgs)) install.packages(new_pkgs)

if (!requireNamespace("htmlwidgets", quietly = TRUE)) {
  install.packages("htmlwidgets")
}
library(htmlwidgets)

# Install FUNGuildR without triggering locked dependency rebuilds
library(remotes)
if (!requireNamespace("FUNGuildR", quietly = TRUE)) {
  remotes::install_github("brendanf/FUNGuildR", upgrade = "never", INSTALL_opts = "--no-lock")
}

# Load core libraries
library(phyloseq)
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(decontam)
library(FUNGuildR)
library(edgeR)
library(ggsci)
library(plotly)
library(scales)
library(scatterpie)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(geosphere)
library(ggrepel)
library(colorspace)
library(rlang)
library(Polychrome)
library(tibble)

# Setup Output Directories
main_out_dir     <- "output_results"
dir_amplicon_bar <- file.path(main_out_dir, "amplicon_barplots")
dir_amplicon_map <- file.path(main_out_dir, "amplicon_maps")
dir_peloton_bar  <- file.path(main_out_dir, "peloton_barplots")
dir_peloton_map  <- file.path(main_out_dir, "peloton_maps")

invisible(lapply(
  c(dir_amplicon_bar, dir_amplicon_map, dir_peloton_bar, dir_peloton_map),
  dir.create, recursive = TRUE, showWarnings = FALSE
))

# --- HELPER FUNCTIONS ---
canon <- function(x) {
  if (is.null(x)) return(character(0))
  x <- as.character(x)
  x <- trimws(x)
  x[x %in% c("", "NA", "Unassigned", "unassigned", "N/A")] <- NA
  return(x)
}

save_barplot_all_formats <- function(ggplot_obj, plotly_obj, base_filename, out_dir) {
  filepath_base <- file.path(out_dir, base_filename)
  saveWidget(plotly_obj, file = paste0(filepath_base, ".html"), selfcontained = TRUE)
  
  ggsave(paste0(filepath_base, ".png"), plot = ggplot_obj, width = 14, height = 10, dpi = 300)
  ggsave(paste0(filepath_base, ".pdf"), plot = ggplot_obj, width = 14, height = 10, device = "pdf")
  ggsave(paste0(filepath_base, ".svg"), plot = ggplot_obj, width = 14, height = 10)
}

save_map_all_formats <- function(map_plot, base_filename, out_dir) {
  filepath_base <- file.path(out_dir, base_filename)
  ggsave(paste0(filepath_base, ".png"),  plot = map_plot, width = 10, height = 8, dpi = 300)
  ggsave(paste0(filepath_base, ".tiff"), plot = map_plot, width = 10, height = 8, dpi = 600, compression = "lzw")
  ggsave(paste0(filepath_base, ".pdf"),  plot = map_plot, width = 10, height = 8, device = "pdf")
  ggsave(paste0(filepath_base, ".svg"),  plot = map_plot, width = 10, height = 8)
}

# Helper to format full MycoType labels for interactive tooltips
format_mycotype_label <- function(x) {
  case_when(
    x == "OMF" ~ "Orchid Mycorrhizal Fungi (OMF)",
    x == "ECM" ~ "Ectomycorrhizal (ECM)",
    x == "ERM" ~ "Ericoid Mycorrhizal (ERM)",
    x == "AM"  ~ "Arbuscular Mycorrhizal (AM)",
    x == "Saprotroph" ~ "Saprotroph",
    TRUE ~ as.character(x)
  )
}

# ==============================================================================
# --- COLOR PALETTES FOR PUBLICATION FIGURES ---
# ==============================================================================

myco_colors <- c(
  "OMF"        = "#2B5C8F",  # Muted Deep Blue
  "ECM"        = "#2D7044",  # Forest Green
  "ERM"        = "#803975",  # Deep Plum/Orchid
  "AM"         = "#D97706",  # Warm Amber/Gold
  "Saprotroph" = "#8B5CF6"   # Muted Purple
)

get_family_colors <- function(family_vec) {
  unique_fams <- sort(unique(family_vec[!is.na(family_vec) & !family_vec %in% c("Unassigned Family", "Other", "", "EXCLUDE_YeastMoldPathogen", "EXCLUDE_UnassignedFamily")]))
  n_fams <- length(unique_fams)
  if (n_fams == 0) return(character(0))
  
  distinct_academic_palette <- c(
    "#332288", "#117733", "#CC6677", "#882255", "#44AA99", 
    "#88CCEE", "#DDCC77", "#AA4499", "#D55E00", "#0072B2", 
    "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#CC79A7", 
    "#125A56", "#A01813", "#4A237A", "#2B5C8F", "#5D8AA8"
  )
  
  if (n_fams <= length(distinct_academic_palette)) {
    pal <- distinct_academic_palette[1:n_fams]
  } else {
    set.seed(42)
    pal <- as.character(Polychrome::createPalette(n_fams, c("#2B5C8F", "#2D7044", "#D97706")))
  }
  
  names(pal) <- unique_fams
  return(pal)
}


# ==============================================================================
# --- FUNGUILD MYCOTYPE MAPPING HELPER ---
# ==============================================================================

exclude_families <- c("Aspergillaceae", "Mrakiaceae", "Trichosporonaceae", 
                      "Trimorphomycetaceae", "Malasseziaceae")
exclude_genera   <- c("Penicillium", "Tausonia", "Apiotrichum", "Saitozyma", 
                      "Malassezia", "Candida", "Cryptococcus", "Rhodotorula")

assign_mycotype_from_guild <- function(guild_str, trophic_str = "", family = "", genus = "") {
  fam <- coalesce(as.character(family), "")
  gen <- coalesce(as.character(genus), "")
  
  # Catch hardcoded exclusions first
  if (fam %in% exclude_families || gen %in% exclude_genera) {
    return("EXCLUDE_YeastMoldPathogen")
  }
  
  g <- paste(coalesce(guild_str, ""), coalesce(trophic_str, ""))
  
  if (g == "" || grepl("NULL|unassigned", g, ignore.case = TRUE) || grepl("Incertae_sedis", fam, ignore.case = TRUE)) {
    return("EXCLUDE_UnassignedFamily")
  }
  
  if (grepl("Orchid Mycorrhizal", g, ignore.case = TRUE)) {
    return("OMF")
  } else if (grepl("Ectomycorrhizal", g, ignore.case = TRUE)) {
    return("ECM")
  } else if (grepl("Ericoid Mycorrhizal", g, ignore.case = TRUE)) {
    return("ERM")
  } else if (grepl("Arbuscular Mycorrhizal", g, ignore.case = TRUE)) {
    return("AM")
  } else if (grepl("Saprotroph", g, ignore.case = TRUE)) {
    return("Saprotroph")
  } else if (grepl("Pathogen|Parasite|Wood Decay|Endophyte", g, ignore.case = TRUE)) {
    return("EXCLUDE_YeastMoldPathogen")
  }
  
  return("EXCLUDE_UnassignedFamily")
}


# --- 1. IMPORT DATA & MERGE ---
ps <- import_biom("otu_table_sintax.biom", refseqfilename = "repseqs.fasta")
meta_data <- import_qiime_sample_data("metadata_w_reads.txt")

ps2 <- merge_phyloseq(ps, meta_data)


# --- 2. TMM NORMALIZATION FUNCTION ---
tmm_normalize_phyloseq <- function(ps) {
  otu_obj <- otu_table(ps)
  is_taxa_rows <- taxa_are_rows(otu_obj)
  
  C <- otu_obj@.Data
  if (!is_taxa_rows) C <- t(C)
  
  C <- as.matrix(C)
  mode(C) <- "numeric"
  
  Y <- edgeR::DGEList(counts = C)
  Y <- edgeR::normLibSizes(Y, method = "TMM")
  eff <- Y$samples$lib.size * Y$samples$norm.factors
  
  mat_tmm    <- t(t(C) / eff) * mean(eff)
  mat_rel    <- t(t(C) / eff)
  mat_cpm    <- edgeR::cpm(Y, normalized.lib.sizes = TRUE)
  mat_logcpm <- edgeR::cpm(Y, normalized.lib.sizes = TRUE, log = TRUE, prior.count = 1)
  
  build_ps <- function(mat) {
    if (!is_taxa_rows) mat <- t(mat)
    new_otu <- otu_table(mat, taxa_are_rows = is_taxa_rows)
    tmp <- ps
    otu_table(tmp) <- new_otu
    return(tmp)
  }
  
  list(
    eff_lib_size  = eff,
    ps_tmm        = build_ps(mat_tmm),
    ps_tmm_rel    = build_ps(mat_rel),
    ps_tmm_cpm    = build_ps(mat_cpm),
    ps_tmm_logcpm = build_ps(mat_logcpm)
  )
}

res <- tmm_normalize_phyloseq(ps2)
ps_tmm       <- res$ps_tmm
ps_tmm_rel   <- res$ps_tmm_rel
eff.lib.size <- res$eff_lib_size


# --- 3. FILTER TO FUNGI ONLY ---
tt_mat <- as.matrix(tax_table(ps_tmm))

find_rank_col <- function(tt, pattern) {
  hits <- sapply(colnames(tt), function(cl) {
    sum(grepl(pattern, tt[, cl], ignore.case = TRUE), na.rm = TRUE)
  })
  if (all(hits == 0)) NA_character_ else names(which.max(hits))
}

kcol <- find_rank_col(tt_mat, '^"?k__'); stopifnot(!is.na(kcol))
kingdom_raw <- tt_mat[, kcol, drop = TRUE]

fungal_taxa <- taxa_names(ps_tmm)[!is.na(kingdom_raw) & grepl('^"?k__Fungi', kingdom_raw, ignore.case = TRUE)]
ps_fungi <- prune_taxa(fungal_taxa, ps_tmm)

cat("Start taxa:", ntaxa(ps_tmm), " -> Fungi-only:", ntaxa(ps_fungi), "\n")


# --- 4. DECONTAM & PREVALENCE FILTERING ---
df <- as.data.frame(sample_data(ps_fungi))
df$LibrarySize <- sample_sums(ps_fungi)
df <- df[order(df$LibrarySize),]
df$Index <- seq(nrow(df))

contamdf.freq <- isContaminant(ps_fungi, method="frequency", conc="reads")
ps.noncontam <- prune_taxa(!contamdf.freq$contaminant, ps_fungi)

otu <- as.matrix(otu_table(ps.noncontam))
if (!taxa_are_rows(ps.noncontam)) otu <- t(otu)
prev <- rowSums(otu > 0)
tot  <- rowSums(otu)

prev_min <- max(2, ceiling(0.02 * nsamples(ps.noncontam)))
tot_min  <- 10

keep <- prev >= prev_min & tot >= tot_min
ps_filt <- prune_taxa(keep, ps.noncontam)

cat("After filters:", ntaxa(ps_filt), "taxa\n")


# --- 5. STANDARDIZE TAXONOMY RANKS ---
tt <- as.matrix(tax_table(ps_filt))

strip_rank_prefix <- function(x) {
  x <- as.character(x)
  x <- gsub('^"|"$', "", x)
  x <- gsub("^[a-z]__", "", x, ignore.case = TRUE)
  x[x %in% c("", "Unassigned", "unassigned")] <- NA
  x
}

rank_patterns <- c(Kingdom='^"?k__', Phylum='^"?p__', Class='^"?c__',
                   Order='^"?o__', Family='^"?f__', Genus='^"?g__', Species='^"?s__')

new_names <- colnames(tt)
for (nm in names(rank_patterns)) {
  col <- find_rank_col(tt, rank_patterns[[nm]])
  if (!is.na(col)) new_names[new_names == col] <- nm
}
colnames(tt) <- new_names
for (cn in intersect(names(rank_patterns), colnames(tt))) tt[, cn] <- strip_rank_prefix(tt[, cn])

tax_table(ps_filt) <- tax_table(tt)


# --- 6. FUNGUILD ASSIGNMENT & MYCOTYPE CLASSIFICATION ---
fung_db <- FUNGuildR::get_funguild_db()

tt_amp <- as.data.frame(tax_table(ps_filt), stringsAsFactors = FALSE)
tt_amp <- tt_amp %>%
  mutate(Taxonomy = paste0(coalesce(Kingdom, ""), ";", coalesce(Phylum, ""), ";", 
                           coalesce(Class, ""), ";", coalesce(Order, ""), ";",
                           coalesce(Family, ""), ";", coalesce(Genus, ""), ";", 
                           coalesce(Species, "")))

funguild_res <- FUNGuildR::funguild_assign(tt_amp, db = fung_db)

tt_amp$Guild       <- funguild_res$guild
tt_amp$trophicMode <- funguild_res$trophicMode
tt_amp$MycoType_clean <- mapply(assign_mycotype_from_guild, 
                                tt_amp$Guild, 
                                tt_amp$trophicMode, 
                                tt_amp$Family, 
                                tt_amp$Genus)

tax_table(ps_filt) <- tax_table(as.matrix(tt_amp))

# Subset to retain valid mycorrhizae & family-resolved saprotrophs
ps_retained <- subset_taxa(ps_filt, !MycoType_clean %in% c("EXCLUDE_YeastMoldPathogen", "EXCLUDE_UnassignedFamily"))


# ==============================================================================
# --- 7. SUBSET TARGET HOST ---
# ==============================================================================

ps_target <- subset_samples(ps_retained, genusspecies == "C-aust")
ps_target <- prune_taxa(taxa_sums(ps_target) > 0, ps_target)


# --- 8. RELATIVE ABUNDANCE & SAVING ---
# PRESERVE FULL UN-TRIMMED PHYLOSEQ OBJECT FOR ALL DOWNSTREAM ANALYSES
ps_rel <- transform_sample_counts(ps_target, function(x) if (sum(x)==0) x else x / sum(x))
ps_rel <- prune_taxa(taxa_sums(ps_rel) > 0, ps_rel)

saveRDS(ps_rel, file = "funguild_Caust_ps_rel.rds")


# --- 9. BARPLOTS (AMPLICON FAMILY & MYCOTYPE) ---
# Filter out families < 1% relative abundance specifically for plotting/maps
plot_data <- psmelt(ps_rel) %>%
  mutate(
    FullTaxonomy  = paste(Kingdom, Phylum, Class, Order, Family, Genus, Species, sep = "; "),
    FullTaxonomy  = gsub("; NA", "", FullTaxonomy),
    MycoType_plot = factor(as.character(MycoType_clean)),
    MycoType_full = format_mycotype_label(MycoType_clean),
    Family_plot   = factor(as.character(Family))
  ) %>%
  filter(!MycoType_plot %in% c("EXCLUDE_YeastMoldPathogen", "EXCLUDE_UnassignedFamily")) %>%
  filter(!is.na(Family_plot) & !Family_plot %in% c("", "NA", "Unassigned", "Unassigned Family")) %>%
  filter(Abundance > 0)

# Filter out families < 1% relative abundance across the dataset for visual plots
fam_means <- plot_data %>%
  group_by(Family_plot) %>%
  summarise(Mean_Abund = mean(Abundance, na.rm = TRUE), .groups = "drop") %>%
  filter(Mean_Abund >= 0.01)

plot_data <- plot_data %>%
  filter(Family_plot %in% fam_means$Family_plot)

family_levels  <- sort(unique(as.character(plot_data$Family_plot)))
family_palette <- get_family_colors(family_levels)

# 9A. Barplot by Family
p_family <- ggplot(plot_data, aes(
  x    = sample,
  y    = Abundance,
  fill = Family_plot,
  text = paste0(
    "<b>Sample:</b> ", sample,
    "<br><b>Abundance:</b> ", round(Abundance, 3),
    "<br><b>Identification:</b> ", FullTaxonomy,
    "<br><b>Family:</b> ", Family_plot,
    "<br><b>MycoType:</b> ", MycoType_full
  )
)) +
  geom_bar(stat = "identity", position = "stack", linewidth = 0) +
  scale_fill_manual(
    values = family_palette, 
    breaks = names(family_palette),
    drop   = TRUE
  ) +
  scale_x_discrete(drop = TRUE) +
  facet_wrap(~county, scales = "free_x", ncol = 4) +
  theme_bw(base_size = 10) +
  theme(
    axis.text.x         = element_blank(), 
    axis.ticks.x        = element_blank(),
    legend.key.size     = unit(0.5, "cm"),
    legend.key.height   = unit(0.5, "cm"),
    legend.key.width    = unit(0.5, "cm"),
    legend.spacing.y    = unit(0.25, "cm"),
    legend.spacing.x    = unit(0.2, "cm"),
    legend.title        = element_text(size = 11, face = "bold"),
    legend.text         = element_text(size = 9.5)
  ) +
  guides(fill = guide_legend(ncol = 3, byrow = TRUE, override.aes = list(alpha = 1), title = "Family")) +
  labs(title = NULL, y = "Relative abundance", x = NULL, fill = "Family")

pp_family <- ggplotly(p_family, tooltip = "text")
save_barplot_all_formats(p_family, pp_family, "fungal_family_barplot", dir_amplicon_bar)

# 9B. Barplot by MycoType
p_mycotype <- ggplot(plot_data, aes(
  x    = sample,
  y    = Abundance,
  fill = MycoType_plot,
  text = paste0(
    "<b>Sample:</b> ", sample,
    "<br><b>Abundance:</b> ", round(Abundance, 3),
    "<br><b>Identification:</b> ", FullTaxonomy,
    "<br><b>MycoType:</b> ", MycoType_full
  )
)) +
  geom_bar(stat = "identity", position = "stack", linewidth = 0) +
  scale_fill_manual(values = myco_colors) +
  scale_x_discrete(drop = TRUE) +
  facet_wrap(~county, scales = "free_x", ncol = 4) +
  theme_bw(base_size = 10) +
  theme(
    axis.text.x         = element_blank(), 
    axis.ticks.x        = element_blank(),
    legend.key.size     = unit(0.5, "cm"),
    legend.key.height   = unit(0.5, "cm"),
    legend.spacing.y    = unit(0.2, "cm"),
    legend.title        = element_text(size = 11, face = "bold"),
    legend.text         = element_text(size = 9.5)
  ) +
  labs(title = NULL, y = "Relative abundance", x = NULL, fill = "MycoType")

pp_mycotype <- ggplotly(p_mycotype, tooltip = "text")
save_barplot_all_formats(p_mycotype, pp_mycotype, "fungal_mycotype_barplot", dir_amplicon_bar)


# ==============================================================================
# --- 10. SPATIAL PIE CHART MAPS & COMPOUND OBJECT GROUPING ---
# ==============================================================================

separate_pies <- function(df_coords, lon_col = "Pie_Lon", lat_col = "Pie_Lat", min_dist = 1.80, origin_clearance = 0.55, step_size = 0.18, max_iter = 4000) {
  df <- df_coords
  if (!all(c(lon_col, lat_col) %in% names(df))) return(df)
  
  n <- nrow(df)
  if (n <= 1) return(df)
  
  for (iter in seq_len(max_iter)) {
    moved <- FALSE
    
    # 1. Repel pie centers from each other
    for (i in 1:(n - 1)) {
      for (j in (i + 1):n) {
        dx <- df[[lon_col]][j] - df[[lon_col]][i]
        dy <- df[[lat_col]][j] - df[[lat_col]][i]
        dist <- sqrt(dx^2 + dy^2)
        
        if (dist < min_dist) {
          moved <- TRUE
          if (dist == 0) {
            dx <- runif(1, -0.03, 0.03)
            dy <- runif(1, -0.03, 0.03)
            dist <- sqrt(dx^2 + dy^2)
          }
          overlap <- min_dist - dist
          off_x <- (dx / dist) * (overlap / 2) * step_size
          off_y <- (dy / dist) * (overlap / 2) * step_size
          
          df[[lon_col]][i] <- df[[lon_col]][i] - off_x
          df[[lat_col]][i] <- df[[lat_col]][i] - off_y
          df[[lon_col]][j] <- df[[lon_col]][j] + off_x
          df[[lat_col]][j] <- df[[lat_col]][j] + off_y
        }
      }
    }
    
    # 2. Maintain strict clearance from original collection dots
    for (i in 1:n) {
      for (k in 1:n) {
        dx_orig <- df[[lon_col]][i] - df$lon[k]
        dy_orig <- df[[lat_col]][i] - df$lat[k]
        dist_orig <- sqrt(dx_orig^2 + dy_orig^2)
        
        if (dist_orig < origin_clearance) {
          moved <- TRUE
          if (dist_orig == 0) {
            dx_orig <- runif(1, 0.01, 0.03)
            dy_orig <- runif(1, 0.01, 0.03)
            dist_orig <- sqrt(dx_orig^2 + dy_orig^2)
          }
          push <- (origin_clearance - dist_orig) * step_size
          df[[lon_col]][i] <- df[[lon_col]][i] + (dx_orig / dist_orig) * push
          df[[lat_col]][i] <- df[[lat_col]][i] + (dy_orig / dist_orig) * push
        }
      }
    }
    
    if (!moved) break
  }
  return(df)
}

build_spatial_map <- function(ps_obj, target_var, base_filename, out_dir, title_text = NULL, color_palette = NULL) {
  
  df_map <- psmelt(ps_obj)
  
  colnames(df_map) <- gsub("^longitude$", "lon", colnames(df_map), ignore.case = TRUE)
  colnames(df_map) <- gsub("^latitude$", "lat", colnames(df_map), ignore.case = TRUE)
  
  # Coordinate hemisphere fix if longitudes were provided as positive numbers
  if (all(df_map$lon > 0, na.rm = TRUE)) df_map$lon <- -abs(df_map$lon)
  
  df_map <- df_map %>%
    mutate(
      Family_plot   = Family,
      MycoType_plot = MycoType_clean
    ) %>%
    filter(!MycoType_plot %in% c("EXCLUDE_YeastMoldPathogen", "EXCLUDE_UnassignedFamily")) %>%
    filter(!is.na(Family_plot) & !Family_plot %in% c("", "NA", "Unassigned", "Unassigned Family"))
  
  # Completely drop families < 1% relative abundance for clear spatial map visuals
  if (target_var != "MycoType_clean") {
    fam_map_means <- df_map %>%
      group_by(Family_plot) %>%
      summarise(Mean_Abund = mean(Abundance, na.rm = TRUE), .groups = "drop") %>%
      filter(Mean_Abund >= 0.01)
    
    df_map <- df_map %>% filter(Family_plot %in% fam_map_means$Family_plot)
  }
  
  actual_target_str <- if (target_var == "MycoType_clean") "MycoType_plot" else target_var
  target_sym        <- rlang::sym(actual_target_str)
  
  pie_long <- df_map %>%
    filter(!is.na(!!target_sym), !!target_sym != "") %>%
    group_by(sample, county, lon, lat, !!target_sym) %>%
    summarise(Abundance = sum(Abundance), .groups = "drop")
  
  pie_data <- pie_long %>%
    tidyr::pivot_wider(
      names_from  = !!target_sym,
      values_from = Abundance,
      values_fill = 0
    ) %>%
    filter(!is.na(lon), !is.na(lat)) %>%
    group_by(county, lon, lat) %>%
    summarise(
      across(where(is.numeric) & !any_of(c("lon", "lat")), ~ mean(.x, na.rm = TRUE)),
      n_samples = n(),
      .groups = "drop"
    )
  
  pie_cols <- setdiff(names(pie_data), c("county", "lon", "lat", "n_samples"))
  
  # CULL SITES/PIES WITH ZERO RETAINED TAXA ABUNDANCE
  pie_data$Total_Abundance <- rowSums(pie_data[, pie_cols, drop = FALSE])
  pie_data <- pie_data %>% filter(Total_Abundance > 0)
  
  pie_data <- pie_data %>%
    mutate(
      Label   = paste0(county, " (", n_samples, ")"),
      fixed_r = 0.40,
      Pie_Lon = lon, 
      Pie_Lat = lat
    )
  
  # Disambiguate overlapping pie positions
  pie_data <- separate_pies(pie_data, lon_col = "Pie_Lon", lat_col = "Pie_Lat", min_dist = 1.80, origin_clearance = 0.55, step_size = 0.18, max_iter = 4000)
  
  # Lock county label atop pie center as a single compound visual unit
  pie_data <- pie_data %>%
    mutate(Label_Lat = Pie_Lat + fixed_r + 0.12)
  
  world     <- ne_countries(scale = "medium", returnclass = "sf")
  states    <- ne_states(country = "United States of America", returnclass = "sf")
  provinces <- ne_states(country = "Canada", returnclass = "sf")
  
  legend_title_text <- ifelse(target_var == "MycoType_clean", "MycoType", "Family")
  
  # Dynamic coordinate limits based on active pie data
  min_lon <- min(pie_data$lon, na.rm = TRUE) - 2.0
  max_lon <- max(pie_data$lon, na.rm = TRUE) + 2.0
  min_lat <- min(pie_data$lat, na.rm = TRUE) - 2.0
  max_lat <- max(pie_data$lat, na.rm = TRUE) + 2.0
  
  map_plot <- ggplot() +
    geom_sf(data = world,     fill = "gray95", color = "black",  linewidth = 0.3) +
    geom_sf(data = states,    fill = NA,       color = "gray60", linewidth = 0.4) +
    geom_sf(data = provinces, fill = NA,       color = "gray60", linewidth = 0.4) +
    # 1. Collection site origin point
    geom_point(data = pie_data, aes(x = lon, y = lat), size = 2.0, color = "black") +
    # 2. Leader segment connecting geographic origin dot to displaced pie chart
    geom_segment(
      data = pie_data, 
      aes(x = lon, y = lat, xend = Pie_Lon, yend = Pie_Lat), 
      color = "gray30", 
      linewidth = 0.4, 
      linetype = "dashed"
    ) +
    # 3. Displaced Pie Chart
    geom_scatterpie(
      aes(x = Pie_Lon, y = Pie_Lat, r = fixed_r), 
      data = pie_data, 
      cols = pie_cols, 
      color = NA, 
      alpha = 0.95
    ) +
    # 4. Locked County Label Header (sharing Pie_Lon coordinate anchor with pie)
    geom_text(
      data = pie_data, 
      aes(x = Pie_Lon, y = Label_Lat, label = Label), 
      size = 2.8, 
      fontface = "bold"
    ) +
    coord_sf(xlim = c(min_lon, max_lon), ylim = c(min_lat, max_lat), expand = FALSE) +
    theme_minimal() +
    theme(
      panel.grid = element_blank(),
      axis.title = element_blank(),
      legend.key.size = unit(0.4, "cm")
    ) +
    labs(title = title_text, fill = legend_title_text)
  
  if (!is.null(color_palette)) {
    map_plot <- map_plot + scale_fill_manual(values = color_palette, breaks = pie_cols)
  }
  
  save_map_all_formats(map_plot, base_filename, out_dir)
  return(map_plot)
}

# 10A. Amplicon Spatial Map by Family
fam_vec_map <- as.data.frame(tax_table(ps_rel))$Family
fam_map_colors <- get_family_colors(fam_vec_map[!is.na(fam_vec_map)])

map_family <- build_spatial_map(
  ps_obj        = ps_rel,
  target_var    = "Family_plot",
  base_filename = "fungal_family_map",
  out_dir       = dir_amplicon_map,
  color_palette = fam_map_colors
)

# 10B. Amplicon Spatial Map by MycoType
map_mycotype <- build_spatial_map(
  ps_obj        = ps_rel,
  target_var    = "MycoType_clean",
  base_filename = "fungal_mycotype_map",
  out_dir       = dir_amplicon_map,
  color_palette = myco_colors
)


# ==============================================================================
# SECTION 11: STANDALONE PELOTON BLAST PIPELINE (STACKED BY REPLICATE)
# ==============================================================================

peloton_csv_path <- "PelotonBLAST_with_metadata.csv"

# 11.1) Import Peloton BLAST Data & Extract Replicate IDs
peloton_raw <- read.csv(peloton_csv_path, stringsAsFactors = FALSE)

peloton_clean <- peloton_raw %>%
  mutate(
    Pct_Identity = as.numeric(str_replace(`X..Pairwise.Identity`, "%", "")),
    Query_ID     = str_trim(str_replace(Query, "Assembly", "")),
    # Extract Replicate ID (e.g., '056aFUNG1' -> '056a')
    Replicate_ID = sub("FUNG.*$", "", Query_ID),
    Latitude     = as.numeric(Latitude),
    Longitude    = as.numeric(Longitude)
  )

# 11.2) UNITE Taxonomy Parsing & Cleaning
peloton_tax <- peloton_clean %>%
  mutate(TaxString = sapply(strsplit(Name, "\\|"), function(x) x[length(x)])) %>%
  separate(
    TaxString,
    into = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"),
    sep  = ";",
    fill = "right"
  ) %>%
  mutate(across(c(Kingdom:Species), ~ {
    x <- as.character(.x)
    x <- gsub('^"|"$', "", x)
    x <- gsub("^[kpcofgs]__", "", x, ignore.case = TRUE)
    x <- gsub("_sp$", "", x, ignore.case = TRUE)
    x <- trimws(x)
    x[x %in% c("", "NA", "Unassigned", "unassigned")] <- NA
    x
  }))

# Apply FUNGuild Assignment to Peloton Dataset
peloton_tax <- peloton_tax %>%
  mutate(Taxonomy = paste0(coalesce(Kingdom, ""), ";", coalesce(Phylum, ""), ";", 
                           coalesce(Class, ""), ";", coalesce(Order, ""), ";",
                           coalesce(Family, ""), ";", coalesce(Genus, ""), ";", 
                           coalesce(Species, "")))

peloton_funguild        <- FUNGuildR::funguild_assign(peloton_tax, db = fung_db)
peloton_tax$Guild       <- peloton_funguild$guild
peloton_tax$trophicMode <- peloton_funguild$trophicMode

peloton_tax$MycoType_clean <- mapply(assign_mycotype_from_guild, 
                                     peloton_tax$Guild, 
                                     peloton_tax$trophicMode, 
                                     peloton_tax$Family, 
                                     peloton_tax$Genus)

# Filter retained taxa
peloton_retained <- peloton_tax %>%
  filter(!MycoType_clean %in% c("EXCLUDE_YeastMoldPathogen", "EXCLUDE_UnassignedFamily")) %>%
  filter(!is.na(Family) & !Family %in% c("", "NA", "Unassigned", "Unassigned Family"))

# 11.3) Summarize Fungi Count Per Replicate
replicate_counts <- peloton_retained %>%
  group_by(Replicate_ID) %>%
  summarise(N_Retained_Fungi = n(), .groups = "drop")

df_pel_stacked <- peloton_retained %>%
  left_join(replicate_counts, by = "Replicate_ID") %>%
  mutate(
    Family_plot   = factor(as.character(Family)),
    MycoType_plot = factor(as.character(MycoType_clean)),
    MycoType_full = format_mycotype_label(MycoType_clean),
    FullTaxonomy  = paste(Kingdom, Phylum, Class, Order, Family, Genus, Species, sep = "; ")
  )

pel_family_palette <- get_family_colors(levels(df_pel_stacked$Family_plot))

# 11.4) Peloton Barplots Stacked by Replicate (Family - Clean Axis)
p_pel_family <- ggplot(df_pel_stacked, aes(
  x    = Replicate_ID,
  fill = Family_plot,
  text = paste0(
    "<b>Replicate:</b> ", Replicate_ID,
    "<br><b>Retained Fungi in Replicate:</b> ", N_Retained_Fungi,
    "<br><b>Identification:</b> ", FullTaxonomy,
    "<br><b>Family:</b> ", Family_plot,
    "<br><b>MycoType:</b> ", MycoType_full
  )
)) +
  geom_bar(position = "stack", linewidth = 0) +
  scale_fill_manual(values = pel_family_palette) +
  scale_x_discrete(drop = TRUE) +
  facet_wrap(~County, scales = "free_x", ncol = 4) +
  theme_bw(base_size = 10) +
  theme(
    axis.text.x         = element_blank(),
    axis.ticks.x        = element_blank(),
    legend.key.size     = unit(0.5, "cm"),
    legend.key.height   = unit(0.5, "cm"),
    legend.key.width    = unit(0.5, "cm"),
    legend.spacing.y    = unit(0.25, "cm"),
    legend.spacing.x    = unit(0.2, "cm"),
    legend.title        = element_text(size = 11, face = "bold"),
    legend.text         = element_text(size = 9.5)
  ) +
  guides(fill = guide_legend(ncol = 2, byrow = TRUE, override.aes = list(alpha = 1))) +
  labs(title = NULL, y = "Number of Fungal OTUs", x = NULL, fill = "Family")

pp_pel_family <- ggplotly(p_pel_family, tooltip = "text")
save_barplot_all_formats(p_pel_family, pp_pel_family, "peloton_family_barplot", dir_peloton_bar)

# 11.5) Peloton Barplots Stacked by Replicate (MycoType - Clean Axis)
p_pel_mycotype <- ggplot(df_pel_stacked, aes(
  x    = Replicate_ID,
  fill = MycoType_plot,
  text = paste0(
    "<b>Replicate:</b> ", Replicate_ID,
    "<br><b>Retained Fungi in Replicate:</b> ", N_Retained_Fungi,
    "<br><b>Identification:</b> ", FullTaxonomy,
    "<br><b>MycoType:</b> ", MycoType_full
  )
)) +
  geom_bar(position = "stack", linewidth = 0) +
  scale_fill_manual(values = myco_colors) +
  scale_x_discrete(drop = TRUE) +
  facet_wrap(~County, scales = "free_x", ncol = 4) +
  theme_bw(base_size = 10) +
  theme(
    axis.text.x         = element_blank(),
    axis.ticks.x        = element_blank(),
    legend.key.size     = unit(0.5, "cm"),
    legend.key.height   = unit(0.5, "cm"),
    legend.spacing.y    = unit(0.2, "cm"),
    legend.title        = element_text(size = 11, face = "bold"),
    legend.text         = element_text(size = 9.5)
  ) +
  labs(title = NULL, y = "Number of Fungal OTUs", x = NULL, fill = "MycoType")

pp_pel_mycotype <- ggplotly(p_pel_mycotype, tooltip = "text")
save_barplot_all_formats(p_pel_mycotype, pp_pel_mycotype, "peloton_mycotype_barplot", dir_peloton_bar)

# --- 11.5.5) CONSTRUCT PHYLOSEQ OBJECT FOR PELOTON MAPS ---
otu_pel <- df_pel_stacked %>%
  count(Replicate_ID, Query_ID) %>%
  pivot_wider(names_from = Replicate_ID, values_from = n, values_fill = 0) %>%
  column_to_rownames("Query_ID") %>%
  as.matrix() %>%
  otu_table(taxa_are_rows = TRUE)

tax_pel <- df_pel_stacked %>%
  distinct(Query_ID, .keep_all = TRUE) %>%
  select(Query_ID, Kingdom, Phylum, Class, Order, Family, Genus, Species, MycoType_clean) %>%
  column_to_rownames("Query_ID") %>%
  as.matrix() %>%
  tax_table()

sam_pel <- df_pel_stacked %>%
  distinct(Replicate_ID, .keep_all = TRUE) %>%
  select(Replicate_ID, County, Latitude, Longitude) %>%
  mutate(sample = Replicate_ID) %>%
  rename(county = County, latitude = Latitude, longitude = Longitude) %>%
  column_to_rownames("Replicate_ID") %>%
  sample_data()

ps_peloton <- phyloseq(otu_pel, tax_pel, sam_pel)

# 11.6) Peloton Spatial Maps
map_pel_family <- build_spatial_map(
  ps_obj        = ps_peloton,
  target_var    = "Family_plot",
  base_filename = "peloton_family_map",
  out_dir       = dir_peloton_map,
  color_palette = pel_family_palette
)

map_pel_mycotype <- build_spatial_map(
  ps_obj        = ps_peloton,
  target_var    = "MycoType_clean",
  base_filename = "peloton_mycotype_map",
  out_dir       = dir_peloton_map,
  color_palette = myco_colors
)

# SAVE THE CORRECT PELOTON PHYLOSEQ OBJECT
saveRDS(ps_peloton, "C:/Users/withe/Desktop/C_aust_integrated_analysis/peloton_Caust_ps_rel.rds")

# ==============================================================================
# GENERATE SPECIES-LEVEL TAXA COUNT & OCCURRENCE TABLE
# ==============================================================================

# Extract tax_table and otu_table from your phyloseq object (ps_rel or ps_target)
tax_df  <- as.data.frame(as(tax_table(ps_rel), "matrix"), stringsAsFactors = FALSE)
otu_mat <- as.matrix(as(otu_table(ps_rel), "matrix"))
if (taxa_are_rows(ps_rel)) otu_mat <- t(otu_mat)

# Calculate sample occurrences (> 0 abundance) and total abundance per ASV
asv_summary <- tax_df %>%
  tibble::rownames_to_column("ASV_ID") %>%
  mutate(
    # Number of host samples in which this ASV was detected
    Sample_Occurrence_Count = colSums(otu_mat[, ASV_ID, drop = FALSE] > 0),
    # Sum of relative abundance across all samples
    Total_Relative_Abundance = colSums(otu_mat[, ASV_ID, drop = FALSE])
  )

# Aggregate to Species / Taxonomy Level
species_taxa_table <- asv_summary %>%
  group_by(Phylum, Class, Order, Family, Genus, Species, MycoType_clean) %>%
  summarise(
    Number_of_ASVs           = n(),
    Total_Sample_Occurrences = sum(Sample_Occurrence_Count),
    Mean_Relative_Abundance  = round(mean(Total_Relative_Abundance), 4),
    Total_Relative_Abundance = round(sum(Total_Relative_Abundance), 4),
    .groups = "drop"
  ) %>%
  arrange(desc(Total_Sample_Occurrences), desc(Total_Relative_Abundance))

# Print top 15 species in console
print(head(species_taxa_table, 15))

# Export to CSV file
write.csv(species_taxa_table, file.path(main_out_dir, "Fungal_Species_Counts_and_Occurrences.csv"), row.names = FALSE)



# Write all loaded packages in your current session directly to a .bib file
knitr::write_bib(names(sessionInfo()$otherPkgs), file = "chapter2_packages.bib")
