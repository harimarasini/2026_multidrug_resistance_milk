# ============================================================
#   AST ANALYSIS — LOAD DIRECTLY FROM CSV FILE
#   Handles both WIDE and LONG format CSV inputs
#   Just change the file path and run!
# ============================================================

# Install packages if needed (run once):
# install.packages(c("ggplot2", "dplyr", "tidyr", "RColorBrewer"))

library(ggplot2)
library(dplyr)
library(tidyr)


# ================================================================
# STEP 1: SET YOUR FILE PATH HERE
# ================================================================

file_path <- "2026_ast_data_amik.csv"   # <-- Change to your actual path
                                          # e.g., "C:/Users/You/Desktop/2026_ast_data_amik.csv"
                                          # Mac/Linux: "/home/you/data/2026_ast_data_amik.csv"


# ================================================================
# STEP 2: DETECT FORMAT AND LOAD
# ================================================================

raw <- read.csv(file_path, stringsAsFactors = FALSE, na.strings = c("", "NA"))

# Clean column names (remove spaces, fix case)
colnames(raw) <- trimws(colnames(raw))

cat("✔ File loaded:", nrow(raw), "rows x", ncol(raw), "columns\n")
cat("  Column names:", paste(colnames(raw), collapse = ", "), "\n\n")

# --- Auto-detect format ---
# Long format: has a column named "Result" or "SIR" or "Susceptibility"
# Wide format: antibiotics are column names

is_long_format <- any(tolower(colnames(raw)) %in% c("result", "sir", "susceptibility", "category"))

if (is_long_format) {
  
  # ── LONG FORMAT HANDLING ──────────────────────────────────
  cat("📋 Detected: LONG format\n\n")
  
  # Rename columns flexibly
  colnames(raw)[tolower(colnames(raw)) == "organism"]       <- "Organism"
  colnames(raw)[tolower(colnames(raw)) == "isolate_id" |
                tolower(colnames(raw)) == "id" |
                tolower(colnames(raw)) == "sample_id"]      <- "Isolate_ID"
  colnames(raw)[tolower(colnames(raw)) == "antibiotic" |
                tolower(colnames(raw)) == "drug" |
                tolower(colnames(raw)) == "antimicrobial"]  <- "Antibiotic"
  colnames(raw)[tolower(colnames(raw)) %in%
                  c("result","sir","susceptibility","category")] <- "Result"
  
  ast_long <- raw %>%
    mutate(Result = toupper(trimws(Result))) %>%
    filter(Result %in% c("S","I","R"))
  
} else {
  
  # ── WIDE FORMAT HANDLING ──────────────────────────────────
  cat("📋 Detected: WIDE format — converting to long...\n\n")
  
  # Identify the Organism and ID columns
  org_col <- colnames(raw)[tolower(colnames(raw)) %in%
               c("organism","species","bacteria","pathogen")][1]
  id_col  <- colnames(raw)[tolower(colnames(raw)) %in%
               c("isolate_id","id","sample_id","no","number","isolate")][1]
  
  if (is.na(org_col)) stop("❌ Could not find 'Organism' column. Please name it 'Organism' in your CSV.")
  
  # All remaining columns are antibiotics
  ab_cols <- setdiff(colnames(raw), c(org_col, id_col))
  
  # Pivot to long format
  ast_long <- raw %>%
    rename(Organism = all_of(org_col)) %>%
    { if (!is.na(id_col)) rename(., Isolate_ID = all_of(id_col)) else mutate(., Isolate_ID = row_number()) } %>%
    pivot_longer(cols = all_of(ab_cols),
                 names_to  = "Antibiotic",
                 values_to = "Result") %>%
    mutate(Result = toupper(trimws(Result))) %>%
    filter(Result %in% c("S","I","R"))
}

# Final clean
ast_long <- ast_long %>%
  mutate(
    Organism   = trimws(Organism),
    Antibiotic = trimws(Antibiotic),
    Result     = factor(Result, levels = c("S","I","R"), ordered = TRUE)
  )

cat("✔ Data ready:", nrow(ast_long), "total observations\n")
cat("  Organisms  :", paste(unique(ast_long$Organism), collapse = " | "), "\n")
cat("  Antibiotics:", paste(unique(ast_long$Antibiotic), collapse = ", "), "\n\n")


# ================================================================
# STEP 3: RESISTANCE SUMMARY TABLE
# ================================================================

cat(strrep("=", 65), "\n")
cat("  RESISTANCE PROFILES PER ORGANISM\n")
cat(strrep("=", 65), "\n")

summary_table <- ast_long %>%
  group_by(Organism, Antibiotic) %>%
  summarise(
    n   = n(),
    S_n = sum(Result == "S"),
    I_n = sum(Result == "I"),
    R_n = sum(Result == "R"),
    S_pct       = round(S_n / n * 100, 1),
    I_pct       = round(I_n / n * 100, 1),
    R_pct       = round(R_n / n * 100, 1),
    NonSusc_pct = round((I_n + R_n) / n * 100, 1),
    .groups = "drop"
  )

# Print nicely per organism
for (org in unique(summary_table$Organism)) {
  sub <- summary_table %>% filter(Organism == org)
  cat("\nOrganism:", org, "| n =", unique(sub$n)[1], "\n")
  cat(sprintf("  %-20s %5s %7s  %5s %7s  %5s %7s  %9s\n",
              "Antibiotic","S","S%","I","I%","R","R%","NonSusc%"))
  cat("  ", strrep("-", 60), "\n", sep = "")
  for (i in 1:nrow(sub)) {
    cat(sprintf("  %-20s %5d %6.1f%%  %5d %6.1f%%  %5d %6.1f%%  %8.1f%%\n",
                sub$Antibiotic[i],
                sub$S_n[i], sub$S_pct[i],
                sub$I_n[i], sub$I_pct[i],
                sub$R_n[i], sub$R_pct[i],
                sub$NonSusc_pct[i]))
  }
}


# ================================================================
# STEP 4: MDR ASSESSMENT (R to >= 3 antibiotics)
# ================================================================

cat("\n", strrep("=", 65), "\n", sep = "")
cat("  MULTI-DRUG RESISTANCE (MDR) ASSESSMENT\n")
cat(strrep("=", 65), "\n")

mdr_table <- ast_long %>%
  group_by(Organism, Isolate_ID) %>%
  summarise(R_count = sum(Result == "R"), .groups = "drop") %>%
  group_by(Organism) %>%
  summarise(
    n       = n(),
    MDR_n   = sum(R_count >= 3),
    MDR_pct = round(MDR_n / n * 100, 1),
    .groups = "drop"
  )

print(mdr_table, row.names = FALSE)


# ================================================================
# STEP 5: STATISTICAL TESTS — FISHER'S EXACT / CHI-SQUARE
#         (Comparing organisms pairwise per antibiotic)
# ================================================================

cat("\n", strrep("=", 65), "\n", sep = "")
cat("  INTER-ORGANISM STATISTICAL COMPARISON\n")
cat("  Fisher's Exact (expected < 5) | Chi-Square (expected >= 5)\n")
cat(strrep("=", 65), "\n")

organisms <- unique(ast_long$Organism)
pairs <- combn(organisms, 2, simplify = FALSE)

all_test_results <- list()

for (pair in pairs) {
  org1 <- pair[1]; org2 <- pair[2]
  d1 <- ast_long %>% filter(Organism == org1)
  d2 <- ast_long %>% filter(Organism == org2)
  shared_abs <- intersect(unique(d1$Antibiotic), unique(d2$Antibiotic))
  
  if (length(shared_abs) == 0) next
  
  cat("\n  ┌──", org1, "vs", org2, "\n")
  cat(sprintf("  %-20s %-15s %-12s %-10s %s\n",
              "Antibiotic","Test","Statistic","p-value","Sig."))
  cat("  ", strrep("-", 62), "\n", sep="")
  
  for (ab in shared_abs) {
    v1 <- factor(d1$Result[d1$Antibiotic == ab], levels = c("S","I","R"))
    v2 <- factor(d2$Result[d2$Antibiotic == ab], levels = c("S","I","R"))
    tab <- rbind(table(v1), table(v2))
    tab <- tab[, colSums(tab) > 0, drop = FALSE]
    
    if (ncol(tab) < 2) {
      cat(sprintf("  %-20s %-15s %-12s %-10s\n", ab, "N/A", "—", "uniform"))
      next
    }
    
    expected <- outer(rowSums(tab), colSums(tab)) / sum(tab)
    
    if (any(expected < 5)) {
      res       <- fisher.test(tab, simulate.p.value = TRUE, B = 10000)
      test_name <- "Fisher's Exact"
      stat_str  <- "—"
      p_val     <- res$p.value
    } else {
      res       <- chisq.test(tab)
      test_name <- "Chi-Square"
      stat_str  <- paste0("X2=", round(res$statistic, 2))
      p_val     <- res$p.value
    }
    
    sig   <- ifelse(p_val < 0.001, "***", ifelse(p_val < 0.01, "**",
             ifelse(p_val < 0.05, "*", "ns")))
    p_str <- ifelse(p_val < 0.001, "<0.001", round(p_val, 3))
    
    cat(sprintf("  %-20s %-15s %-12s %-10s %s\n",
                ab, test_name, stat_str, p_str, sig))
    
    all_test_results[[paste(org1, org2, ab, sep="|")]] <- 
      data.frame(Org1=org1, Org2=org2, Antibiotic=ab,
                 Test=test_name, p_value=p_val, Sig=sig)
  }
}


# ================================================================
# STEP 6: KRUSKAL-WALLIS — Within-organism antibiotic comparison
# ================================================================

cat("\n", strrep("=", 65), "\n", sep = "")
cat("  KRUSKAL-WALLIS TEST (Within-organism, across antibiotics)\n")
cat("  Ordinal scoring: S=0, I=1, R=2\n")
cat(strrep("=", 65), "\n\n")

ast_long <- ast_long %>%
  mutate(Score = as.integer(Result) - 1)  # S=0, I=1, R=2

for (org in unique(ast_long$Organism)) {
  sub <- ast_long %>% filter(Organism == org)
  ab_list <- unique(sub$Antibiotic)
  
  # Only keep antibiotics with variation
  groups <- lapply(ab_list, function(ab) sub$Score[sub$Antibiotic == ab])
  vary   <- sapply(groups, function(g) length(unique(g)) > 1)
  groups <- groups[vary]
  labels <- ab_list[vary]
  
  if (length(groups) < 2) {
    cat(sprintf("  %-22s: Insufficient variation\n", org))
    next
  }
  
  res   <- kruskal.test(groups)
  sig   <- ifelse(res$p.value < 0.001, "***", ifelse(res$p.value < 0.01, "**",
           ifelse(res$p.value < 0.05, "*", "ns")))
  p_str <- ifelse(res$p.value < 0.001, "<0.001",
                  formatC(res$p.value, format="f", digits=4))
  
  cat(sprintf("  %-22s: H=%.2f, df=%d, p=%s %s\n",
              org, res$statistic, res$parameter, p_str, sig))
  cat(sprintf("  %s Tested: %s\n\n", strrep(" ", 24), paste(labels, collapse=", ")))
}


# ================================================================
# STEP 7: VISUALISATION — Stacked Bar Chart
# ================================================================

plot_data <- summary_table %>%
  select(Organism, Antibiotic, S_pct, I_pct, R_pct) %>%
  pivot_longer(cols = c(S_pct, I_pct, R_pct),
               names_to  = "Category",
               values_to = "Percentage") %>%
  mutate(Category = recode(Category,
                           "S_pct" = "S",
                           "I_pct" = "I",
                           "R_pct" = "R"),
         Category = factor(Category, levels = c("R","I","S")))

sir_colors <- c("S" = "#2ecc71", "I" = "#f39c12", "R" = "#e74c3c")

p <- ggplot(plot_data, aes(x = Antibiotic, y = Percentage, fill = Category)) +
  geom_bar(stat = "identity", position = "stack", color = "white", linewidth = 0.3) +
  facet_wrap(~Organism, scales = "free_x", ncol = 2) +
  scale_fill_manual(values = sir_colors,
                    labels = c("S"="Sensitive","I"="Intermediate","R"="Resistant"),
                    name   = "Susceptibility") +
  labs(
    title    = "Antimicrobial Susceptibility Profiles",
    subtitle = paste("Data source:", basename(file_path)),
    x = "Antibiotic", y = "Percentage of Isolates (%)"
  ) +
  theme_bw(base_size = 11) +
  theme(
    axis.text.x      = element_text(angle = 40, hjust = 1, size = 9),
    strip.background = element_rect(fill = "#2c3e50"),
    strip.text       = element_text(color = "white", face = "bold"),
    legend.position  = "bottom",
    plot.title       = element_text(face = "bold", size = 13)
  )

print(p)
# To save the plot:
# ggsave("AST_Plot.png", p, width = 12, height = 9, dpi = 300)


# ================================================================
# STEP 8: EXPORT RESULTS TO CSV
# ================================================================

write.csv(summary_table, "AST_Resistance_Summary.csv", row.names = FALSE)
write.csv(mdr_table,     "AST_MDR_Summary.csv",        row.names = FALSE)

if (length(all_test_results) > 0) {
  test_df <- do.call(rbind, all_test_results)
  write.csv(test_df, "AST_Statistical_Tests.csv", row.names = FALSE)
}

cat("\n✔ Results exported:\n")
cat("  → AST_Resistance_Summary.csv\n")
cat("  → AST_MDR_Summary.csv\n")
cat("  → AST_Statistical_Tests.csv\n\n")
cat(strrep("=", 65), "\n")
cat("Significance: *** p<0.001 | ** p<0.01 | * p<0.05 | ns = not significant\n")
cat(strrep("=", 65), "\n")
