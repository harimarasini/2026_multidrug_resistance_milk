# ============================================================
#   ANTIMICROBIAL SUSCEPTIBILITY TEST (AST) - FULL ANALYSIS
#   Statistical Tests: Fisher's Exact | Chi-Square | Kruskal-Wallis
#   For use in RStudio
# ============================================================

# Install packages if needed (uncomment first time):
# install.packages(c("ggplot2", "reshape2", "RColorBrewer", "dplyr"))

library(ggplot2)
library(reshape2)
library(dplyr)


# ================================================================
# SECTION 1: DATA ENTRY
# ================================================================

s_uberis <- data.frame(
  ID          = 1:18,
  Penicillin  = rep("S", 18),
  Tetracycline= c("I","R","I","R","I","R","I","R","I","R","I","R","I","R","I","R","I","R"),
  Cotrimoxazole=c("S","I","S","R","S","S","I","S","S","I","S","R","I","S","S","I","S","R"),
  Ceftriaxone = rep("S", 18),
  Gentamycin  = rep("S", 18),
  Doxycycline = c("S","S","I","S","S","S","S","I","S","S","S","I","S","S","I","S","S","I"),
  stringsAsFactors = FALSE
)

s_aureus <- data.frame(
  ID          = 1:11,
  Penicillin  = rep("R", 11),
  Tetracycline= c("R","R","R","R","I","R","I","R","R","I","R"),
  Cotrimoxazole=c("S","S","I","R","S","S","I","R","S","S","I"),
  Ceftriaxone = rep("S", 11),
  Gentamycin  = rep("S", 11),
  Doxycycline = c("S","I","S","I","S","S","S","I","S","S","I"),
  stringsAsFactors = FALSE
)

kec <- data.frame(
  ID          = 1:7,
  Tetracycline= c("S","R","R","I","R","I","R"),
  Cotrimoxazole=c("S","I","R","R","I","S","R"),
  Doxycycline = c("S","I","R","I","S","R","I"),
  Ceftriaxone = rep("S", 7),
  Gentamicin  = rep("S", 7),
  Florfenicol = rep("S", 7),
  stringsAsFactors = FALSE
)

e_coli <- data.frame(
  ID          = 1:8,
  Tetracycline= c("S","R","R","R","R","I","R","I"),
  Cotrimoxazole=c("S","R","R","I","R","I","R","S"),
  Doxycycline = c("S","I","S","I","S","I","R","S"),
  Ceftriaxone = rep("S", 8),
  Gentamicin  = rep("S", 8),
  Florfenicol = c("S","S","S","S","S","S","I","S"),
  stringsAsFactors = FALSE
)

pseudomonas <- data.frame(
  ID            = 1:6,
  Ciprofloxacin = c("S","R","R","S","R","S"),
  Levofloxacin  = c("S","S","R","S","R","R"),
  Gentamicin    = c("S","S","S","R","R","S"),
  Amikacin      = c("S","S","S","S","S","R"),
  Ceftriaxone   = rep("R", 6),
  Florfenicol   = rep("R", 6),
  stringsAsFactors = FALSE
)


# ================================================================
# SECTION 2: HELPER FUNCTIONS
# ================================================================

# Convert SIR to ordered factor
to_sir <- function(x) factor(toupper(x), levels = c("S","I","R"), ordered = TRUE)

# SIR to numeric (S=0, I=1, R=2) for Kruskal-Wallis
sir_numeric <- function(x) {
  as.integer(to_sir(x)) - 1
}

# Resistance summary table for one organism
resistance_summary <- function(df, organism) {
  ab_cols <- setdiff(colnames(df), "ID")
  n <- nrow(df)
  cat("\n", strrep("=", 65), "\n", sep = "")
  cat(" Organism:", organism, "| n =", n, "\n")
  cat(strrep("=", 65), "\n", sep = "")
  cat(sprintf("%-18s %6s %7s  %6s %7s  %6s %7s  %8s\n",
              "Antibiotic", "S", "S%", "I", "I%", "R", "R%", "NonSusc%"))
  cat(strrep("-", 65), "\n")
  for (ab in ab_cols) {
    vals <- toupper(df[[ab]])
    s <- sum(vals == "S"); i <- sum(vals == "I"); r <- sum(vals == "R")
    cat(sprintf("%-18s %6d %6.1f%%  %6d %6.1f%%  %6d %6.1f%%  %7.1f%%\n",
                ab, s, s/n*100, i, i/n*100, r, r/n*100, (i+r)/n*100))
  }
}

# MDR check (resistant to >= 3 drugs)
mdr_check <- function(df, organism) {
  ab_cols <- setdiff(colnames(df), "ID")
  r_count <- rowSums(sapply(ab_cols, function(ab) toupper(df[[ab]]) == "R"))
  mdr_n <- sum(r_count >= 3)
  n <- nrow(df)
  cat(sprintf("  %-22s: %d/%d (%.1f%%) are MDR (R to ≥3 antibiotics)\n",
              organism, mdr_n, n, mdr_n/n*100))
}

# Compare two organisms across shared antibiotics
# Automatically selects Fisher's Exact or Chi-Square based on expected cell counts
compare_organisms <- function(df1, df2, name1, name2) {
  ab1 <- setdiff(colnames(df1), "ID")
  ab2 <- setdiff(colnames(df2), "ID")
  shared <- intersect(ab1, ab2)

  if (length(shared) == 0) {
    cat("  No shared antibiotics between", name1, "and", name2, "\n")
    return(invisible(NULL))
  }

  cat("\n  ┌──", name1, "vs", name2,
      paste0("(n=", nrow(df1), " vs n=", nrow(df2), ")"), "\n")
  cat(sprintf("  %-18s %-15s %-12s %-10s %s\n",
              "Antibiotic", "Test Used", "Statistic", "p-value", "Sig."))
  cat("  ", strrep("-", 62), "\n", sep = "")

  results <- list()
  for (ab in shared) {
    v1 <- factor(toupper(df1[[ab]]), levels = c("S","I","R"))
    v2 <- factor(toupper(df2[[ab]]), levels = c("S","I","R"))
    tab <- rbind(table(v1), table(v2))
    # Remove zero-sum columns
    keep <- colSums(tab) > 0
    tab  <- tab[, keep, drop = FALSE]

    if (ncol(tab) < 2) {
      cat(sprintf("  %-18s %-15s %-12s %-10s %s\n",
                  ab, "N/A", "—", "uniform", ""))
      next
    }

    expected <- outer(rowSums(tab), colSums(tab)) / sum(tab)
    use_fisher <- any(expected < 5)

    if (use_fisher) {
      res <- fisher.test(tab, simulate.p.value = TRUE, B = 10000)
      test_name <- "Fisher's Exact"
      stat_str  <- ifelse(!is.null(res$estimate),
                          paste0("OR=", round(res$estimate, 2)), "—")
      p_val <- res$p.value
    } else {
      res <- chisq.test(tab)
      test_name <- "Chi-Square"
      stat_str  <- paste0("X2=", round(res$statistic, 2))
      p_val <- res$p.value
    }

    sig <- ifelse(p_val < 0.001, "***",
           ifelse(p_val < 0.01,  "**",
           ifelse(p_val < 0.05,  "*", "ns")))
    p_str <- ifelse(p_val < 0.001, "<0.001", round(p_val, 3))
    cat(sprintf("  %-18s %-15s %-12s %-10s %s\n",
                ab, test_name, stat_str, p_str, sig))

    results[[ab]] <- list(test = test_name, p = p_val, sig = sig)
  }
  invisible(results)
}

# Kruskal-Wallis test across antibiotics within one organism
kruskal_intra <- function(df, organism) {
  ab_cols <- setdiff(colnames(df), "ID")
  groups  <- lapply(ab_cols, function(ab) sir_numeric(df[[ab]]))
  # Only keep non-uniform groups
  vary    <- sapply(groups, function(g) length(unique(g)) > 1)
  groups  <- groups[vary]
  labels  <- ab_cols[vary]
  if (length(groups) < 2) {
    cat(sprintf("  %-22s: Insufficient variation for test\n", organism))
    return(invisible(NULL))
  }
  res <- kruskal.test(groups)
  sig <- ifelse(res$p.value < 0.001, "***",
         ifelse(res$p.value < 0.01,  "**",
         ifelse(res$p.value < 0.05,  "*", "ns")))
  p_str <- ifelse(res$p.value < 0.001, "<0.001",
                  formatC(res$p.value, format = "f", digits = 4))
  cat(sprintf("  %-22s: H = %.2f, df = %d, p = %s %s\n",
              organism, res$statistic, res$parameter, p_str, sig))
  cat(sprintf("  %s Tested antibiotics: %s\n\n",
              strrep(" ", 24), paste(labels, collapse = ", ")))
}


# ================================================================
# SECTION 3: RUN ANALYSIS
# ================================================================

cat("\n\n")
cat(strrep("#", 70), "\n")
cat("#       ANTIMICROBIAL SUSCEPTIBILITY TEST — STATISTICAL ANALYSIS\n")
cat(strrep("#", 70), "\n")

## --- 3A. Resistance Profiles ---
cat("\n\n▌ SECTION A: RESISTANCE PROFILES PER ORGANISM\n")
resistance_summary(s_uberis,    "S. uberis")
resistance_summary(s_aureus,    "S. aureus")
resistance_summary(kec,         "K.E.C")
resistance_summary(e_coli,      "E. coli")
resistance_summary(pseudomonas, "Pseudomonas spp.")

## --- 3B. MDR ---
cat("\n\n▌ SECTION B: MULTI-DRUG RESISTANCE (MDR) ASSESSMENT\n")
cat("  Definition: Resistant (R) to ≥ 3 antibiotics\n\n")
mdr_check(s_uberis,    "S. uberis")
mdr_check(s_aureus,    "S. aureus")
mdr_check(kec,         "K.E.C")
mdr_check(e_coli,      "E. coli")
mdr_check(pseudomonas, "Pseudomonas spp.")

## --- 3C. Inter-organism comparisons ---
cat("\n\n▌ SECTION C: INTER-ORGANISM COMPARISONS\n")
cat("  Test selection: Fisher's Exact if any expected cell < 5; else Chi-Square\n")
compare_organisms(s_uberis, s_aureus, "S. uberis", "S. aureus")
compare_organisms(kec,      e_coli,   "K.E.C",     "E. coli")

## --- 3D. Intra-organism Kruskal-Wallis ---
cat("\n\n▌ SECTION D: INTRA-ORGANISM — ANTIBIOTIC COMPARISON (Kruskal-Wallis)\n")
cat("  Ordinal scoring: S=0, I=1, R=2 | Tests resistance distribution across antibiotics\n\n")
kruskal_intra(s_uberis,    "S. uberis")
kruskal_intra(s_aureus,    "S. aureus")
kruskal_intra(kec,         "K.E.C")
kruskal_intra(e_coli,      "E. coli")
kruskal_intra(pseudomonas, "Pseudomonas spp.")


# ================================================================
# SECTION 4: VISUALISATION (ggplot2)
# ================================================================

# Function to build stacked bar chart data for one organism
make_plot_data <- function(df, organism) {
  ab_cols <- setdiff(colnames(df), "ID")
  n <- nrow(df)
  do.call(rbind, lapply(ab_cols, function(ab) {
    vals <- toupper(df[[ab]])
    data.frame(
      Organism   = organism,
      Antibiotic = ab,
      Category   = c("S","I","R"),
      Percentage = c(sum(vals=="S"), sum(vals=="I"), sum(vals=="R")) / n * 100,
      stringsAsFactors = FALSE
    )
  }))
}

plot_data <- bind_rows(
  make_plot_data(s_uberis,    "S. uberis"),
  make_plot_data(s_aureus,    "S. aureus"),
  make_plot_data(kec,         "K.E.C"),
  make_plot_data(e_coli,      "E. coli"),
  make_plot_data(pseudomonas, "Pseudomonas")
)

plot_data$Category <- factor(plot_data$Category, levels = c("R","I","S"))

sir_colors <- c("S" = "#2ecc71", "I" = "#f39c12", "R" = "#e74c3c")

p <- ggplot(plot_data, aes(x = Antibiotic, y = Percentage, fill = Category)) +
  geom_bar(stat = "identity", position = "stack", color = "white", linewidth = 0.3) +
  facet_wrap(~Organism, scales = "free_x", ncol = 2) +
  scale_fill_manual(values = sir_colors,
                    labels = c("S" = "Sensitive", "I" = "Intermediate", "R" = "Resistant"),
                    name = "Susceptibility") +
  labs(
    title    = "Antimicrobial Susceptibility Profiles by Organism",
    subtitle = "S = Sensitive | I = Intermediate | R = Resistant",
    x        = "Antibiotic",
    y        = "Percentage of Isolates (%)"
  ) +
  theme_bw(base_size = 11) +
  theme(
    axis.text.x     = element_text(angle = 40, hjust = 1, size = 9),
    strip.background= element_rect(fill = "#2c3e50"),
    strip.text      = element_text(color = "white", face = "bold"),
    legend.position = "bottom",
    plot.title      = element_text(face = "bold", size = 13),
    plot.subtitle   = element_text(size = 10, color = "gray40")
  )

print(p)
# Save plot:
# ggsave("AST_Susceptibility_Plot.png", p, width = 12, height = 9, dpi = 300)

cat("\n", strrep("=", 70), "\n")
cat("Significance codes: *** p<0.001 | ** p<0.01 | * p<0.05 | ns not significant\n")
cat(strrep("=", 70), "\n")
