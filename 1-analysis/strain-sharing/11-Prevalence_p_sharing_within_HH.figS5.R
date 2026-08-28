# This script is to investigate the probability of strain-sharing within households with members positive for strains. 
source(file.path(here::here(),"0-config.R"))
# Load data
R21_sample_list_detection = read.csv(file.path(here::here(), "data/R21_sample_list_detection.csv"))
load(file.path(here::here(), "data/HH.rates.RData"))

# genus | number of HH with detection | num of HH with strain-sharing | percentage | upper CI? | lower CI
get_prob_table <- function(sample_df, sharing_df, target_detect, area_name, target_label) {

  # Initialize output table
  prob_table <- data.frame(
    target = character(),
    area = character(),
    num_hh_detect = numeric(),
    numb_hh_sharing = numeric(),
    percentage = numeric(),
    upper = numeric(),
    lower = numeric(),
    stringsAsFactors = FALSE
  )

  # Get list of households with at least one detection
  hh.list <- sample_df %>%
    filter(.data[[target_detect]] == "Yes", area == area_name) %>%
    pull(hhid) %>%
    unique()

  # Get list of households with strain-sharing
  hh.sharing.list <- sharing_df %>%
    filter(HH_type == "Within", sharing.pair > 0, area == area_name) %>%
    pull(HH1) %>%
    unique()

  # Compute confidence interval
  ci <- binom.confint(length(hh.sharing.list), length(hh.list), method = "wilson")

  # Fill in results
  prob_table[1, ] <- c(
    target_label,
    area_name,
    length(hh.list),
    length(hh.sharing.list),
    round(length(hh.sharing.list) / length(hh.list) * 100, 2),
    round(ci$upper * 100, 2),
    round(ci$lower * 100, 2)
  )

  return(prob_table)
}
#

# Define target info as a single table
library(tibble)

target_info <- tibble(
  label = c("Escherichia", "Campylobacter", "Klebsiella", "Enterococcus",
            "Enterobacter", "Bacteroides", "Bifidobacterium"),
  detect_col = c("ecoli_detect", "camp_detect", "kleb_detect", "enterco_detect",
                 "enterob_detect", "bacteroides_detect", "bifidobacterium_detect"),
  sharing_df = list(
    hh.comb.ML.Esch, hh.comb.ML.Camp.mag, hh.comb.ML.Kleb, hh.comb.ML.Entc,
    hh.comb.ML.Entb, hh.comb.ML.Bact, hh.comb.ML.Bifi
  )
)

# Helper to run for one area
get_area_prob_table <- function(area_name) {
  do.call(rbind, lapply(1:nrow(target_info), function(i) {
    get_prob_table(
      sample_df   = R21_sample_list_detection,
      sharing_df  = target_info$sharing_df[[i]],
      target_detect = target_info$detect_col[i],
      area_name   = area_name,
      target_label = target_info$label[i]
    )
  }))
}

# Generate tables for both areas
prob_table_rural <- get_area_prob_table("rural")
prob_table_urban <- get_area_prob_table("urban")

# Combine all results
prob_table_all <- rbind(prob_table_rural, prob_table_urban)

prob_table_num_cols <- c(
  "num_hh_detect",
  "numb_hh_sharing",
  "percentage",
  "upper",
  "lower"
)

prob_table_all[, prob_table_num_cols] <- lapply(prob_table_all[, prob_table_num_cols], function(col) {
  as.numeric(as.character(unlist(col)))
})

# Plot
prob_table_all$target = factor(prob_table_all$target, levels = c("Escherichia", "Campylobacter", "Klebsiella", "Enterococcus",
                                          "Enterobacter", "Bacteroides", "Bifidobacterium"))



pd <- position_dodge(width = 0.7)  # one dodge to rule them all

prob.within.all.fig <- ggplot(prob_table_all, aes(x = target, y = percentage, fill = area)) +
  geom_col(position = pd, width = 0.7) +                       # geom_col = stat="identity"
  geom_errorbar(
    aes(ymin = lower, ymax = upper, group = area),             # ensure same grouping
    position = pd,
    width = 0.2                                                # cap width only
  ) +
  scale_fill_brewer(palette = "Set2") +
  labs(
    title = "Probability of strain-sharing within households with members positive for strains.",
    x = NULL,
    y = "Probability of strain-sharing (%)",
    fill = "study"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 16),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
    axis.text.y = element_text(size = 12),
    strip.text = element_text(face = "bold", size = 14),
    legend.title = element_text(face = "bold"),
    legend.position = "top"
  )

#Save into a file
ggsave(file.path(here::here(),"data/figures/prob.within.all.fig.pdf"), dpi = 300, scale = 0.3, width = 500, height = 500, units = "mm", plot = prob.within.all.fig)




