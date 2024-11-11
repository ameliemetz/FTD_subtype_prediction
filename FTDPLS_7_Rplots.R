# this code creates a bunch of figures
# requires PLS results / longitudinal patterns as input (see matlab code)
# 1, PLS results: cognition patterns (which cognitive scores contribute how much to PLS latent variables)
# 2, Linear regression results: differences in individual brain/cognition scores (based on PLS latent variables) between FTD subtypes
# 3, longitudinal progression of brain/cognition scores (based on PLS latent variables) of individual participants per FTD subgroup
# 4, longitudinal progression of brain/cognition scores (based on PLS latent variables) of individual participants in all participants

library(ggplot2)
library(ggpubr)
library(dplyr)
library(tidyr)
library(cowplot)
setwd("/data/dadmah/metame/NIFD_PLSpaper")

#######################################################################################################
# 1
# plot cognition pattern for each LV
# maximal model
beh_plot_max <- function(df) {
    ggplot(data=df) +
    geom_bar(aes(y = labels, x = -1*mean, fill=sig), stat = 'identity') +
    xlim(-1.2,1.2)+
    geom_errorbar(aes(y = labels, xmin = -1*lower, xmax = -1*upper), width = 0.2,
                  position = position_dodge(width = 0.9)) +
    theme_minimal()+
    scale_fill_manual(values=c('grey','darkred','darkblue'))+
    scale_y_discrete(limits = rev(levels(df$labels)))+
    labs(x='',y='')+
    theme(legend.position = "none",axis.text.y = element_text(size = 15))
}

ci_files <- paste0("NIFD_PLS_csvfiles/CI_beh", 1:4, ".csv")
ci_data <- lapply(ci_files, function(file) {
  df <- read.csv(file, stringsAsFactors = TRUE)
  df$labels <- factor(df$labels, levels = unique(df$labels))
  df$sig <- as.factor(df$sig)
  return(df)
})

plots <- lapply(ci_data, beh_plot_max)

plot_grid(plotlist = plots, nrow = 4)

# minimal model
beh_plot_min <- function(df) {
    ggplot(data = df) +
    geom_bar(aes(y = labels, x = mean, fill=sig), stat = 'identity') +
    xlim(-1.2,1.2)+
    geom_errorbar(aes(y = labels, xmin = lower, xmax = upper), width = 0.2,
                  position = position_dodge(width = 0.9)) +
    theme_minimal()+
    scale_fill_manual(values=c('grey','darkblue','darkred'))+
    scale_y_discrete(limits = rev(levels(df$labels)))+
    labs(x='',y='')+
    theme(legend.position = "none",axis.text.y = element_text(size = 15))
}

ci_files <- paste0("NIFD_PLS_csvfiles/CI_beh", 1:3, "m.csv")
ci_data <- lapply(ci_files, function(file) {
  df <- read.csv(file, stringsAsFactors = TRUE)
  df$labels <- factor(df$labels, levels = unique(df$labels))
  df$sig <- as.factor(df$sig)
  return(df)
})

plots <- lapply(ci_data, beh_plot_min)
plot_grid(plotlist = plots, nrow = 3)


#######################################################################################################
# 2
# plot progression of brain and cognition patterns per DX group

file_paths <- list("NIFD_PLS_csvfiles/scores_LVI.csv", "NIFD_PLS_csvfiles/scores_LVII.csv")
#file_paths <- list("NIFD_PLS_csvfiles/scores_LVIII.csv", "NIFD_PLS_csvfiles/scores_LVIV.csv")

LV_list <- lapply(file_paths, function(path) {
  read.csv(path) %>%
    mutate(DX = case_when(
      dx == "BV" ~ "bvFTD",
      dx == "SV" ~ "svPPA",
      dx == "PNFA" ~ "nfvPPA",
      TRUE ~ "control"
    ))
})

col <- c("bvFTD" = "#1982C4", "svPPA" = "#8AC926", "nfvPPA" = "#FF595E", "control" = "#6A4C93")

# Function to calculate adjusted R^2 and p-value
get_model_stats <- function(data) {
  model <- lm(usc ~ vsc * DX, data = data)
  adj_r2 <- summary(model)$adj.r.squared
  p_val <- summary(aov(model))[[1]]$`Pr(>F)`[1]
  list(adj_r2 = adj_r2, p_val = p_val)
}

# Function to create PLS plot for a given latent variable
create_pls_plot <- function(data, title) {
  stats <- get_model_stats(data)
  label <- sprintf("adj.R2 = %.3f, p < 0.001", stats$adj_r2)
  
  ggplot(data, aes(x = vsc, y = usc, group = DX, color = DX, fill = DX, linetype = DX)) +
    geom_point(aes(shape = DX), size=2) +
    geom_smooth(method = "lm", formula = y ~ x, alpha = 0.2, fullrange = TRUE) +
    xlim(-7, 7) + ylim(-7, 7) +
    scale_color_manual(values = col) +
    scale_fill_manual(values = col) +
    scale_linetype_manual(values = c(3,2,1)) +
    labs(title = title, x = "Brain Score", y = "Cognitive Score", color = "Diagnosis", fill = "Diagnosis",shape = "Diagnosis", linetype = "Diagnosis") +
    annotate("text", x = 3, y = -5, label = label, size = 5) +
    theme_minimal() +
    guides(color = guide_legend(override.aes = list(shape = 16, linetype = 1))) +
    theme(
      plot.title = element_text(size = 16, face = "bold"),    # Title font size
      axis.title = element_text(size = 16),                   # Axis labels font size
      axis.text = element_text(size = 14),                    # Axis tick labels font size
      legend.title = element_text(size = 14),                 # Legend title font size
      legend.text = element_text(size = 14)                   # Legend text font size
    ) 
}



# Generate plots for Latent Variables I and II
plots <- mapply(create_pls_plot, LV_list, c("Latent Variable I", "Latent Variable II"), SIMPLIFY = FALSE)
#plots <- mapply(create_pls_plot, LV_list, c("Latent Variable III", "Latent Variable IV"), SIMPLIFY = FALSE)

ggarrange(plotlist = plots, common.legend = TRUE, legend = "right")


#######################################################################################################
# 3
# plot progression of brain and cognition patterns per FTD subgroup

col <- c("bvFTD" = "#1982C4", "svPPA" = "#8AC926", "nfvPPA" = "#FF595E", "control" = "#6A4C93")
my_comparisons <- list(c("svPPA", "bvFTD"), c("svPPA", "nfvPPA"), c("nfvPPA", "bvFTD"))

delta <- delta %>%
  mutate(DX = case_when(
    DX == "BV" ~ "bvFTD",
    DX == "SV" ~ "svPPA",
    DX == "PNFA" ~ "nfvPPA",
    TRUE ~ "control"
  ))

# Function to create boxplot with jitter and comparisons
create_boxplot <- function(data, y_value, title, y_label) {
  ggplot(data, aes(x = DX, y = !!sym(y_value), color = DX)) +
    geom_jitter(size = 2) +  # Adjust jitter point size
    geom_boxplot(aes(fill = DX, col = DX), width = 0.2, alpha = 0.2, outlier.shape = NA) +
    scale_colour_manual(values = col, guide = "none") +
    scale_fill_manual(values = col, guide = "none") +
    labs(title = title, x = "", y = y_label) +
    ylim(-3, 4) +
    theme_minimal() +
    stat_compare_means(comparisons = my_comparisons, method = "t.test", paired = FALSE, label = 'p.signif', hide.ns = TRUE) +
    theme(
      plot.title = element_text(size = 16, face = "bold"),    # Title font size
      axis.title = element_text(size = 14),                   # Axis labels font size
      axis.text = element_text(size = 12),                    # Axis tick labels font size
      legend.title = element_text(size = 12),                 # Legend title font size
      legend.text = element_text(size = 10)                   # Legend text font size
    )
}

# Create individual plots LV I and II
a <- create_boxplot(delta, "delta_1", "Yearly rate of change in cognitive scores
Latent variable I", 'delta')
b <- create_boxplot(delta, "delta_2", " 
Latent variable II", 'delta')
c <- create_boxplot(delta, "delta_5", "Yearly rate of change in brain scores
Latent variable I", 'delta')
d <- create_boxplot(delta, "delta_6", " 
Latent variable II", 'delta')

ggarrange(a, b, c, d, common.legend = TRUE, legend = "none")

# Create individual plots LV III and IV
a <- create_boxplot(delta, "delta_3", "Yearly rate of change in cognitive scores
Latent variable III", 'delta')
b <- create_boxplot(delta, "delta_4", " 
Latent variable IV", 'delta')
c <- create_boxplot(delta, "delta_7", "Yearly rate of change in brain scores
Latent variable III", 'delta')
d <- create_boxplot(delta, "delta_8", " 
Latent variable IV", 'delta')

ggarrange(a, b, c, d, common.legend = TRUE, legend = "none")


#######################################################################################################
# 4
# plot progression of brain and cognition patterns in all participants

compare_base <- read.csv("NIFD_PLS_csvfiles/compare_base.csv", header=TRUE)
compare_lng <- read.csv("NIFD_PLS_csvfiles/compare_lng.csv", header=TRUE)

compare_base <- compare_base |>
  mutate(time = 'baseline')
compare_lng <- compare_lng |>
  mutate(time = '1-year follow-up')
compare <- rbind(compare_base, compare_lng)

my_comparisons <- list(c('baseline','1-year follow-up'))
c2=c("#CC79A7","#D55E00")

create_lv_plot <- function(data, y_var, title, y_title) {
  ggplot(data = data, aes(x = factor(time, level=c('baseline', '1-year follow-up')), y = .data[[y_var]], color = time)) +
    geom_boxplot() +
    geom_jitter(shape = 16, alpha = 0.9, height = 0.1) +
    labs(title = title) +
    scale_color_manual(values = c2, guide = "none") +
    labs(title = title, x = "", y = y_title) +
    ylim(-6, 8) +
    theme_minimal() +
    stat_compare_means(comparisons = my_comparisons, method = "t.test", paired = TRUE, label = 'p.signif', hide.ns = TRUE) +
    theme(
      plot.title = element_text(size = 16, face = "bold"),    # Title font size
      axis.title = element_text(size = 14),                   # Axis labels font size
      axis.text = element_text(size = 12),                    # Axis tick labels font size
    )
}

behavior_titles <- paste0("Latent variable ", c("I","II","III","IV"))
brain_titles <- paste0("Latent variable ", c("I","II","III","IV"))
behavior_y <- paste0("behavior scores")
brain_y <- paste0("brain scores")

behavior_plots <- lapply(1:4, function(i) create_lv_plot(compare, paste0("lng_", i), behavior_titles[i], behavior_y))
brain_plots <- lapply(5:8, function(i) create_lv_plot(compare, paste0("lng_", i), brain_titles[i - 4], brain_y))

ggarrange(plotlist = behavior_plots, ncol = 4,common.legend = TRUE)
ggarrange(plotlist = brain_plots, ncol = 4, common.legend = TRUE)
