library(ggplot2)
library(ggpubr)
library(dplyr)
library(cowplot)

# plot cognition pattern for each LV
# maximal model
beh_plot_max <- function(df) {
    ggplot(data=CI) +
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


#plot progression of brain and cognition patterns per DX group

delta <- read.csv("NIFD_PLS_csvfiles/delta.csv", stringsAsFactors=TRUE)
col=c('#377eb8','#ff7f00','#4daf4a')
my_comparisons <- list( c("SV", "BV"), c("SV", "PNFA"), c("PNFA", "BV") )

plot_delta <- function(variable) {
  ggplot(delta, aes(x = DX, y = .data[[variable]], color = DX)) +
  geom_jitter()+
  geom_boxplot(width = 0.2,)+
  scale_colour_manual(values = col, aesthetics = c("fill", "color"))+
  labs(x="",y='')+
  theme_light()+ 
  stat_compare_means(comparisons = my_comparisons,method = "t.test",paired=FALSE,label='p.signif',hide.ns=TRUE)
}

plots <- lapply(paste0("delta_", 1:8), plot_delta)
ggarrange(plots[[1]], plots[[2]], plots[[5]], plots[[6]], 
          common.legend = TRUE, legend = 'right')
ggarrange(plots[[3]], plots[[4]], plots[[7]], plots[[8]], 
          common.legend = TRUE, legend = 'right')

#plot progression of brain/cognition patterns across all groups
compare_base <- read.csv("NIFD_PLS_csvfiles/compare_base.csv", header=TRUE)
compare_lng <- read.csv("NIFD_PLS_csvfiles/compare_lng.csv", header=TRUE)

compare_base <- compare_base |>
  mutate(time = 'baseline')
compare_lng <- compare_lng |>
  mutate(time = '1-year follow-up')
compare <- rbind(compare_base, compare_lng)

c2=c("#CC79A7","#D55E00")

create_lv_plot <- function(data, y_var, title) {
  ggplot(data = data, aes(x = time, y = .data[[y_var]], color = time)) +
    geom_boxplot() +
    geom_jitter(shape = 16, alpha = 0.9, height = 0.1) +
    labs(title = title) +
    scale_color_manual(values = c2) +
    theme_classic() +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.title.x = element_blank(),
      axis.title.y = element_blank()
    )
}

# Create plot titles
behavior_titles <- paste0("LV-", 1:4, " behavior pattern")
brain_titles <- paste0("LV-", 1:4, " brain pattern")

# Generate plots for behavior and brain patterns
behavior_plots <- lapply(1:4, function(i) create_lv_plot(compare, paste0("lng_", i), behavior_titles[i]))
brain_plots <- lapply(5:8, function(i) create_lv_plot(compare, paste0("lng_", i), brain_titles[i - 4]))

# Arrange the plots into two grids
ggarrange(plotlist = behavior_plots, ncol = 4, common.legend = TRUE)
ggarrange(plotlist = brain_plots, ncol = 4, common.legend = TRUE)
