# Library, data loading and functions ------------------------------------------------
suppressPackageStartupMessages({
  library(tidyverse) # tidying data
  library(readxl) # loading data
  library(rstatix) # to remove outliers
  library(patchwork) # patchwork plots
  library(ggpubr) # to run compare_means in t-tests
  library(oaColors) # nice color palettes
})

## Reading data ####
setwd("~/Documents/PhD/Papers/Paper I/Manuscript/submission")
df <-
  read_excel('metabolic_genes-pfaffl.xlsx', sheet = 'Pfaffl_summary_long')

df <-
  # correcting gene names to HUGO nomenclature
  df %>%
  mutate_at('gene', str_replace, 'C-myc', 'MYC') %>%
  mutate_at('gene', str_replace, 'ATG/ULK', 'ULK1') %>%
  mutate_at('gene', str_replace, 'Cathepsin b', 'CATB') %>%
  mutate_at('gene', str_replace, 'IL-1b', 'IL1B') %>%
  mutate_at('gene', str_replace, 'Glut1', 'SLC2A1') %>%
  mutate_at('gene', str_replace, 'iNOS', 'NOS2') %>%
  mutate_at('gene', str_replace, 'mTOR', 'MTOR') %>%
  mutate_at('gene', str_replace, 'H1F1A', 'HIF1A') %>%
  arrange(.$gene)

## Defining variable levels ##
df$gene <-
  factor(
    df$gene,
    levels = c(
      'ULK1',
      'MTOR',
      'HIF1A',
      'IL1B',
      'NOS2',
      'CATB',
      'SIX1',
      'MYC',
      'SLC2A1'
    )
  )
df$treatment <-
  factor(
    df$treatment,
    levels = c(
      'Control',
      'Astragalus',
      'Hyaluronic acid',
      'Imiquimod',
      'Poly I:C'
    )
  )
df$day <- factor(df$day)
###

# Tidying data ------------------------------------------------------------
## Identifying and removing outliers using rstatix ####
# Values above Q3 + 3xIQR or below Q1 - 3xIQR are considered as extreme outliers.
identifying_outliers <- df %>%
  group_by(treatment, day, gene) %>%
  identify_outliers(ratio)  # identifying outliers in 'ratio' column

identifying_outliers <-
  full_join(df, identifying_outliers)  # joining outlier information with dataframe

identifying_outliers["is.extreme"][is.na(identifying_outliers["is.extreme"])] <-
  FALSE  # converting NAs from extreme column into FALSE

outliers_removed <-
  subset(identifying_outliers,
         identifying_outliers$is.extreme == FALSE)  # subsetting all non-extreme outliers/values into a new dataframe
###

## Testing normality ####
outliers_removed <-
  select(outliers_removed, -is.outlier, -is.extreme)

normality_test <- outliers_removed %>%
  nest_by(gene, treatment, day, .key = 'nested') %>%
  mutate(sw = shapiro.test(nested$ratio)$p.value) %>%
  mutate(normal =
           case_when(sw >= .05 ~ 'Yes',
                     sw <= .05 ~ 'No')) %>%
  print(n = Inf)
###

## Summarising data ####
no_outliers_summary <- summarise(
  # summarising data
  group_by(outliers_removed, gene, day, treatment),
  avg = mean(ratio, na.rm = T),
  stdev = sd(ratio, na.rm = T),
  se = sd(ratio, na.rm = T) /
    sqrt(length(ratio))
)
no_outliers_summary %>%
  filter(., treatment != 'Control') %>%
  arrange(., treatment, day) -> plotting_summary  # removing Control data for plotting
###


# Statistics and plotting -------------------------------------------------
stat.test <-
  compare_means(
    ratio ~ treatment,
    data = outliers_removed,
    group.by = c('gene', 'day'),
    method = 't.test',
    ref.group = 'Control',
    p.adjust.method = 'fdr',
  )

plotting_summary <- plotting_summary %>% rename(group2 = treatment)

statistics_summary <-
  inner_join(plotting_summary, stat.test)

statistics_summary <-
  as_tibble(statistics_summary) %>% select(group1, group2, gene, day, avg, se, p.signif)

statistics_summary <-
  statistics_summary %>% rename(control = group1, treatment = group2)

statistics_summary <-
  statistics_summary %>% mutate(
    single_asterisk = case_when(
      p.signif == '*' |
        p.signif == '**' |
        p.signif == '***' |
        p.signif == '****' ~ '*'
    )
  )

palette <-
  oaPalette(numColors = 5, alpha = 0.9) # creating palette for color reference

# setting limits for errorbars
limits <- aes(ymax = diff + (diff > 0) * se,
              ymin = diff - (diff < 0) * se)


## Astragalus ####
no_outliers_summary %>%
  filter((treatment == 'Astragalus') | (treatment == 'Control')) %>%
  group_by(day) %>%
  arrange(gene) %>%
  mutate(diff = avg - lag(avg, default = first(avg))) %>%
  filter(treatment != 'Control') -> plotting_summary

statistics_summary %>% 
  select(treatment, gene, day, avg, single_asterisk) %>% 
  inner_join(plotting_summary, .) -> astragalus_plots

astragalus_plots %>% 
  mutate(y = ifelse(diff < 0, diff - se - 0.4, diff +
                                         se + 0.2)) -> astragalus_plots

p5 <- astragalus_plots %>%
  ggplot(aes(x = day, y = diff, fill = gene)) +
  geom_bar(
    stat = 'identity',
    position = 'dodge2',
    colour = 'black',
    size = .2
  ) +
  geom_errorbar(limits,
                width = 0,
                position = position_dodge(.9)) +
  theme_linedraw(base_size = 10) +
  scale_y_continuous(limits = c(-2, 3),
                     breaks = scales::pretty_breaks(n = 5)) +
  xlab('Days post-immunostimulant bath') +
  ylab('Relative expression') +
  geom_vline(
    xintercept = c(1.5, 2.5, 3.5, 4.5, 5.5, 6.5),
    colour = 'black',
    linetype = 'dotted',
    lwd = 0.2
  ) +
  geom_hline(yintercept = 0,
             colour = 'grey',
             lwd = 0.2) +
  geom_text(
    data = astragalus_plots %>% filter(treatment == 'Astragalus'),
    aes(x = day,
        y = y,
        label = single_asterisk),
    position = position_dodge(width = 0.9),
    na.rm = T,
    size = 4
  ) +
  scale_fill_manual(
    values = c(
      "#FF0000",
      "#EB4E32",
      "#5361FF",
      "#3E6EFF",
      "#54C0FF",
      "#00950080",
      "#F0EE76",
      "#EAFF03",
      "#FF00FF80"
    )
  ) +
  ggtitle('Astragalus') +
  theme(
    legend.position = 'bottom',
    legend.title = element_blank(),
    legend.direction = 'horizontal',
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    # axis.ticks.x = element_blank(),
    axis.title.x = element_blank(),
  )


## Hyaluronic acid ####
no_outliers_summary %>%
  filter((treatment == 'Hyaluronic acid') |
           (treatment == 'Control')) %>%
  group_by(day) %>%
  arrange(gene) %>%
  mutate(diff = avg - lag(avg, default = first(avg))) %>%
  filter(treatment != 'Control') -> plotting_summary

statistics_summary %>% 
  select(treatment, gene, day, avg, single_asterisk) %>% 
  inner_join(plotting_summary, .) -> hyaluronica_plots

hyaluronica_plots %>% 
  mutate(y = ifelse(diff < 0, diff - se - 0.6, diff +
                                          se + 0.2)) -> hyaluronica_plots

p6 <- hyaluronica_plots %>%
  ggplot(aes(x = day, y = diff, fill = gene)) +
  geom_bar(
    stat = 'identity',
    position = 'dodge2',
    colour = 'black',
    size = .2
  ) +
  geom_errorbar(limits,
                width = 0,
                position = position_dodge(.9)) +
  theme_linedraw(base_size = 10) +
  scale_y_continuous(limits = c(-2, 6.5),
                     breaks = scales::pretty_breaks(n = 4)) +
  xlab('Days post-immunostimulant bath') +
  ylab('Relative expression') +
  geom_vline(
    xintercept = c(1.5, 2.5, 3.5, 4.5, 5.5, 6.5),
    colour = 'black',
    linetype = 'dotted',
    lwd = 0.2
  ) +
  geom_hline(yintercept = 0,
             colour = 'grey',
             lwd = 0.2) +
  geom_text(
    data = hyaluronica_plots %>% filter(treatment == 'Hyaluronic acid'),
    aes(x = day,
        y = y,
        label = single_asterisk),
    position = position_dodge(width = 0.9),
    na.rm = T,
    size = 4
  ) +
  scale_fill_manual(
    values = c(
      "#FF0000",
      "#EB4E32",
      "#5361FF",
      "#3E6EFF",
      "#54C0FF",
      "#00950080",
      "#F0EE76",
      "#EAFF03",
      "#FF00FF80"
    )
  ) +
  ggtitle('Hyaluronic acid') +
  theme(
    legend.position = 'bottom',
    legend.title = element_blank(),
    legend.direction = 'horizontal',
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    # axis.ticks.x = element_blank(),
    axis.title.x = element_blank(),
  )


## Imiquimod ####
no_outliers_summary %>%
  filter((treatment == 'Imiquimod') | (treatment == 'Control')) %>%
  group_by(day) %>%
  arrange(gene) %>%
  mutate(diff = avg - lag(avg, default = first(avg))) %>%
  filter(treatment != 'Control') -> plotting_summary

statistics_summary %>% 
  select(treatment, gene, day, avg, single_asterisk) %>% 
  inner_join(plotting_summary, .) -> imiquimod_plots

imiquimod_plots %>% mutate(y = ifelse(diff < 0, diff - se - 0.6, diff +
                                        se + 0.2)) -> imiquimod_plots


p7 <- imiquimod_plots %>%
  ggplot(aes(x = day, y = diff, fill = gene)) +
  geom_bar(
    stat = 'identity',
    position = 'dodge2',
    colour = 'black',
    size = .2
  ) +
  geom_errorbar(limits,
                width = 0,
                position = position_dodge(.9)) +
  theme_linedraw(base_size = 10) +
  scale_y_continuous(limits = c(-2, 7.5),
                     breaks = scales::pretty_breaks(n = 4)) +
  xlab('Days post-immunostimulant bath') +
  ylab('Relative expression') +
  geom_vline(
    xintercept = c(1.5, 2.5, 3.5, 4.5, 5.5, 6.5),
    colour = 'black',
    linetype = 'dotted',
    lwd = 0.2
  ) +
  geom_hline(yintercept = 0,
             colour = 'grey',
             lwd = 0.2) +
  geom_text(
    data = imiquimod_plots %>% filter(treatment == 'Imiquimod'),
    aes(x = day,
        y = y,
        label = single_asterisk),
    position = position_dodge(width = 0.9),
    na.rm = T,
    size = 4
  ) +
  scale_fill_manual(
    values = c(
      "#FF0000",
      "#EB4E32",
      "#5361FF",
      "#3E6EFF",
      "#54C0FF",
      "#00950080",
      "#F0EE76",
      "#EAFF03",
      "#FF00FF80"
    )
  ) +
  ggtitle('Imiquimod') +
  theme(
    legend.position = 'bottom',
    legend.title = element_blank(),
    legend.direction = 'horizontal',
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    # axis.ticks.x = element_blank(),
    axis.title.x = element_blank(),
  )


## Poly I:C ####
no_outliers_summary %>%
  filter((treatment == 'Poly I:C') | (treatment == 'Control')) %>%
  group_by(day) %>%
  arrange(gene) %>%
  mutate(diff = avg - lag(avg, default = first(avg))) %>%
  filter(treatment != 'Control') -> plotting_summary

statistics_summary %>% 
  select(treatment, gene, day, avg, single_asterisk) %>% 
  inner_join(plotting_summary, .) -> polyic_plots

polyic_plots %>% 
  mutate(y = ifelse(diff < 0, diff - se - 0.6, diff + se +
                                     0.2)) -> polyic_plots


p8 <- polyic_plots %>%
  ggplot(aes(x = day, y = diff, fill = gene)) +
  geom_bar(
    stat = 'identity',
    position = 'dodge2',
    colour = 'black',
    size = .2
  ) +
  geom_errorbar(limits,
                width = 0,
                position = position_dodge(.9)) +
  theme_linedraw(base_size = 10) +
  scale_y_continuous(limits = c(-2, 6),
                     breaks = scales::pretty_breaks(n = 4)) +
  xlab('Days post immunostimulant bath') +
  ylab('Relative expression') +
  geom_vline(
    xintercept = c(1.5, 2.5, 3.5, 4.5, 5.5, 6.5),
    colour = 'black',
    linetype = 'dotted',
    lwd = 0.2
  ) +
  geom_hline(yintercept = 0,
             colour = 'grey',
             lwd = 0.2) +
  geom_text(
    data = polyic_plots %>% filter(treatment == 'Poly I:C'),
    aes(x = day,
        y = y,
        label = single_asterisk),
    position = position_dodge(width = 0.9),
    na.rm = T,
    size = 4
  ) +
  scale_fill_manual(
    values = c(
      "#FF0000",
      "#EB4E32",
      "#5361FF",
      "#3E6EFF",
      "#54C0FF",
      "#00950080",
      "#F0EE76",
      "#EAFF03",
      "#FF00FF80"
    )
  ) +
  ggtitle('Poly I:C') +
  theme(
    legend.position = 'bottom',
    legend.title = element_blank(),
    legend.direction = 'horizontal',
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    # axis.ticks.x = element_blank(),
    # axis.title.x = element_blank(),
  )


## Patchwork ####
patchwork <- (p5 / p6 / p7 / p8)

patchwork + 
  plot_annotation(title = 'Metabolic genes',
                  theme = theme(plot.title = element_text(hjust = 0.5, vjust = 2))) +
  plot_layout(guides = 'collect') &
  theme(
    legend.position = 'bottom',
    legend.text = element_text(margin = margin(r = .5, unit = 'cm')),
    text = element_text(size = 14, family = 'serif'),
    axis.text = element_text(face = 'bold')
  )
