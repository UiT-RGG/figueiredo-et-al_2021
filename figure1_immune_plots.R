# Libraries, functions and data loading ------------------------------------------------
suppressPackageStartupMessages({
  library(NCmisc)
  library(grDevices)
  library(tidyverse)
  library(readxl)
  library(rstatix)
  library(patchwork)
  library(ggpubr)
})

## Pastellizing colors function ##
hsv2rgb <- function(x) {
  h <- x[1, 1]
  s <- x[2, 1]
  v <- x[3, 1]
  
  C <- s * v
  
  hdash <- h * 6
  X <- C * (1 - abs(hdash %% 2 - 1))
  
  if (0 <= hdash & hdash <= 1)
    RGB1 <- c(C, X, 0)
  if (1 <= hdash & hdash <= 2)
    RGB1 <- c(X, C, 0)
  if (2 <= hdash & hdash <= 3)
    RGB1 <- c(0, C, X)
  if (3 <= hdash & hdash <= 4)
    RGB1 <- c(0, X, C)
  if (4 <= hdash & hdash <= 5)
    RGB1 <- c(X, 0, C)
  if (5 <= hdash & hdash <= 6)
    RGB1 <- c(C, 0, X)
  
  RGB1 + (v - C)
}

pastellize <- function(x, p) {
  # x is a colour
  # p is a number in [0,1]
  # p = 1 will give no pastellization
  
  # convert hex or letter names to rgb
  if (is.character(x))
    x <- col2rgb(x) / 255
  
  # convert vector to rgb
  if (is.numeric(x))
    x <- matrix(x, nr = 3)
  
  col <- rgb2hsv(x, maxColorValue = 1)
  col[2, 1] <- col[2, 1] * p
  col <- hsv2rgb(col)
  
  # return in convenient format for plots
  rgb(col[1], col[2], col[3])
}

## Pastellizing colours ##
pastel_1 <- pastellize('#00CCFF', 0.8)
pastel_2 <- pastellize('#3399FF', 0.8)
pastel_3 <- pastellize('#0066CC', 0.8)
pastel_4 <- pastellize('#00FF00', 0.8)
pastel_5 <- pastellize('#009900', 0.8)
pastel_6 <- pastellize('#FF3300', 0.8)
pastel_7 <- pastellize('#CC3300', 0.8)
pastel_8 <- pastellize('#993399', 0.8)
pastel_9 <- pastellize('#663399', 0.8)
pastel_10 <- pastellize('#6600CC', 0.8)
###

## Reading data ####
setwd("~/Documents/PhD/Papers/Paper I/Raw data")
df <-
  read_excel('immune_genes-pfaffl.xlsx', sheet = 'Pfaffl_summary_long')

df <-
  # correcting gene names to HUGO nomenclature
  df %>%
  mutate_at('gene', str_replace, 'IFNa', 'IFN\u03B1') %>%
  mutate_at('gene', str_replace, 'MX', 'MX1') %>%
  mutate_at('gene', str_replace, 'Viperin', 'RSAD2') %>%
  arrange(.$gene)

## Defining variable levels ##
df$gene <-
  factor(
    df$gene,
    levels = c(
      'TLR3',
      'TLR7',
      'MDA5',
      'STAT1',
      'PKR',
      'IFN\u03B1',
      'IFNc',
      'ISG15',
      'MX1',
      'RSAD2'
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
  group_by(outliers_removed, gene, day, treatment),
  avg = mean(ratio, na.rm = T),
  se = sd(ratio, na.rm = T) /
    sqrt(length(ratio))
)
no_outliers_summary %>% 
  filter(., treatment != 'Control') %>% 
  arrange(., treatment, day) -> plotting_summary  # removing Control data for plotting
###

# Statistics and plotting --------------------------------------------------------------
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
###

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
  mutate(y = ifelse(diff < 0, diff - se - 0.5, diff +
                                         se + 0.2)) -> astragalus_plots

p1 <- astragalus_plots %>%
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
      pastel_1,
      pastel_2,
      pastel_3,
      pastel_4,
      pastel_5,
      pastel_6,
      pastel_7,
      pastel_8,
      pastel_9,
      pastel_10
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
plotting_summary
statistics_summary %>% select(treatment, gene, day, avg, single_asterisk) %>% inner_join(plotting_summary, .) -> hyaluronica_plots
hyaluronica_plots %>% mutate(y = ifelse(diff < 0, diff - se - 0.5, diff +
                                          se + 0.2)) -> hyaluronica_plots

p2 <- hyaluronica_plots %>%
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
      pastel_1,
      pastel_2,
      pastel_3,
      pastel_4,
      pastel_5,
      pastel_6,
      pastel_7,
      pastel_8,
      pastel_9,
      pastel_10
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

imiquimod_plots %>% 
  mutate(y = ifelse(diff < 0, diff - se - 0.6, diff +
                                        se + 0.3)) -> imiquimod_plots


p3 <- imiquimod_plots %>%
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
  scale_y_continuous(limits = c(-2, 8),
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
      pastel_1,
      pastel_2,
      pastel_3,
      pastel_4,
      pastel_5,
      pastel_6,
      pastel_7,
      pastel_8,
      pastel_9,
      pastel_10
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
  mutate(y = ifelse(diff < 0, diff - se - 0.2, diff + se +
                                     0.1)) -> polyic_plots

p4 <- polyic_plots %>%
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
  scale_y_continuous(limits = c(-1, 3),
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
      pastel_1,
      pastel_2,
      pastel_3,
      pastel_4,
      pastel_5,
      pastel_6,
      pastel_7,
      pastel_8,
      pastel_9,
      pastel_10
    )
  ) +
  ggtitle('Poly I:C') +
  theme(
    legend.position = 'bottom',
    legend.title = element_blank(),
    legend.direction = 'horizontal',
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank()
  )


## Patchwork ####
patchwork <- (p1 / p2 / p3 / p4)
patchwork + 
  plot_annotation(title = 'Immune genes',
                  theme = theme(plot.title = element_text(hjust = 0.5, vjust = 2))) +
  plot_layout(guides = 'collect') &
  theme(
    legend.position = 'bottom',
    legend.text = element_text(margin = margin(r = .8, unit = 'cm')),
    text = element_text(size = 14, family = 'Palatino'),
    axis.text = element_text(face = 'bold')
  )


# Checking package usage and references --------------------------------------------------
# getwd()
# p <-
#   NCmisc::list.functions.in.file(
#     '/Users/ffi007/Documents/PhD/Papers/Paper I/Manuscript/scripts/manuscript_immune_plots.R'
#   )
# summary(p)
# p$'package:grDevices'

## Citations ####
# if(nchar(system.file(package="tidyverse"))) citation("tidyverse") # Wickam et al., 2019
# if(nchar(system.file(package="ggpubr"))) citation("ggpubr") # Kassambara 2020
# if(nchar(system.file(package="grDevices"))) citation("grDevices") # R Core Team 2020
# if(nchar(system.file(package="NCmisc"))) citation("NCmisc") # Cooper 2018
# if(nchar(system.file(package="patchwork"))) citation("patchwork") # Pedersen 2020
# if(nchar(system.file(package="rstatix"))) citation("rstatix") # Kassambara 2021
# if(nchar(system.file(package="stats"))) citation("stats") # R Core Team 2020
# if(nchar(system.file(package="utils"))) citation("utils") # R Core Team 2020
# if(nchar(system.file(package="oaColors"))) citation("oaColors") # Waddell 2015
# if(nchar(system.file(package="readxl"))) citation("readxl") # Wickam and Bryan 2019
# if(nchar(system.file(package="Hmisc"))) citation("Hmisc") # Harrell Jr. et al., 2021
