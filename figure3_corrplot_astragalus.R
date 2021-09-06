# Library and data loading ------------------------------------------------
suppressPackageStartupMessages({
  library(tidyverse)
  library(readxl)
  library(rstatix)
  library(corrr)  # creating correlograms
  library(patchwork)
  library(Hmisc)  # for rcorr
  library(oaColors)
})

## Reading and selecting data ####
setwd("~/Documents/PhD/Papers/Paper I/Raw data")
df <-
  read_excel('correlations_immune-vs-metabolic.xlsx', sheet = 'immune vs metabolic')

df <-
  # correcting gene names to HUGO nomenclature
  df %>% filter(treatment == 'Astragalus') %>%
  mutate_at('gene', str_replace, 'HIF1a', 'HIF1A') %>%
  mutate_at('gene', str_replace, 'IFNa', 'IFN\u03B1') %>%
  mutate_at('gene', str_replace, 'C-myc', 'MYC') %>%
  mutate_at('gene', str_replace, 'ATG/ULK', 'ULK1') %>%
  mutate_at('gene', str_replace, 'Cathepsin b', 'CATB') %>%
  mutate_at('gene', str_replace, 'IL-1b', 'IL1B') %>%
  mutate_at('gene', str_replace, 'MX', 'MX1') %>%
  mutate_at('gene', str_replace, 'mTOR', 'MTOR') %>%
  mutate_at('gene', str_replace, 'iNOS', 'NOS2') %>%
  mutate_at('gene', str_replace, 'Glut1', 'SLC2A1') %>%
  mutate_at('gene', str_replace, 'Viperin', 'RSAD2') %>%
  arrange(.$gene)
###

# Tidying data ------------------------------------------------------------
## Identifying and removing outliers using rstatix ####
# Values above Q3 + 3xIQR or below Q1 - 3xIQR are considered as extreme outliers
astragalus_identifying_outliers <- df %>%
  group_by(treatment, day, gene) %>%
  identify_outliers(ratio)  # identifying outliers in 'ratio' column

astragalus_identifying_outliers <-
  full_join(df, astragalus_identifying_outliers)  # joining outlier information with dataframe

# converting NAs from extreme column into FALSE
astragalus_identifying_outliers["is.extreme"][is.na(astragalus_identifying_outliers["is.extreme"])] <-
  FALSE

# subsetting all non-extreme outliers/values into a new dataframe
astragalus_outliers_removed <-
  subset(
    astragalus_identifying_outliers,
    astragalus_identifying_outliers$is.extreme == FALSE
  )
###

## Organizing data by gene and day ####
astragalus_correlation <- astragalus_outliers_removed %>%
  group_by(gene, day) %>%
  summarise_at(vars(ratio),
               list(ratio = mean))
###

## Preparing dataframe for correlation ####
# reordering columns to immune + metabolic
col_order <-
  c(
    'TLR3',
    'TLR7',
    'MDA5',
    'STAT1',
    'PKR',
    'IFN\u03B1',
    'IFNc',
    'ISG15',
    'MX1',
    'RSAD2',
    'ULK1',
    'CATB',
    'MYC',
    'SLC2A1',
    'HIF1A',
    'IL1B',
    'NOS2',
    'MTOR',
    'SIX1'
  )

astragalus_wide <- astragalus_correlation %>%
  pivot_wider(names_from = gene,
              values_from = ratio) %>%
  column_to_rownames(., var = 'day') %>%
  .[, col_order]
###

# Correlations -------------------------------------------------------------
## Creating functions and tidying data ####
cors <-
  function(df) {
    # function that formats a dataframe for corr and turns each of the three elements of the list into separate dataframes
    M <- Hmisc::rcorr(as.matrix(df), type = 'spearman')
    Mdf <- map(M, ~ data.frame(.x))
  }


formatted_cors <-
  function(df) {
    # formatting correlations into single dataframe
    cors(df) %>%
      map( ~ select(
        .x,
        'ULK1',
        'CATB',
        'MYC',
        'SLC2A1',
        'HIF1A',
        'IL1B',
        'NOS2',
        'MTOR',
        'SIX1'
      )) %>%
      map(~ slice(.x, 1:10)) %>%
      map(~ rownames_to_column(.x, var = "immune")) %>%
      map(~ pivot_longer(.x, -immune, "metabolic")) %>%
      bind_rows(.id = "id") %>%
      pivot_wider(names_from = id, values_from = value) %>%
      mutate(
        sig_p = ifelse(P < .05, T, F),
        p_if_sig = ifelse(P < .05, P, NA),
        r_if_sig = ifelse(P < .05, r, NA),
        above_80 = ifelse(r > 0.8, TRUE, NA)
      )
  }

astragalus_formatted <- formatted_cors(astragalus_wide)

# Plotting ####
P1 <- astragalus_formatted %>%
  ggplot(aes(immune, metabolic, col = r)) +
  geom_tile(col = "black", fill = "white") +
  geom_point(aes(size = abs(r)), shape = 15) +
  geom_text(
    aes(
      x = immune,
      y = metabolic,
      label = ifelse(sig_p == 'TRUE', '*', '')
    ),
    vjust = 0.8,
    size = 6,
    colour = 'black'
  ) +
  theme_linedraw() +
  scale_color_gradient2(
    mid = "white",
    low = "#009500FF",
    high = "#FF0000FF",
    limits = c(-1, 1)
  ) +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  scale_size(range = c(0, 11), guide = NULL) +
  theme(
    panel.border = element_rect(color = 'black', fill = NA),
    legend.position = 'bottom',
    legend.direction = 'horizontal',
    legend.title = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title = element_blank()
  ) + ggtitle('Astragalus')
