library(tidyverse)
library(patchwork)
library(here)
library(sessioninfo)

set.seed(0)

plot_dir = here('plots', '13_tripod_trios', '17_trio_cartoon') 
n_per_panel <- 30

dir.create(plot_dir, showWarnings = FALSE)

tf_levels <- tibble(
  tf_group = factor(c("Low TF", "Mid TF", "High TF"),
    levels = c("Low TF", "Mid TF", "High TF")
  ),
  tf_value = c(0.25, 0.55, 0.85),
  slope = c(0.55, 0.9, 1.5),
  intercept = c(0.12, 0.08, 0.04)
)

cartoon_df <- tf_levels |>
  mutate(data = pmap(
    list(tf_group, tf_value, slope, intercept),
    function(tf_group, tf_value, slope, intercept) {
      peak_accessibility <- rbeta(n_per_panel, shape1 = 2.6, shape2 = 2.2)
      noise_sd <- 0.14

      tibble(
        tf_group = tf_group,
        tf_expression = tf_value,
        peak_accessibility = peak_accessibility,
        gene_expression = intercept +
          slope * peak_accessibility +
          rnorm(n_per_panel, mean = 0, sd = noise_sd)
      ) |>
        mutate(gene_expression = pmin(pmax(gene_expression, 0), 1.45))
    }
  )) |>
  select(data) |>
  unnest(data)

p <- ggplot(
  cartoon_df,
  aes(
    x = peak_accessibility,
    y = gene_expression,
    color = tf_expression
  )
) +
  geom_point(alpha = 0.75, size = 2) +
  geom_abline(
    intercept = 0,
    slope = 1,
    linetype = "dotted",
    linewidth = 0.5,
    color = "grey55"
  ) +
  geom_smooth(
    method = "lm",
    se = FALSE,
    linewidth = 0.9,
    color = "black"
  ) +
  facet_wrap(~tf_group, nrow = 1) +
  scale_color_viridis_c() +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1.45), expand = FALSE) +
  labs(
    x = "Peak accessibility",
    y = "Gene expression",
    color = "TF expression"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    strip.background = element_rect(fill = "grey95"),
    legend.position = "bottom"
  )

ggsave(
  filename = file.path(plot_dir, "trio_cartoon_panel.pdf"),
  plot = p,
  width = 5,
  height = 3
)

session_info()
