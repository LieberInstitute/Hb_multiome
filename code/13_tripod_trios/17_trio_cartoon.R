library(tidyverse)
library(patchwork)
library(here)
library(sessioninfo)

set.seed(0)

plot_dir = here('plots', '13_tripod_trios', '17_trio_cartoon') 
n_per_panel <- 20

dir.create(plot_dir, showWarnings = FALSE)

tf_levels <- tibble(
  tf_group = factor(c("Low TF", "Mid TF", "High TF"),
    levels = c("Low TF", "Mid TF", "High TF")
  ),
  tf_value = c(0.2, 0.6, 1),
  slope = c(0.5, 1.75, 3),
  intercept = c(0.12, 0.08, 0.04)
)

peak_accessibility_shared <- rbeta(n_per_panel, shape1 = 2.6, shape2 = 2.2)

cartoon_df <- tf_levels |>
  mutate(data = pmap(
    list(tf_group, tf_value, slope, intercept),
    function(tf_group, tf_value, slope, intercept) {
      noise_sd <- 0.14

      tibble(
        tf_group = tf_group,
        tf_expression = tf_value,
        peak_accessibility = peak_accessibility_shared,
        gene_expression = intercept +
          slope * peak_accessibility_shared +
          rnorm(n_per_panel, mean = 0, sd = noise_sd)
      ) |>
        mutate(gene_expression = pmax(gene_expression, 0))
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
  geom_point(size = 3) +
  geom_smooth(
    aes(group = tf_group),
    method = "lm",
    se = FALSE,
    linewidth = 0.9
  ) +
  scale_color_viridis_c() +
  labs(
    x = "Peak accessibility", y = "Gene expression", color = "TF expression"
  ) +
  theme_minimal(base_size = 18) +
  theme(
    legend.position = "bottom"
  )

ggsave(
  filename = file.path(plot_dir, "trio_cartoon_single_panel.pdf"), plot = p
)

session_info()
