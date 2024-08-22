### Load Packages ----
library(tidyverse)
library(formattable)
library(fmsb)

### Load Data ----
result_data <- read_csv("data/result.csv")
load("data/metrics.rds")

### Result Table and Graph ----
# alphabetically rearrange rows by tool name and round to 3 decimal places
table_metrics <- metrics %>%
  arrange(Tool) %>%
  mutate(across(where(is.numeric), round, 3))

# result table
# color by type
formattable(table_metrics,
            align = c("l","c","c","c","c"),
            list(Accuracy = color_tile("white", "#28c752"),
                 Sensitivity = color_tile("white", "#28c752"),
                 Specificity = color_tile("white", "#28c752"),
                 PPV = color_tile("white", "#28c752"),
                 NPV = color_tile("white", "#28c752"),
                 MCC = color_tile("white", "#28c752")))

# standardized color (0 to 1)
formattable(table_metrics,
            align = c("l","c","c","c","c"),
            list(Accuracy = formatter("span", 
                                      style = function(x){
                                        style(display            = "block",
                                              padding            = "0 4px",
                                              `border-radius`    = "4px",
                                              `background-color` = csscolor(gradient(
                                                as.numeric(c(c(0,1), unlist(table_metrics[[2]]))),
                                                "white", "#28c752"))[-(1:2)]
                                        )}),
                 Sensitivity = formatter("span", 
                                         style = function(x){
                                           style(display            = "block",
                                                 padding            = "0 4px",
                                                 `border-radius`    = "4px",
                                                 `background-color` = csscolor(gradient(
                                                   as.numeric(c(c(0,1), unlist(table_metrics[[3]]))),
                                                   "white", "#28c752"))[-(1:2)]
                                           )}),
                 Specificity = formatter("span", 
                                         style = function(x){
                                           style(display            = "block",
                                                 padding            = "0 4px",
                                                 `border-radius`    = "4px",
                                                 `background-color` = csscolor(gradient(
                                                   as.numeric(c(c(0,1), unlist(table_metrics[[4]]))),
                                                   "white", "#28c752"))[-(1:2)]
                                           )}),
                 PPV = formatter("span", 
                                 style = function(x){
                                   style(display            = "block",
                                         padding            = "0 4px",
                                         `border-radius`    = "4px",
                                         `background-color` = csscolor(gradient(
                                           as.numeric(c(c(0,1), unlist(table_metrics[[5]]))),
                                           "white", "#28c752"))[-(1:2)]
                                   )}),
                 NPV = formatter("span", 
                                 style = function(x){
                                   style(display            = "block",
                                         padding            = "0 4px",
                                         `border-radius`    = "4px",
                                         `background-color` = csscolor(gradient(
                                           as.numeric(c(c(0,1), unlist(table_metrics[[6]]))),
                                           "white", "#28c752"))[-(1:2)]
                                   )}),
                 MCC = formatter("span", 
                                 style = function(x){
                                   style(display            = "block",
                                         padding            = "0 4px",
                                         `border-radius`    = "4px",
                                         `background-color` = csscolor(gradient(
                                           as.numeric(c(c(0,1), unlist(table_metrics[[7]]))),
                                           "white", "#28c752"))[-(1:2)]
                                   )})))

# tidy table_metrics tibble for visualization
metrics_tidy <- table_metrics %>%
  pivot_longer(cols = c(Accuracy, Sensitivity, Specificity, PPV, NPV, MCC),
               names_to = "Performance Metric",
               values_to = "Value")

# result barplot faceted by performance metric
ggplot(data = metrics_tidy, mapping = aes(x = Tool, 
                                          y = Value)) +
  geom_col() +
  facet_wrap(~`Performance Metric`,
             scales = "free_y")

# result barplot colored by tool type
ggplot(data = metrics_tidy, mapping = aes(x = `Performance Metric`, 
                                          y = Value)) +
  geom_col(aes(fill = Tool),
           position = "dodge")

# number of pathogenic and benign variants for each cancer type
num_var <- result_data %>%
  group_by(cancer) %>%
  summarize(count_total = n(),
            count_path = sum(label),
            count_ben = count_total - count_path) %>%
  mutate(cancer = str_to_title(cancer)) %>%
  pivot_longer(cols = c(count_path, count_ben))

# barplot of variant count per cancer type
ggplot(num_var, aes(x = reorder(cancer, -count_total), y = value)) +
  geom_col(aes(fill = name)) +
  labs(title = "Number of Variants in Each Cancer Type",
       x = "Cancer Type",
       y = "Number of Variants") +
  scale_fill_manual(values = c("#6bb3f5", "#fe7575"),
                    name = "Pathogenic/Benign",
                    labels = c("Benign", "Pathogenic")) +
  theme_minimal() +
  theme(axis.title = element_text(size = 15),
        axis.text = element_text(size = 11))

# radar chart
max_min <- table_metrics %>%
  add_row(Tool = "max", Accuracy = 1, Sensitivity = 1, Specificity = 1, 
          PPV = 1.0, NPV = 1, MCC = 1, .before = 1) %>%
  add_row(Tool = "min", Accuracy = 0, Sensitivity = 0, Specificity = 0, 
          PPV = 0, NPV = 0, MCC = 0, .before = 2) %>%
  column_to_rownames('Tool')

create_beautiful_radarchart <- function(data, color = "#00AFBB", axistype = 0,
                                        vlabels = colnames(data), vlcex = 0.7,
                                        caxislabels = NULL, title = NULL, ...){
  radarchart(
    data, axistype = axistype,
    # Customize the polygon
    pcol = color, pfcol = scales::alpha(color, 0.5), plwd = 2, plty = 1,
    # Customize the grid
    cglcol = "grey", cglty = 1, cglwd = 0.8,
    # Customize the axis
    axislabcol = "grey",
    # Variable labels
    vlcex = vlcex, vlabels = vlabels,
    caxislabels = caxislabels, title = title, ...
  )
}

opar <- par()
par(mar=rep(1.2,4))
par(mfrow=c(2,4))

for (i in 1:7){
  create_beautiful_radarchart(
    data = max_min[c(1, 2, i+2), ], 
    title = rownames(max_min)[i+2],
    axistype = 1,
    vlabels = c("Accuracy", "Sens", "Spec", "PPV", "NPV", "MCC"),
    caxislabels = c("0.00", "0.25", "0.50", "0.75", "1.00"),
    vlcex = 1.3,
    calcex = 1.1,
    cex.main = 1.7,
    centerzero = TRUE
  )
}

par <- par(opar)
