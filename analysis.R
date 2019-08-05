# This R code reproduces some results presented in Jevsenak J., 2019. 
# Daily climate data reveal stronger climate-growth relationships for 
# an extended European tree-ring network. Quaternary Science Reviews.

# Open Table S1
calculations <- read.csv("Table_S1.csv") 

# In table calculations, there are all calculated correlation coefficients for 
# the three climate variables: average temperature, sum of precipitation and SPEI. 

# Load required R packages
library("lubridate")
library("dplyr")
library("ggplot2")
library("reshape2")

##########################
# Results in Section 3.1 #
##########################

# Range of elevations 
calculations %>% summarise(elevation_min = min(elevation, na.rm = TRUE),
                           elevation_max = max(elevation, na.rm = TRUE))

# Share of Frits Schweingruber's chronologies
calculations %>% distinct(key_code, .keep_all = TRUE) %>%
  mutate(Schweingruber = grepl("Schweingrub", contributor)) %>% summarise(mean(Schweingruber))

# Species representativnes
calculations %>% distinct(key_code, .keep_all = TRUE) %>%
  dplyr::group_by(species_code) %>% summarise(N = n()) %>% arrange(-N)

# Types of measurements with shares
calculations %>% filter(boot == "Without Bootstrap") %>%
  dplyr::group_by(measurement_type) %>% 
        summarise(share = n()/nrow(dplyr::filter(calculations, boot == "Without Bootstrap"))) %>% 
        arrange(-share)

# Conifers/broadleaves
calculations %>% filter(boot == "Without Bootstrap") %>%
  dplyr::group_by(brd_con) %>% 
  summarise(share = n()/nrow(dplyr::filter(calculations, boot == "Without Bootstrap"))) %>%
  arrange(-share)

# The number of years included in analysis
calculations %>% filter(boot == "Without Bootstrap") %>%
  mutate(n_years = end_year - start_year + 1) %>% 
  summarise(min(n_years), 
            max(n_years),
            mean(n_years))

##################################
# Results in Section 3.2 and 3.3 #
##################################

# Are signs egual? New variable is added: equal 
calculations <- calculations %>%
  mutate(equal = if_else((calculation_daily < 0 & calculation_monthly < 0) | (calculation_daily > 0 & calculation_monthly > 0), 1,0))

# I calculate the overlap of daily and monthly time windows. New variable is added: overlap
calculations$overlap <- NA

# Calculate the overlapp of daily and monthly periods
for (i in 1:nrow(calculations)){
  
  daily_d <- seq(from = calculations[i, "onset_daily"], to = calculations[i, "end_daily"], by = 1)
  monthly_d <- seq(from = calculations[i, "onset_monthly"], to = calculations[i, "end_monthly"], by = 1)
  
  overlap <- sum(daily_d %in% monthly_d)
  
  calculations[i, "overlap"] <- overlap 
}

# The number and share of examples, where signs are not equal and/or overlap is less than 6
dplyr::filter(calculations, equal == FALSE & overlap < 7) %>% 
  filter(boot == "Without Bootstrap") %>% 
  group_by(clim_var) %>% 
  summarise(share_daily = n() / 1860)

# Comparisson of daily and monthly correlations 
# (separately for bootstrapped and nonbootstrapped correlations)
dplyr::filter(calculations, equal == TRUE & overlap > 6) %>% 
  group_by(boot, clim_var) %>% 
  summarise(n_daily = n(),
            
            daily_avg = mean(abs(calculation_daily)),
            daily_sd = sd(abs(calculation_daily)),
            
            monthly_avg = mean(abs(calculation_monthly)),
            monthly_sd = sd(abs(calculation_monthly))) %>% 
  
  mutate(diff = daily_avg - monthly_avg)

# Comparisson of daily and monthly correlations - grouped by type of proxies
# (separately for bootstrapped and nonbootstrapped correlations)
dplyr::filter(calculations, equal == TRUE & overlap > 6) %>% 
  group_by(measurement_type, clim_var, boot) %>% 
  summarise(n_daily = n(),
            
            daily_avg = mean(abs(calculation_daily)),
            daily_sd = sd(abs(calculation_daily)),
            
            monthly_avg = mean(abs(calculation_monthly)),
            monthly_sd = sd(abs(calculation_monthly))) %>% 
  
  mutate(diff = daily_avg - monthly_avg)

# What is the share of chronologies used for bootstrapp?
calculations %>% 
  group_by(boot) %>% 
  summarise(n = n(),
            share_daily = n / (1860 * 3))

# What is the share of bootstrapp correlations?
(nrow(dplyr::filter(calculations, boot == "Bootstrap")))/(nrow(dplyr::filter(calculations, boot == "Without Bootstrap")))

# What are differences reported for Rsquared
dplyr::filter(calculations, equal == TRUE & overlap > 6) %>% 
  group_by(boot, clim_var) %>%
  mutate(diff = abs(calculation_daily) - abs(calculation_monthly),
         diff_squared = calculation_daily^2 - calculation_monthly^2
         ) %>%
  summarise(
    mean(diff),
    mean(diff_squared)*100)
  
# Histograms of daily and monthly approach
hist_df <- dplyr::select(calculations, key_code, calculation_monthly, calculation_daily, clim_var, boot)
hist_df <- melt(hist_df)
hist_df <- mutate(hist_df, variable = ifelse(variable == "calculation_monthly", "Monthly Approach", "Daily Approach"),
                clim_var = ifelse(clim_var == "avg_temperatures", "Temperature Data", 
                                  ifelse(clim_var == "SPEI", "SPEI Data", "Precipitation Data")),
                clim_var = factor(clim_var, levels = c("Temperature Data", "Precipitation Data", "SPEI Data"))
                
                )


ggplot(hist_df, aes(x = abs(value))) +
  geom_histogram(bins = 50) + facet_grid(variable + boot ~ clim_var) + theme_bw() +
  xlab("Correlation Coefficient") +
  geom_vline(data=aggregate(abs(hist_df[5]), hist_df[c(2,3,4)], mean), 
             mapping=aes(xintercept=value), color="red")+
  theme(strip.text = element_text(size = 10))

ggsave(filename="Figure_5_histograms.png", width = 6, height = 7)

#########################################################################

# Histograms of differences
calculations <- dplyr::mutate(calculations, diff = abs(calculation_daily) - abs(calculation_monthly),
                              diff_squared = abs(calculation_daily)^2 - abs(calculation_monthly)^2)

ggplot((dplyr::filter(calculations, boot == "Without Bootstrap") %>% filter(equal == TRUE & overlap > 6)), aes(diff)) + geom_histogram() +
  facet_grid(. ~ clim_var) + theme_bw() +
  xlab("Differences between absolute day-wise and month-wise\n aggregated correlation coefficient ")

ggsave(filename="Figure_3_Histograms_differences.png", width = 6, height = 4)

# Summary statistic of differences
dplyr::filter(calculations, boot == "Without Bootstrap") %>% 
  filter(equal == TRUE & overlap > 6) %>%  
  group_by(clim_var) %>%
  summarise(n(), mean(diff), median(diff), round(sd(diff), 6), max(diff), min(diff),
            mean(diff_squared), median(diff_squared), round(sd(diff_squared), 6), 
            max(diff_squared), min(diff_squared))

# Comparison of confidence intervals for daily and monthly approach
calculations <- arrange(calculations, -abs(calculation_daily))
calculations <- mutate(calculations, clim_var = ifelse(clim_var == "avg_temperatures", 
                                                       "Temperature Data", 
                                    ifelse(clim_var == "SPEI", "SPEI Data", 
                                           "Precipitation Data")),
                  clim_var = factor(clim_var, 
                                    levels = c("Temperature Data", "Precipitation Data", 
                                               "SPEI Data"))
                  
)


# This function ensures the same order of factors on plot as in table
spre <- function(x) factor(x, levels = rev(unique(x)))

ggplot(dplyr::filter(calculations, equal == TRUE & overlap > 6), aes(y = spre(key_code), yend = spre(key_code), 
                                    x = abs(upper_CI_daily), xend = abs(lower_CI_daily))) + 
  geom_segment(col = "blue", alpha = 0.5) +
  geom_segment(aes(y = spre(key_code), yend = spre(key_code), 
                   x = abs(upper_CI_monthly), xend = abs(lower_CI_monthly)),
               col = "red", alpha = 0.5) +
  theme_bw() +
  xlab("95 % Confidence Intervals for Absolute Correlation Coefficients") +
  theme(axis.text.y=element_blank(),axis.ticks.y=element_blank(), axis.title.y = element_blank(),
        panel.grid.minor = element_blank(), panel.grid.major = element_blank()) + 
  facet_grid(.~clim_var) +
  scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1))

ggsave(filename="Figure_6_OverlapPlot.png")

# In the folowing step, we will calculate the overlap of CI for monthly and daily approach
calculations$share_overlap_daily <- NA
calculations$share_overlap_monthly <- NA

for (i in 1:nrow(calculations)){

  
  if (is.na(calculations[i, "upper_CI_daily"])) {
    next
  }

  sekvenca <- c(
    
    a <- abs(calculations[i, "upper_CI_daily"]),
    b <- abs(calculations[i, "lower_CI_daily"]),
    c <- abs(calculations[i, "upper_CI_monthly"]),
    d <- abs(calculations[i, "lower_CI_monthly"])
  )
  
  sekvenca <- sort(sekvenca, decreasing = FALSE)
  
  first_dist <- sekvenca[3] - sekvenca[1]
  second_dist <- sekvenca[4] - sekvenca[2] 
  third_dist <- sekvenca[3] - sekvenca[2] 
  
  daily_dist <- abs(abs(calculations[i, "upper_CI_daily"]) - abs(calculations[i, "lower_CI_daily"]))
  monthly_dist <- abs(abs(calculations[i, "upper_CI_monthly"]) - abs(calculations[i, "lower_CI_monthly"]))

  share_overlap_daily <-  third_dist / daily_dist
  share_overlap_monthly <-  third_dist / monthly_dist
  
  calculations[i, "share_overlap_daily"] <- share_overlap_daily
  calculations[i, "share_overlap_monthly"] <- share_overlap_monthly
  
}

# Plot shares of overlap
ggplot(calculations, aes(share_overlap_daily)) + geom_histogram(alpha = 0.5, col = "blue") + theme_bw() +
  geom_histogram(aes(share_overlap_monthly), col = "red", alpha = 0.5) +
  xlab("Share of overlap") 

ggsave(filename="Figure_7_plot_shares.png")

# In how many cases, overlap is less then 0.5
# absolute numbers
sum(calculations$share_overlap_daily < 0.5, na.rm = TRUE)
sum(calculations$share_overlap_monthly < 0.5, na.rm = TRUE)

# shares
sum(calculations$share_overlap_daily < 0.5, na.rm = TRUE) / sum(!is.na(calculations$share_overlap_daily)) * 100
sum(calculations$share_overlap_monthly < 0.5, na.rm = TRUE) / sum(!is.na(calculations$share_overlap_monthly)) *100

# Share of overlap greater than 60 %
sum(calculations$share_overlap_daily >= 0.6, na.rm = TRUE) / sum(!is.na(calculations$share_overlap_daily)) * 100
sum(calculations$share_overlap_monthly >= 0.6, na.rm = TRUE) / sum(!is.na(calculations$share_overlap_monthly)) *100

# Differences between daily and monthly correlations by proxies

diff_together <- dplyr::filter(calculations, equal == TRUE & overlap > 6) %>%
                group_by(boot, clim_var) %>% 
                   summarise(n_daily = n(),
                   daily_avg = mean(abs(calculation_daily)),
                   daily_sd = sd(abs(calculation_daily)),
                   daily_min = min(abs(calculation_daily)),
                   daily_max = max(abs(calculation_daily)),
                   monthly_avg = mean(abs(calculation_monthly)),
                   monthly_sd = sd(abs(calculation_monthly)),
                   monthly_min = min(abs(calculation_monthly)),
                   monthly_max = max(abs(calculation_monthly))) %>% 
              mutate(diff = daily_avg - monthly_avg)

diff_together <- tibble::add_column(diff_together, measurement_type = "All proxies together", .before = 2)

diff_proxies <- dplyr::filter(calculations, equal == TRUE & overlap > 6) %>% 
  group_by(boot, clim_var, measurement_type) %>% 
  summarise(n_daily = n(),
            daily_avg = mean(abs(calculation_daily)),
            daily_sd = sd(abs(calculation_daily)),
            daily_min = min(abs(calculation_daily)),
            daily_max = max(abs(calculation_daily)),
            monthly_avg = mean(abs(calculation_monthly)),
            monthly_sd = sd(abs(calculation_monthly)),
            monthly_min = min(abs(calculation_monthly)),
            monthly_max = max(abs(calculation_monthly))) %>% 
  mutate(diff = daily_avg - monthly_avg)

diffs <- rbind(data.frame(diff_together), data.frame(diff_proxies))

diffs <- dplyr::arrange(diffs, clim_var, measurement_type, boot)

##########################
# Results in Section 3.4 #
##########################

# Data frame for plot
wind_diffs <- filter(calculations, equal == TRUE & overlap > 6) %>%
    filter(boot == "Without Bootstrap") %>% mutate(
  Onset = onset_daily - onset_monthly,
  End = end_daily - end_monthly, 
  Length = length_daily - length_monthly
)

wind_diffs <- dplyr::select(wind_diffs, "key_code","clim_var", "measurement_type", "Onset", "End", "Length")
wind_diffs <- melt(wind_diffs)

ggplot(wind_diffs, aes(value)) + geom_histogram() +
  facet_grid(variable ~ clim_var) + 
  xlab("Differences Between Daily and Monthly Approach") + 
  theme_bw()

ggsave(filename = "Figure_9_win_diffs.png")

# what are sd, mean median?
wind_diffs %>% group_by(clim_var, variable) %>% 
  summarise(mean(value), median(value), sd(value), max(value), min(value))

# Here I create a table for window distributions
wind_dist <- dplyr::filter(calculations, equal == TRUE & overlap > 6) %>% 
  dplyr::filter(boot == "Without Bootstrap") %>%
  dplyr::select("key_code", "length_daily",  "onset_daily", "end_daily", 
                "onset_monthly", "end_monthly", "length_monthly",
                "clim_var", "measurement_type")

wind_dist <- melt(wind_dist)
wind_dist <- mutate(wind_dist, period = ifelse(grepl("daily", variable), "Daily", "Monthly"))

wind_dist %>% group_by(measurement_type, clim_var, variable, period) %>% 
  summarise(average = median(value))

# This function converts the first letter to uppercase
firstup <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}


wind_dist <- dplyr::mutate(wind_dist, var_hold = gsub("_daily","",variable), 
                           var_hold2 = gsub("_monthly","",var_hold), 
                           var_hold2 = firstup(var_hold2),
                           var_hold2 = factor(var_hold2, levels = c("Onset", "End", "Length"))
                           
)

wind_dist <- dplyr::mutate(wind_dist, var_hold = gsub("_daily","",variable), 
                                      var_hold2 = gsub("_monthly","",var_hold), 
                           var_hold2 = firstup(var_hold2),
                           var_hold2 = factor(var_hold2, levels = c("Onset", "End", "Length"))
                           
                           )

ggplot(wind_dist, aes(value,color = period, fill = period)) + 
  geom_density( alpha = 0.5) + 
  theme_bw()  +
  guides(fill=guide_legend(ncol=2)) +
  theme(legend.position = "bottom") +
  facet_grid(clim_var ~ var_hold2) +
  scale_fill_manual(name = "",values = c("red", "blue")) +
  scale_color_manual(values = c("red", "blue"), guide = "none") +
  geom_vline(data=aggregate(abs(wind_dist[5]), wind_dist[c(2,8, 6)], median), 
             mapping=aes(xintercept=value), 
             col = c(rep("red", 9), rep("blue", 9)))

ggsave(filename="Figure_8_Density_plot_win.png")


# What about time windows? Do they differ between bootstrapped and simple approach?
dplyr::filter(calculations, equal == TRUE & overlap > 6) %>% 
  group_by(boot, clim_var) %>% 
  summarise(n_observations = n(),
            onset_daily = median(onset_daily),
            end_daily = median(end_daily),
            length_daily = median(length_daily),
            onset_monthly = median(onset_monthly),
            end_monthly = median(end_monthly),
            length_monthly = median(length_monthly)
            
            ) %>%
  
            arrange(clim_var)
