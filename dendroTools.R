# The following R code was used for the analysis of day-wise and month-wise aggregated correlation 
# coefficients, published by Jevsenak J., 2019. Daily climate data reveal stronger climate-growth 
# relationships for an extended European tree-ring network. Quaternary Science Reviews.

# Part 1/2: Temperatures and precipitation

# Load required R packages
library("dplR")
library("dplyr")
library("lubridate")
library("reshape2")
library("dendroTools")
library("stringi")
library("stringr")

# Transforms E-OBS daily data into matrix suitable for daily_response()
daily_transform <- function(input){
  
  input$date <- ymd(input$date)
  input$year <- year(input$date)
  input$doy <- yday(input$date)
  input$lat <- NULL
  input$date <- NULL
  input$lon <- NULL
  colnames(input)[1] <- "value"
  daily_matrix <- dcast(year ~ doy, data = input, value.var = "value")
  row.names(daily_matrix) <- daily_matrix$year
  daily_matrix$year <- NULL
  daily_matrix
}

# Transforms E-OBS daily data into matrix suitable for monthly_response()
transform_to_monthly <- function(input, fun.aggregate = sum){
  
  input$date <- ymd(input$date)
  input$year <- year(input$date)
  input$month <- month(input$date)
  input$lat <- NULL
  input$date <- NULL
  input$lon <- NULL
  colnames(input)[1] <- "value"
  daily_matrix <- dcast(year ~ month, data = input, value.var = "value", fun.aggregate = fun.aggregate, na.rm = TRUE)
  row.names(daily_matrix) <- daily_matrix$year
  daily_matrix$year <- NULL
  daily_matrix
}

# Extracts information from analysis and converts into a txt file 
step2 <- function(input){
  
  analysed_period_start <- list()
  analysed_period_end <- list()
  calculation <- list()
  
  calculation_lower <- list()
  calculation_upper <- list()
  
  w_length <- list()
  opti_start<- list()
  previous_year <- list()
  
  for (i in 1:length(input)){
    
    splited <- strsplit(input[[i]]$plot_extreme$labels$title, split = "\n")
    analysed_period_start[[i]] <- str_extract_all(splited[[1]][1],"\\(?[0-9,.]+\\)?")[[1]][1]
    analysed_period_end[[i]] <- str_extract_all(splited[[1]][1],"\\(?[0-9,.]+\\)?")[[1]][2]
    
    calculation[[i]] <- as.numeric(as.character(summary(input[[i]])[5, 2]))
    calculation_lower[[i]] <- as.numeric(as.character(summary(input[[i]])[6, 2]))
    calculation_upper[[i]]  <- as.numeric(as.character(summary(input[[i]])[7, 2]))
    
    w_length[[i]] <- as.numeric(as.character(summary(input[[i]])[11, 2]))
    opti_start[[i]] <- str_extract_all(splited[[1]][5],"\\(?[0-9,.]+\\)?")[[1]][1]
    previous_year[[i]] <- grepl("Previous", splited[[1]][5])
    
  }
  
  final_df <- data.frame(matrix(as.numeric(unlist(analysed_period_start))),
                         matrix(as.numeric(unlist(analysed_period_end))),
                         matrix(as.numeric(unlist(calculation))),
                         
                         matrix(as.numeric(unlist(calculation_lower))),
                         matrix(as.numeric(unlist(calculation_upper))),
                         
                         matrix(as.numeric(unlist(w_length))),
                         matrix(as.numeric(unlist(opti_start))),
                         matrix(unlist(previous_year))
  )
  
  colnames(final_df) <- c("start", "end", "calculation", "lower_bound", "upper_bound" ,"w_length", "opti_start", "previous_year")
  
  final_df <- dplyr::mutate(final_df, opti_end = opti_start + w_length - 1)
  
  final_df <- dplyr::select(final_df, start, end, calculation, lower_bound, upper_bound , w_length, opti_start, opti_end, previous_year)
  
  final_df
}

######################################################
# Open chronologies, temperatures and precipitation ##
######################################################

path # path of your meta file
meta_file <- read.table(path, header = TRUE)

###########################################################################
#T daily
results_Tavg_daily <- list()
results_Tavg_monthly <- list()

results_P_daily <- list()
results_P_monthly <- list()

# SPEI is not analysed here, check r file SPEI.R

for (m in 1:nrow(meta_file)){
  
  if (meta_file[m,"source"] == "BACI ISO"){
    
    chronology <- read.table(paste0("C:/Users/JernejJ/Desktop/GDS analysis revision/chronologies/spline32/", meta_file[m,"file_code"], ".txt"))
    chronology[,2] <- NULL
    
  } else {
    
    chronology <- read.crn(paste0("C:/Users/JernejJ/Desktop/GDS analysis revision/chronologies/spline32/", meta_file[m,"file_code"], ".crn"))
    chronology[,2] <- NULL
    
  }
  
  temperatureAVG <- read.table(paste0("C:/Users/JernejJ/Desktop/GDS analysis revision/Climate data/T_", meta_file[m,"key_clim"], ".txt"), header = TRUE)
  temperaturesAVG_daily <- daily_transform(temperatureAVG)
  temperaturesAVG_monthly <- transform_to_monthly(temperatureAVG, fun.aggregate = mean)
  
  precipitation <- read.table(paste0("C:/Users/JernejJ/Desktop/GDS analysis revision/Climate data/P_", meta_file[m,"key_clim"], ".txt"), header = TRUE)
  precipitations_daily <- daily_transform(precipitation)
  precipitations_monthly <- transform_to_monthly(precipitation, fun.aggregate = sum)
  
  results_Tavg_daily[[1]] <- daily_response(response = chronology, env_data = temperaturesAVG_daily, 
                                         lower_limit = 21, upper_limit = 365, 
                                         method = "cor", row_names_subset = TRUE, fixed_width = 0, alpha = 0.5,
                                         remove_insignificant = FALSE, previous_year = FALSE, boot = TRUE, boot_n = 1000)
  
  results_Tavg_monthly[[1]] <- monthly_response(response = chronology, env_data = temperaturesAVG_monthly, 
                                             method = "cor", row_names_subset = TRUE, alpha = 0.5,
                                             remove_insignificant = FALSE, previous_year = FALSE, boot = TRUE, boot_n = 1000)

  
  results_P_daily[[1]] <- daily_response(response = chronology, env_data = precipitations_daily, 
                                         lower_limit = 21, upper_limit = 365, 
                                         method = "cor", row_names_subset = TRUE, alpha = 0.5,
                                         remove_insignificant = FALSE, previous_year = FALSE, boot = TRUE, boot_n = 1000)
  
  results_P_monthly[[1]] <- monthly_response(response = chronology, env_data = precipitations_monthly, 
                                             method = "cor", row_names_subset = TRUE, alpha = 0.5,
                                             remove_insignificant = FALSE, previous_year = FALSE, boot = TRUE, boot_n = 1000)
  


  Tavg_daily <- step2(results_Tavg_daily)
  colnames(Tavg_daily) <- paste0(colnames(Tavg_daily), "_daily")
  
  Tavg_monthly <- step2(results_Tavg_monthly)
  colnames(Tavg_monthly) <- paste0(colnames(Tavg_monthly), "_monthly")
  
  temperatures_Tavg <- cbind(Tavg_daily, Tavg_monthly)
  temperatures_Tavg$variable <- "avg_temperatures"
  
  
  
  P_daily <- step2(results_P_daily)
  colnames(P_daily) <- paste0(colnames(P_daily), "_daily")
  
  P_monthly <- step2(results_P_monthly)
  colnames(P_monthly) <- paste0(colnames(P_monthly), "_monthly")
  
  precipitation <- cbind(P_daily, P_monthly)
  precipitation$variable <- "precipitation"
  
  
  final_results <- rbind(temperatures_Tavg, precipitation)
  
  final_results <- tibble::add_column(final_results, key_code = meta_file[m,"file_code"], .before = 1)
  
  write.table(final_results, gsub(".crn","",paste0("C:/Users/JernejJ/Desktop/GDS analysis revision/calculations_yesBOOT/TP_", meta_file[m,"file_code"], ".txt")))
  
  print(as.character(meta_file[m,"file_code"]))

  }
