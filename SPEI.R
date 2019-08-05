# The following R code was used for the analysis of day-wise and month-wise aggregated correlation 
# coefficients, published by Jevsenak J., 2019. Daily climate data reveal stronger climate-growth 
# relationships for an extended European tree-ring network. Quaternary Science Reviews.

# Part 2/2: SPEI

library("boot")
library("dplR")
library("dplyr")
library("lubridate")
library("reshape2")
library("dendroTools")
library("stringi")
library("stringr")
library("SPEI")
library("lmom")
library("zoo")
library("boot")

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

# Extracts information from analysis and converts into a txt file 
step2_SPEI <- function(input){
  
  analysed_period_start <- list()
  analysed_period_end <- list()
  calculation <- list()
  
  calculation_lower <- list()
  calculation_upper <- list()
  
  w_length <- list()
  opti_start<- list()
  previous_year <- list()
  
  analysed_period_start[[i]] <- substr(as.character(input[4 ,2]), 1, 4)
  analysed_period_end[[i]] <- substr(as.character(input[4 ,2]), 8, 11)
  
  calculation[[i]] <- as.numeric(as.character(input[5, 2]))
  calculation_lower[[i]] <- as.numeric(as.character(input[6, 2]))
  calculation_upper[[i]]  <- as.numeric(as.character(input[7, 2]))
  
  w_length[[i]] <-  as.numeric(as.character(input[11, 2]))
  opti_start[[i]] <- as.numeric(str_extract_all(as.character(input[8, 2]),"\\(?[0-9,.]+\\)?")[[1]][1])
  previous_year[[i]] <- grepl("Previous", as.character(input[8, 2]))
  
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

# This function is copy pasted from GGitHub repository
# https://github.com/sbegueria/SPEI
# It is slightly modified to fit in the monthly_response()
# All credits goes to authors
SPEI_monthly <- function(SPEI_source_file, scale = 1){
  
  SPIEa <- SPEI_f[,-c(1,2) ]
  SPIEa = zooreg(SPIEa, start=c(1950,1))
  SPIEa= as.ts(SPIEa) 
  data = SPIEa
  
  kernel=list(type='rectangular',shift=0)
  distribution='log-Logistic' 
  fit='ub-pwm'
  na.rm=FALSE
  ref.start=NULL
  ref.end=NULL
  x=FALSE
  params=NULL
  
  scale <- as.numeric(scale)
  na.rm <- as.logical(na.rm)
  x <- as.logical(x)
  #if (!exists("data",inherits=F) | !exists("scale",inherits=F)) {
  #	stop('Both data and scale must be provided')
  #}
  if (!(distribution %in% c('log-Logistic', 'Gamma', 'PearsonIII'))) {
    stop('Distrib must be one of "log-Logistic", "Gamma" or "PearsonIII"')
  }
  if (!(fit %in% c('max-lik', 'ub-pwm', 'pp-pwm'))) {
    stop('Method must be one of "ub-pwm" (default), "pp-pwm" or "max-lik"')
  }
  if ( (!is.null(ref.start) && length(ref.start)!=2) | (!is.null(ref.end) && length(ref.end)!=2) ) {
    stop('Start and end of the reference period must be a numeric vector of length two.')
  }
  
  if (!is.ts(data)) {
    data <- ts(as.matrix(data), frequency = 12, start = c(1950, 1))
  } else {
    data <- ts(as.matrix(data), frequency=frequency(data), start=start(data))
  }
  m <- ncol(data)
  fr <- frequency(data)
  
  
  coef = switch(distribution,
                "Gamma" = array(NA,c(2,m,fr),list(par=c('alpha','beta'),colnames(data),NULL)),
                "log-Logistic" = array(NA,c(3,m,fr),list(par=c('xi','alpha','kappa'),colnames(data),NULL)),
                "PearsonIII" = coef <- array(NA,c(3,m,fr),list(par=c('mu','sigma','gamma'),colnames(data),NULL))
  )
  
  dim_one = ifelse(distribution == "Gamma", 2, 3)
  
  if (!is.null(params)) {
    if (dim(params)[1]!=dim_one | dim(params)[2]!=m | dim(params)[3]!=12) {
      stop(paste('parameters array should have dimensions (', dim_one, ', ', m, ', 12)',sep=' '))
    }
  }
  
  # Loop through series (columns in data)
  if (!is.null(ref.start) && !is.null(ref.end)) {
    data.fit <- window(data,ref.start,ref.end)	
  } else {
    data.fit <- data
  }
  std <- data*NA
  for (s in 1:m) {
    # Cumulative series (acu)
    acu <- data.fit[,s]
    acu.pred <- data[,s]
    if (scale>1) {
      wgt <- kern(scale,kernel$type,kernel$shift)
      acu[scale:length(acu)] <- rowSums(embed(acu,scale)*wgt,na.rm=na.rm)
      acu[1:(scale-1)] <- NA
      acu.pred[scale:length(acu.pred)] <- rowSums(embed(acu.pred,scale)*wgt,na.rm=na.rm)
      acu.pred[1:(scale-1)] <- NA
    }
    
    # Loop through the months
    for (c in (1:fr)) {
      # Filter month m, excluding NAs
      f <- which(cycle(acu)==c)
      f <- f[!is.na(acu[f])]
      ff <- which(cycle(acu.pred)==c)
      ff <- ff[!is.na(acu.pred[ff])]
      
      # Monthly series, sorted
      month <- sort.default(acu[f], method="quick")
      
      if (length(month)==0) {
        std[f] <- NA
        next()
      }
      
      if (is.null(params)) {
        month_sd = sd(month,na.rm=TRUE)
        if (is.na(month_sd) || (month_sd == 0)) {
          std[f] <- NA
          next
        }
        
        if(distribution != "log-Logistic"){
          pze <- sum(month==0)/length(month)
          month = month[month > 0]
        }
        
        # Stop early and assign NAs if month's data is length < 4
        if(length(month) < 4){
          std[ff,s] = NA
          coef[,s,c] <- NA
          next
        }
        
        # Calculate probability weighted moments based on fit with lmomco or TLMoments
        pwm = switch(fit,
                     "pp-pwm" = pwm.pp(month,-0.35,0, nmom=3),
                     #pwm.ub(month, nmom=3)
                     TLMoments::PWM(month, order=0:2)
        )
        
        # Check L-moments validity
        lmom <- pwm2lmom(pwm)
        if ( !are.lmom.valid(lmom) || anyNA(lmom[[1]]) || any(is.nan(lmom[[1]])) ){
          next
        }
        
        # lmom fortran functions need specific inputs L1, L2, T3
        # this is handled by lmomco internally with lmorph
        fortran_vec = c(lmom$lambdas[1:2], lmom$ratios[3])
        
        # Calculate parameters based on distribution with lmom then lmomco
        f_params = switch(distribution,
                          "log-Logistic" = tryCatch(lmom::pelglo(fortran_vec), error = function(e){ parglo(lmom)$para }),
                          "Gamma" = tryCatch(lmom::pelgam(fortran_vec), error = function(e){ pargam(lmom)$para }),
                          "PearsonIII" = tryCatch(lmom::pelpe3(fortran_vec), error = function(e){ parpe3(lmom)$para })
        )
        
        # Adjust if user chose log-Logistic and max-lik
        if(distribution == 'log-Logistic' && fit=='max-lik'){
          f_params = parglo.maxlik(month, f_params)$para
        }
      } else {
        
        f_params = as.vector(params[,s,c])
        
      }
      
      # Calculate cdf based on distribution with lmom
      cdf_res = switch(distribution,
                       "log-Logistic" = lmom::cdfglo(acu.pred[ff], f_params),
                       "Gamma" = lmom::cdfgam(acu.pred[ff], f_params),
                       "PearsonIII" = lmom::cdfpe3(acu.pred[ff], f_params)				  				
      )
      
      std[ff,s] = qnorm(cdf_res)
      coef[,s,c] <- f_params
      
      # Adjust if user chose Gamma or PearsonIII
      if(distribution != 'log-Logistic'){ 
        std[ff,s] = qnorm(pze + (1-pze)*pnorm(std[ff,s]))
      }
      
    } # next c (month)
  } # next s (series)
  colnames(std) <- colnames(data)
  
  z <- list(call=match.call(expand.dots=FALSE),
            fitted=std,coefficients=coef,scale=scale,kernel=list(type=kernel$type,
                                                                 shift=kernel$shift,values=kern(scale,kernel$type,kernel$shift)),
            distribution=distribution,fit=fit,na.action=na.rm)
  if (x) z$data <- data
  if (!is.null(ref.start)) z$ref.period <- rbind(ref.start,ref.end)
  
  SPEI_source_file[,3] <- z$fitted
  
  SPEI_source_file
}

# This function is copy pasted from GitHub repository
# https://github.com/sbegueria/SPEI
# It is slightly modified to fit in the daily_response()
# All credits goes to authors
SPEI_daily <- function(SPEI_source_file, scale = 90){
  
  SPIEa <- SPEI_source_file[,-c(1,2, 3) ]
  SPIEa = zooreg(SPIEa, start=as.Date("1950-01-01"))
  SPIEa= as.ts(SPIEa) 
  data = SPIEa
  
  kernel=list(type='rectangular',shift=0)
  distribution='log-Logistic' 
  fit='ub-pwm'
  na.rm=FALSE
  ref.start=NULL
  ref.end=NULL
  x=FALSE
  params=NULL
  
  scale <- as.numeric(scale)
  na.rm <- as.logical(na.rm)
  x <- as.logical(x)
  #if (!exists("data",inherits=F) | !exists("scale",inherits=F)) {
  #	stop('Both data and scale must be provided')
  #}
  if (!(distribution %in% c('log-Logistic', 'Gamma', 'PearsonIII'))) {
    stop('Distrib must be one of "log-Logistic", "Gamma" or "PearsonIII"')
  }
  if (!(fit %in% c('max-lik', 'ub-pwm', 'pp-pwm'))) {
    stop('Method must be one of "ub-pwm" (default), "pp-pwm" or "max-lik"')
  }
  if ( (!is.null(ref.start) && length(ref.start)!=2) | (!is.null(ref.end) && length(ref.end)!=2) ) {
    stop('Start and end of the reference period must be a numeric vector of length two.')
  }
  
  if (!is.ts(data)) {
    data <- ts(as.matrix(data), frequency = 12, start = c(1950, 1))
  } else {
    data <- ts(as.matrix(data), frequency=frequency(data), start=start(data))
  }
  m <- ncol(data)
  fr <- frequency(data)
  
  
  coef = switch(distribution,
                "Gamma" = array(NA,c(2,m,fr),list(par=c('alpha','beta'),colnames(data),NULL)),
                "log-Logistic" = array(NA,c(3,m,fr),list(par=c('xi','alpha','kappa'),colnames(data),NULL)),
                "PearsonIII" = coef <- array(NA,c(3,m,fr),list(par=c('mu','sigma','gamma'),colnames(data),NULL))
  )
  
  dim_one = ifelse(distribution == "Gamma", 2, 3)
  
  if (!is.null(params)) {
    if (dim(params)[1]!=dim_one | dim(params)[2]!=m | dim(params)[3]!=12) {
      stop(paste('parameters array should have dimensions (', dim_one, ', ', m, ', 12)',sep=' '))
    }
  }
  
  # Loop through series (columns in data)
  if (!is.null(ref.start) && !is.null(ref.end)) {
    data.fit <- window(data,ref.start,ref.end)	
  } else {
    data.fit <- data
  }
  std <- data*NA
  for (s in 1:m) {
    # Cumulative series (acu)
    acu <- data.fit[,s]
    acu.pred <- data[,s]
    if (scale>1) {
      wgt <- kern(scale,kernel$type,kernel$shift)
      acu[scale:length(acu)] <- rowSums(embed(acu,scale)*wgt,na.rm=na.rm)
      acu[1:(scale-1)] <- NA
      acu.pred[scale:length(acu.pred)] <- rowSums(embed(acu.pred,scale)*wgt,na.rm=na.rm)
      acu.pred[1:(scale-1)] <- NA
    }
    
    # Loop through the months
    for (c in (1:fr)) {
      # Filter month m, excluding NAs
      f <- which(cycle(acu)==c)
      f <- f[!is.na(acu[f])]
      ff <- which(cycle(acu.pred)==c)
      ff <- ff[!is.na(acu.pred[ff])]
      
      # Monthly series, sorted
      month <- sort.default(acu[f], method="quick")
      
      if (length(month)==0) {
        std[f] <- NA
        next()
      }
      
      if (is.null(params)) {
        month_sd = sd(month,na.rm=TRUE)
        if (is.na(month_sd) || (month_sd == 0)) {
          std[f] <- NA
          next
        }
        
        if(distribution != "log-Logistic"){
          pze <- sum(month==0)/length(month)
          month = month[month > 0]
        }
        
        # Stop early and assign NAs if month's data is length < 4
        if(length(month) < 4){
          std[ff,s] = NA
          coef[,s,c] <- NA
          next
        }
        
        # Calculate probability weighted moments based on fit with lmomco or TLMoments
        pwm = switch(fit,
                     "pp-pwm" = pwm.pp(month,-0.35,0, nmom=3),
                     #pwm.ub(month, nmom=3)
                     TLMoments::PWM(month, order=0:2)
        )
        
        # Check L-moments validity
        lmom <- pwm2lmom(pwm)
        if ( !are.lmom.valid(lmom) || anyNA(lmom[[1]]) || any(is.nan(lmom[[1]])) ){
          next
        }
        
        # lmom fortran functions need specific inputs L1, L2, T3
        # this is handled by lmomco internally with lmorph
        fortran_vec = c(lmom$lambdas[1:2], lmom$ratios[3])
        
        # Calculate parameters based on distribution with lmom then lmomco
        f_params = switch(distribution,
                          "log-Logistic" = tryCatch(lmom::pelglo(fortran_vec), error = function(e){ parglo(lmom)$para }),
                          "Gamma" = tryCatch(lmom::pelgam(fortran_vec), error = function(e){ pargam(lmom)$para }),
                          "PearsonIII" = tryCatch(lmom::pelpe3(fortran_vec), error = function(e){ parpe3(lmom)$para })
        )
        
        # Adjust if user chose log-Logistic and max-lik
        if(distribution == 'log-Logistic' && fit=='max-lik'){
          f_params = parglo.maxlik(month, f_params)$para
        }
      } else {
        
        f_params = as.vector(params[,s,c])
        
      }
      
      # Calculate cdf based on distribution with lmom
      cdf_res = switch(distribution,
                       "log-Logistic" = lmom::cdfglo(acu.pred[ff], f_params),
                       "Gamma" = lmom::cdfgam(acu.pred[ff], f_params),
                       "PearsonIII" = lmom::cdfpe3(acu.pred[ff], f_params)				  				
      )
      
      std[ff,s] = qnorm(cdf_res)
      coef[,s,c] <- f_params
      
      # Adjust if user chose Gamma or PearsonIII
      if(distribution != 'log-Logistic'){ 
        std[ff,s] = qnorm(pze + (1-pze)*pnorm(std[ff,s]))
      }
      
    } # next c (month)
  } # next s (series)
  colnames(std) <- colnames(data)
  
  z <- list(call=match.call(expand.dots=FALSE),
            fitted=std,coefficients=coef,scale=scale,kernel=list(type=kernel$type,
                                                                 shift=kernel$shift,values=kern(scale,kernel$type,kernel$shift)),
            distribution=distribution,fit=fit,na.action=na.rm)
  if (x) z$data <- data
  if (!is.null(ref.start)) z$ref.period <- rbind(ref.start,ref.end)
  
  SPEI_source_file[,4] <- z$fitted
  
  SPEI_source_file
}

#########################################################################
#########################################################################

meta_file <- read.table("C:/Users/Dendro/Desktop/GDS/TableS1_faza3.txt", header = TRUE)

for (m in 1:nrow(meta_file)){
  
  boot_ci_type = "norm"
  boot_conf_int = 0.95
  boot_n = 1000
  boot = TRUE
  cor_method = "pearson"
  upper_limit = 1
  lower_limit = 1
  cor_method = "pearson"
  reference_window = "start"
  method = "cor"
  K = 1
  
  if (meta_file[m,"source"] == "BACI ISO"){
    
    chronology <- read.table(paste0("C:/Users/Dendro/Desktop/GDS/chronologies/", meta_file[m,"file_code"], ".txt"))
    chronology[,2] <- NULL
    
  } else {
    
    chronology <- read.crn(paste0("C:/Users/Dendro/Desktop/GDS/chronologies/", meta_file[m,"file_code"], ".crn"))
    chronology[,2] <- NULL
    
  }
  
  SPEI_f <- read.table(paste0("C:/Users/Dendro/Desktop/GDS/Climate data/SPEI_", meta_file[m,"key_clim"], ".txt"), header = TRUE)
  
  if (sum(is.na(SPEI_f[, 4]))>0){
    
    SPEI_f <- tidyr::fill(SPEI_f, SPEI)
    
  }
  
  temporal_matrix_list <- list()
  temporal_matrix_upper_list <- list()
  temporal_matrix_lower_list <- list()
  
  response = chronology
  
  b = 1
  
  for (i in (21:365)){
    
    SPEI_temp <- SPEI_daily(SPEI_source_file = SPEI_f, scale = i)
    
    vector <- SPEI_temp[, 4]
    nNA <- sum(is.na(vector))
    
    if (nNA > i - 1){
      
      stop("NA problem!!!!")
      
    }
    
    SPEI_temp <- SPEI_temp[c(1:(nrow(SPEI_temp)-nNA)),c(1,2,3)]
    
    SPEI_temp$values <- vector[(nNA+1):length(vector)]
    
    SPEIX_daily <- daily_transform(SPEI_temp)

    env_data = SPEIX_daily
    
    ncol_response <- ncol(response)
    
    colnames_response <- colnames(response)
    
    env_data$yearABC <- row.names(env_data)
    response$yearABC <- row.names(response)
    
    temporal_data <- merge(response, env_data, by = "yearABC")
    
    response <- data.frame(temporal_data[, c(2:(1 + ncol_response))],
                           row.names = temporal_data$yearABC)
    colnames(response) <- colnames_response
    
    env_data <- data.frame(temporal_data[, c((1 + ncol_response + 1):
                                               ncol(temporal_data))],
                           row.names = temporal_data$yearABC)
    
    
    temporal_matrix <- matrix(NA, nrow = 1,
                              ncol = (ncol(env_data) - 21) + 1)
    # Here I create two additional temporal matrices, which will be used to store
    # lower and upper limits of bootstrap estimates
    temporal_matrix_lower <- temporal_matrix
    temporal_matrix_upper <- temporal_matrix
    
    for (j in 0: (ncol(env_data) - i)) {
      
      x <- env_data[, (j + 1)]
      
      x <- matrix(x, nrow = nrow(env_data), ncol = 1)
      
      if (boot == FALSE){
        temporal_correlation <- cor(response[, 1], x[, 1], method = cor_method)
        temporal_lower <- NA
        temporal_upper <- NA
      } else if (boot == TRUE){
        temp_df_boot <- cbind(response[, 1], x[, 1])
        
        temp_df_boot <- temp_df_boot[complete.cases(temp_df_boot), ]
        
        calc <- boot(temp_df_boot, boot_f, fun = "cor", cor.type = cor_method, R = boot_n)
        
        temporal_correlation <- colMeans(calc$t)[1]
        
        ci_int <- boot.ci(calc, conf = boot_conf_int, type = boot_ci_type)
        temporal_lower <- ci_int$norm[2]
        temporal_upper <- ci_int$norm[3]
      } else {
        print(paste0("boot should be TRUE or FALSE, instead it is ", boot))
      }
      
      temporal_matrix[1, j + 1] <- temporal_correlation
      temporal_matrix_lower[1, j + 1] <- temporal_lower
      temporal_matrix_upper[1, j + 1] <- temporal_upper
      
    }
    
    temporal_matrix_list[[b]] <- temporal_matrix
    temporal_matrix_lower_list[[b]] <- temporal_matrix_lower
    temporal_matrix_upper_list[[b]] <- temporal_matrix_upper
    b = b +1
  
  }
  
  temporal_matrix <- do.call(rbind, temporal_matrix_list)
  temporal_matrix_lower <- do.call(rbind, temporal_matrix_lower_list)
  temporal_matrix_upper <- do.call(rbind, temporal_matrix_upper_list)
  
  temporal_rownames <- as.vector(seq(from = 21, to = 365,
                                     by = 1))
  row.names(temporal_matrix) <- temporal_rownames
  row.names(temporal_matrix_lower) <- temporal_rownames
  row.names(temporal_matrix_upper) <- temporal_rownames
  
  temporal_colnames <- as.vector(seq(from = 1,
                                     to = ncol(temporal_matrix), by = 1))
  colnames(temporal_matrix) <- temporal_colnames
  colnames(temporal_matrix_lower) <- temporal_colnames
  colnames(temporal_matrix_upper) <- temporal_colnames
  
  #########################################################################################
  #########################################################################################
  
  if(is.finite(mean(temporal_matrix, na.rm = TRUE)) == FALSE){
    stop("All calculations are insignificant! Change the alpha argument!")
  }
  
  overall_max <- max(temporal_matrix, na.rm = TRUE)
  overall_min <- min(temporal_matrix, na.rm = TRUE)
  
  # absolute vales of overall_maximum and overall_minimum are compared and
  # one of the following two if functions is used
  # There are unimportant warnings produced:
  # no non-missing arguments to max; returning -Inf
  
  if ((abs(overall_max) > abs(overall_min)) == TRUE) {
    
    # maximum value is located. Row indeces are needed to query information
    # about the window width used to calculate the maximum. Column name is
    # needed to query the starting day.
    max_result <- suppressWarnings(which.max(apply(temporal_matrix,
                                                   MARGIN = 2, max,
                                                   na.rm = TRUE)))
    plot_column <- max_result
    max_index <- which.max(temporal_matrix[, names(max_result)])
    row_index <- row.names(temporal_matrix)[max_index]
  }
  
  if ((abs(overall_max) < abs(overall_min)) == TRUE) {
    
    min_result <- suppressWarnings(which.min(apply(temporal_matrix,
                                                   MARGIN = 2, min,
                                                   na.rm = TRUE)))
    plot_column <- min_result
    min_index <- which.min(temporal_matrix[, names(min_result)])
    row_index <- row.names(temporal_matrix)[min_index]
  }
  
  # Element 5
  # Here we create the fifth element of the final list: Analysed period in the
  # form of min(year) - max(year), e.g. 1950 - 2015
  min_env_data <- min(as.numeric(row.names(env_data)))
  min_response <- min(as.numeric(row.names(response)))
  
  max_env_data <- max(as.numeric(row.names(env_data)))
  max_response <- max(as.numeric(row.names(response)))
  
  min_together <- min(min_env_data, min_response)
  max_together <- min(max_env_data, max_response)
  
  
  analysed_period <- paste(as.character(min_together),
                           as.character(max_together),
                           sep = " - ")
  if (nchar(analysed_period) < 9) {
    analysed_period <- NA
  }
  
  
  # String for titles
  if (method == "cor"){
    title_string <- "Correlation Coefficients"
  } else if (method == "lm"){
    title_string <- "Linear Regression"
  } else if (method == "brnn"){
    title_string <- "ANN With Bayesian Regularization"
  } else (print("The selection of method is not correct"))
  
  
  
  
  final_list <- list(calculations = temporal_matrix, method = method,
                     metric = cor_method, analysed_period = analysed_period,
                     optimized_return = NA,
                     optimized_return_all = NA,
                     transfer_function = NA, temporal_stability = NA,
                     cross_validation = NA,
                     plot_heatmap = NA,
                     plot_extreme = NA,
                     plot_specific = NA,
                     PCA_output = NA,
                     type = "daily",
                     reference_window = "start",
                     boot_lower = temporal_matrix_lower,
                     boot_upper = temporal_matrix_upper)
  
  #####################################################################
  #####################################################################
  #####################################################################
  
  # A) Extracting a matrix from a list and converting it into a data frame
  result_daily_response <- final_list
  object <- final_list
  
  type <- data.frame(object[[14]])
  
  result_daily_element1 <- data.frame(object[[1]])
  
  reference_window <- object[[15]]
  
  # To keep RCMD check happy:
  
  # With the following chunk, overall_maximum and overall_minimum values of
  # result_daily_element1 matrix are calculated.
  overall_max <- max(result_daily_element1, na.rm = TRUE)
  overall_min <- min(result_daily_element1, na.rm = TRUE)
  
  # absolute vales of overall_maximum and overall_minimum are compared and
  # one of the following two if functions is used
  # There are unimportant warnings produced:
  # no non-missing arguments to max; returning -Inf
  # Based on the answer on the StackOverlow site:
  # https://stackoverflow.com/questions/24282550/no-non-missing-arguments-warning-when-using-min-or-max-in-reshape2
  # Those Warnings could be easily ignored
  if ((abs(overall_max) > abs(overall_min)) == TRUE) {
    
    # maximum value is located. Row indeces are needed to query information
    # about the window width used to calculate the maximum. Column name is
    # needed to query the starting day.
    max_result <- suppressWarnings(which.max(apply(result_daily_element1,
                                                   MARGIN = 2, max, na.rm = TRUE)))
    plot_column <- max_result
    plot_column_source <- plot_column
    max_index <- which.max(result_daily_element1[, names(max_result)])
    row_index <- row.names(result_daily_element1)[max_index]
    temporal_vector <- unlist(result_daily_element1[max_index, ])
    temporal_vector <- data.frame(temporal_vector)
    calculated_metric <- round(max(temporal_vector, na.rm = TRUE), 3)
    
    lower_bound <- final_list$boot_lower[max_index, as.numeric(max_result)]
    upper_bound <- final_list$boot_upper[max_index, as.numeric(max_result)]
    
    # Here we remove missing values at the end of the temporal_vector.
    # It is important to remove missing values only at the end of the
    # temporal_vector!
    row_count <- nrow(temporal_vector)
    delete_rows <- 0
    while (is.na(temporal_vector[row_count, ] == TRUE)){
      delete_rows <- delete_rows + 1
      row_count <-  row_count - 1
    }
    # To check if the last row is a missing value
    if (is.na(temporal_vector[nrow(temporal_vector), ] == TRUE)) {
      temporal_vector <-  temporal_vector[-c((row_count + 1):(row_count +
                                                                delete_rows)), ]
    }
    temporal_vector <- data.frame(temporal_vector)
  }
  
  if ((abs(overall_max) < abs(overall_min)) == TRUE) {
    
    # minimum value is located. Row indeces are needed to query information
    # about the window width used to calculate the minimum. Column name is
    # needed to query the starting day.
    min_result <- suppressWarnings(which.min(apply(result_daily_element1,
                                                   MARGIN = 2, min, na.rm = TRUE)))
    plot_column <- min_result
    plot_column_source <- plot_column
    min_index <- which.min(result_daily_element1[, names(min_result)])
    row_index <- row.names(result_daily_element1)[min_index]
    temporal_vector <- unlist(result_daily_element1[min_index, ])
    temporal_vector <- data.frame(temporal_vector)
    calculated_metric <- round(min(temporal_vector, na.rm = TRUE), 3)
    
    lower_bound <- final_list$boot_lower[min_index, as.numeric(min_result)]
    upper_bound <- final_list$boot_upper[min_index, as.numeric(min_result)]
    # Here we remove missing values
    # We remove missing values at the end of the temporal_vector.
    # It is important to remove missing values only at the end of the
    # temporal_vector!
    
    row_count <- nrow(temporal_vector)
    delete_rows <- 0
    while (is.na(temporal_vector[row_count, ] == TRUE)){
      delete_rows <- delete_rows + 1
      row_count <-  row_count - 1
    }
    # To check if the last row is a missing value
    if (is.na(temporal_vector[nrow(temporal_vector), ] == TRUE)) {
      temporal_vector <-  temporal_vector[-c((row_count + 1):(row_count +
                                                                delete_rows)), ]
    }
    temporal_vector <- data.frame(temporal_vector)
  }
  
  # In case of previous_year == TRUE, we calculate the day of a year
  # (plot_column), considering 366 days of previous year.
  if (nrow(temporal_vector) > 366 & plot_column > 366) {
    previous_year = TRUE
    plot_column_extra <- plot_column %% 366
  } else {
    previous_year = FALSE
    plot_column_extra <- plot_column
  }
  
  
  if (nrow(temporal_vector) > 366) {
    previous_year <- TRUE
  } else {
    previous_year <- FALSE
  }
  
  # Here we define a data frame of dates and corresponing day of year (doi). Later
  # this dataframe will be used to describe tht optimal sequence od days
  doy <- seq(1:730)
  date <- seq(as.Date('2013-01-01'),as.Date('2014-12-31'), by = "+1 day")
  # date[366] <- as.Date('2015-12-31')
  date <- format(date, "%b %d")
  date_codes <- data.frame(doy = doy, date = date)
  
  # Here, there is a special check if optimal window width is divisible by 2 or not.
  if (as.numeric(row_index)%%2 == 0){
    adjustment_1 = 0
    adjustment_2 = 1
  } else {
    adjustment_1 = 1
    adjustment_2 = 2
  }
  
  if (reference_window == "start"){
    Optimal_string <- paste(as.character(date_codes[plot_column_extra, 2]),"-",
                            as.character(date_codes[plot_column_extra + as.numeric(row_index) - 1, 2]))
  } else if (reference_window == "end") {
    Optimal_string <- paste(as.character(date_codes[plot_column_extra - as.numeric(row_index) + 1, 2]),"-",
                            as.character(date_codes[plot_column_extra, 2]))
  } else if (reference_window == "middle") {
    Optimal_string <- paste(as.character(date_codes[(round2((plot_column_extra - as.numeric(row_index)/2)) - adjustment_1), 2]),"-",
                            as.character(date_codes[(round2((plot_column_extra + as.numeric(row_index)/2)) - adjustment_2), 2]))
  }
  
  # Here we define titles. They differ importantly among methods and arguments
  # in the final output list from daily_response() function
  if (result_daily_response[[2]] == "cor"){
    y_lab <- NA
  } else if (result_daily_response[[2]] == "pcor"){
    y_lab <- NA
  } else if (result_daily_response[[3]] == "r.squared"){
    y_lab <- "Explained Variance"
  } else if (result_daily_response[[3]] == "adj.r.squared"){
    y_lab <- "Adjusted Explained Variance"
  }
  
  if (reference_window == 'start' &&  plot_column > 366 && nrow(temporal_vector) > 366){
    reference_string <- paste0("Starting Day of Optimal Window Width: Day ",
                               plot_column_extra, " of Current Year")}
  
  if (reference_window == 'start' &&  plot_column <= 366 && nrow(temporal_vector) > 366){
    reference_string <- paste0("Starting Day of Optimal Window Width: Day ",
                               plot_column_extra, " of Previous Year")}
  
  if (reference_window == 'start' &&  plot_column <=  366 && nrow(temporal_vector) <=  366){
    reference_string <- paste0("Starting Day of Optimal Window Width: Day ",
                               plot_column_extra)}
  
  if (reference_window == 'end' &&  plot_column > 366 && nrow(temporal_vector) > 366){
    reference_string <- paste0("Ending Day of Optimal Window Width: Day ",
                               plot_column_extra, " of Current Year")}
  
  if (reference_window == 'end' &&  plot_column  <= 366 && nrow(temporal_vector) > 366){
    reference_string <- paste0("Ending Day of Optimal Window Width: Day ",
                               plot_column_extra, " of Previous Year")}
  
  if (reference_window == 'end' &&  plot_column  <=  366 && nrow(temporal_vector) <=  366){
    reference_string <- paste0("Ending Day of Optimal Window Width: Day ",
                               plot_column_extra)}
  
  
  if (reference_window == 'middle' &&  plot_column > 366 && nrow(temporal_vector) > 366){
    reference_string <- paste0("Middle Day of Optimal Window Width: Day ",
                               plot_column_extra, " of Current Year")}
  
  if (reference_window == 'middle' &&  plot_column  <= 366 && nrow(temporal_vector) > 366){
    reference_string <- paste0("Middle Day of Optimal Window Width: Day ",
                               plot_column_extra, " of Previous Year")}
  
  if (reference_window == 'middle' &&  plot_column  <=  366 && nrow(temporal_vector) <=  366){
    reference_string <- paste0("Middle Day of Optimal Window Width: Day ",
                               plot_column_extra)}
  
  optimal_window_string <- paste0("Optimal Window Width: ", as.numeric(row_index),
                                  " Days")
  
  optimal_calculation <- paste0("The Highest ", y_lab,": " , calculated_metric)
  
  period_string <- paste0("Analysed Period: ", result_daily_response[[4]])
  
  if (result_daily_response[[2]] == 'cor'){
    method_string <- paste0("Correlation Coefficient (", result_daily_response[[3]], ")")
    
  } else if (result_daily_response[[2]] == 'pcor'){
    method_string <- paste0("Partial Correlation Coefficient (", result_daily_response[[3]], ")")
    
  } else if (result_daily_response[[2]] == 'lm'){
    method_string <- paste0("Linear Regression")
  } else if (result_daily_response[[2]] == 'brnn'){
    method_string <- paste0("ANN with Bayesian Regularization")
  }
  
  if (type == "monthly"){
    
    # Plural or singular?
    if (as.numeric(row_index) == 1){
      month_string <- " Month"
      
    } else {
      month_string <- " Months"
    }
    
    # In case of previous_year == TRUE, we calculate the day of a year
    # (plot_column), considering 366 days of previous year.
    
    if (ncol(result_daily_response[[1]]) == 24 & plot_column > 12) {
      plot_column_extra <- plot_column %% 12
    } else {
      plot_column_extra <- plot_column
    }
    
    
    if (ncol(result_daily_response[[1]]) == 24) {
      previous_year <- TRUE
    } else {
      previous_year <- FALSE
    }
    
    
    if (reference_window == 'start' &&  plot_column > 12 && ncol(result_daily_response[[1]]) == 24){
      reference_string <- paste0("Starting Month of Optimal Window Width: Month ",
                                 plot_column_extra, " of Current Year")}
    
    if (reference_window == 'start' &&  plot_column <= 12 && ncol(result_daily_response[[1]]) == 24){
      reference_string <- paste0("Starting Month of Optimal Window Width: Month ",
                                 plot_column_extra, " of Previous Year")}
    
    if (reference_window == 'start' &&  plot_column <=  12 && ncol(result_daily_response[[1]]) <=  12){
      reference_string <- paste0("Starting Month of Optimal Window Width: Month ",
                                 plot_column_extra)}
    
    
    
    optimal_window_string <- paste0("Optimal Window Width: ", as.numeric(row_index),
                                    month_string)
    
    # Here we define a data frame of months. Later
    # this dataframe will be used to describe tht optimal sequence od days
    
    if (ncol(result_daily_response[[1]]) == 24){
      date_codes <- c("Jan*", "Feb*", "Mar*", "Apr*", "May*", "Jun*", "Jul*", "Aug*", "Sep*", "Oct*", "Nov*", "Dec*",
                      "Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")
      
    } else if (ncol(result_daily_response[[1]]) == 12){
      
      date_codes <- c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")
      
    }
    
    
    if (reference_window == "start"){
      Optimal_string <- paste(as.character(date_codes[plot_column_source]),"-",
                              as.character(date_codes[plot_column_source + as.numeric(row_index) - 1]))
    } else if (reference_window == "end") {
      Optimal_string <- paste(as.character(date_codes[plot_column_source - as.numeric(row_index) + 1]),"-",
                              as.character(date_codes[plot_column_source]))
    } else if (reference_window == "middle") {
      Optimal_string <- paste(as.character(date_codes[(round2((plot_column_source - as.numeric(row_index)/2)) - adjustment_1)]),"-",
                              as.character(date_codes[(round2((plot_column_source + as.numeric(row_index)/2)) - adjustment_2)]))
    }
    
    if (as.numeric(row_index == 1)){
      Optimal_string <- substr(Optimal_string, 1, nchar(Optimal_string)-6)
    }
    
  }
  
  output_df <-  data.frame(Variable = c("approach",
                                        "method",
                                        "metric",
                                        "analysed_years",
                                        "maximal_calculated_metric",
                                        "lower_ci",
                                        "upper_ci",
                                        "reference_window",
                                        "analysed_previous_year",
                                        "optimal_time_window",
                                        "optimal_time_window_length"),
                           
                           Value = c(result_daily_response[[14]],
                                     method_string,
                                     y_lab,
                                     result_daily_response[[4]],
                                     calculated_metric,
                                     round(lower_bound, 3),
                                     round(upper_bound, 3),
                                     reference_string,
                                     previous_year,
                                     Optimal_string,
                                     as.numeric(row_index)))
  
######################################### Extract results 
  
  SPEIa_daily <- step2_SPEI(output_df)
  colnames(SPEIa_daily) <- paste0(colnames(SPEIa_daily), "_daily")
  
  
#########################################################

# Here I start with monthly SPEI calculations

  boot_ci_type = "norm"
  boot_conf_int = 0.95
  boot_n = 1000
  boot = TRUE
  cor_method = "pearson"
  upper_limit = 1
  lower_limit = 1
  cor_method = "pearson"
  reference_window = "start"
  method = "cor"
  K = 1

  if (meta_file[m,"source"] == "BACI ISO"){
    
    chronology <- read.table(paste0("C:/Users/Dendro/Desktop/GDS/chronologies/", meta_file[m,"file_code"], ".txt"))
    chronology[,2] <- NULL
    
  } else {
    
    chronology <- read.crn(paste0("C:/Users/Dendro/Desktop/GDS/chronologies/", meta_file[m,"file_code"], ".crn"))
    chronology[,2] <- NULL
    
  }
  
  SPEI_f <- read.table(paste0("C:/Users/Dendro/Desktop/GDS/Climate data/SPEI_", meta_file[m,"key_clim"], ".txt"), header = TRUE)
  
  if (sum(is.na(SPEI_f[, 4]))>0){
    
    SPEI_f <- tidyr::fill(SPEI_f, SPEI)
    
  }
  
  SPEI_f <- SPEI_f[,-c(1,2) ]
  
  SPEI_f <- dplyr::mutate(SPEI_f, date = ymd(date),
                          year = year(date),
                          month = month(date))
  
  SPEI_f <- dplyr::group_by(SPEI_f, year, month) %>% summarise(value = mean(SPEI, na.rm = TRUE))
  
  SPEI_f <- data.frame(SPEI_f)
  
  temporal_matrix_list <- list()
  temporal_matrix_upper_list <- list()
  temporal_matrix_lower_list <- list()
  
  response = chronology
  
  b = 1
  
  for (i in (1:9)){
    
    SPEI_temp <- SPEI_monthly(SPEI_source_file = SPEI_f, scale = i)
    colnames(SPEI_temp)[3] <- "values"
    
    vector <- SPEI_temp[, 3]
    nNA <- sum(is.na(vector))
    
    if (nNA > i - 1){
      
      stop("NA problem!!!!")
      
    }
    
    SPEI_temp <- SPEI_temp[c(1:(nrow(SPEI_temp)-nNA)),c(1,2,3)]
    
    SPEI_temp$values <- vector[(nNA+1):length(vector)]
    
    SPEIX_monthly <- dcast(year ~ month, data = SPEI_temp, value.var = "values", na.rm = TRUE)
    
    row.names(SPEIX_monthly) <- SPEIX_monthly$year
    SPEIX_monthly$year <- NULL
    
    env_data = SPEIX_monthly
    
    ncol_response <- ncol(response)
    
    colnames_response <- colnames(response)
    
    env_data$yearABC <- row.names(env_data)
    response$yearABC <- row.names(response)
    
    temporal_data <- merge(response, env_data, by = "yearABC")
    
    response <- data.frame(temporal_data[, c(2:(1 + ncol_response))],
                           row.names = temporal_data$yearABC)
    colnames(response) <- colnames_response
    
    env_data <- data.frame(temporal_data[, c((1 + ncol_response + 1):
                                               ncol(temporal_data))],
                           row.names = temporal_data$yearABC)
    
    temporal_matrix <- matrix(NA, nrow = 1,
                              ncol = (ncol(env_data) - 1) + 1)
    # Here I create two additional temporal matrices, which will be used to store
    # lower and upper limits of bootstrap estimates
    temporal_matrix_lower <- temporal_matrix
    temporal_matrix_upper <- temporal_matrix
    
    for (j in 0: (ncol(env_data) - i)) {
      
      x <- env_data[, (j + 1)]
      
      x <- matrix(x, nrow = nrow(env_data), ncol = 1)
      
      if (boot == FALSE){
        temporal_correlation <- cor(response[, 1], x[, 1], method = cor_method)
        temporal_lower <- NA
        temporal_upper <- NA
      } else if (boot == TRUE){
        temp_df_boot <- cbind(response[, 1], x[, 1])
        
        temp_df_boot <- temp_df_boot[complete.cases(temp_df_boot), ]
        
        calc <- boot(temp_df_boot, boot_f, fun = "cor", cor.type = cor_method, R = boot_n)
        
        temporal_correlation <- colMeans(calc$t)[1]
        
        ci_int <- boot.ci(calc, conf = boot_conf_int, type = boot_ci_type)
        temporal_lower <- ci_int$norm[2]
        temporal_upper <- ci_int$norm[3]
      } else {
        print(paste0("boot should be TRUE or FALSE, instead it is ", boot))
      }
      
      temporal_matrix[1, j + 1] <- temporal_correlation
      temporal_matrix_lower[1, j + 1] <- temporal_lower
      temporal_matrix_upper[1, j + 1] <- temporal_upper
      
    }
    
    temporal_matrix_list[[b]] <- temporal_matrix
    temporal_matrix_lower_list[[b]] <- temporal_matrix_lower
    temporal_matrix_upper_list[[b]] <- temporal_matrix_upper
    b = b +1
    
  }
  
  temporal_matrix <- do.call(rbind, temporal_matrix_list)
  temporal_matrix_lower <- do.call(rbind, temporal_matrix_lower_list)
  temporal_matrix_upper <- do.call(rbind, temporal_matrix_upper_list)
  
  temporal_rownames <- as.vector(seq(from = 1, to = 9,
                                     by = 1))
  row.names(temporal_matrix) <- temporal_rownames
  row.names(temporal_matrix_lower) <- temporal_rownames
  row.names(temporal_matrix_upper) <- temporal_rownames
  
  temporal_colnames <- as.vector(seq(from = 1,
                                     to = ncol(temporal_matrix), by = 1))
  colnames(temporal_matrix) <- temporal_colnames
  colnames(temporal_matrix_lower) <- temporal_colnames
  colnames(temporal_matrix_upper) <- temporal_colnames
  
  #####################################################################
  
  if(is.finite(mean(temporal_matrix, na.rm = TRUE)) == FALSE){
    stop("All calculations are insignificant! Change the alpha argument!")
  }
  
  overall_max <- max(temporal_matrix, na.rm = TRUE)
  overall_min <- min(temporal_matrix, na.rm = TRUE)
  
  # absolute vales of overall_maximum and overall_minimum are compared and
  # one of the following two if functions is used
  # There are unimportant warnings produced:
  # no non-missing arguments to max; returning -Inf
  
  if ((abs(overall_max) > abs(overall_min)) == TRUE) {
    
    # maximum value is located. Row indeces are needed to query information
    # about the window width used to calculate the maximum. Column name is
    # needed to query the starting day.
    max_result <- suppressWarnings(which.max(apply(temporal_matrix,
                                                   MARGIN = 2, max,
                                                   na.rm = TRUE)))
    plot_column <- max_result
    max_index <- which.max(temporal_matrix[, names(max_result)])
    row_index <- row.names(temporal_matrix)[max_index]
  }
  
  if ((abs(overall_max) < abs(overall_min)) == TRUE) {
    
    min_result <- suppressWarnings(which.min(apply(temporal_matrix,
                                                   MARGIN = 2, min,
                                                   na.rm = TRUE)))
    plot_column <- min_result
    min_index <- which.min(temporal_matrix[, names(min_result)])
    row_index <- row.names(temporal_matrix)[min_index]
  }
  
  # Element 5
  # Here we create the fifth element of the final list: Analysed period in the
  # form of min(year) - max(year), e.g. 1950 - 2015
  min_env_data <- min(as.numeric(row.names(env_data)))
  min_response <- min(as.numeric(row.names(response)))
  
  max_env_data <- max(as.numeric(row.names(env_data)))
  max_response <- max(as.numeric(row.names(response)))
  
  min_together <- min(min_env_data, min_response)
  max_together <- min(max_env_data, max_response)
  
  
  analysed_period <- paste(as.character(min_together),
                           as.character(max_together),
                           sep = " - ")
  if (nchar(analysed_period) < 9) {
    analysed_period <- NA
  }
  
  # String for titles
  if (method == "cor"){
    title_string <- "Correlation Coefficients"
  } else if (method == "lm"){
    title_string <- "Linear Regression"
  } else if (method == "brnn"){
    title_string <- "ANN With Bayesian Regularization"
  } else (print("The selection of method is not correct"))
  
  final_list <- list(calculations = temporal_matrix, method = method,
                     metric = cor_method, analysed_period = analysed_period,
                     optimized_return = NA,
                     optimized_return_all = NA,
                     transfer_function = NA, temporal_stability = NA,
                     cross_validation = NA,
                     plot_heatmap = NA,
                     plot_extreme = NA,
                     plot_specific = NA,
                     PCA_output = NA,
                     type = "monthly",
                     reference_window = "start",
                     boot_lower = temporal_matrix_lower,
                     boot_upper = temporal_matrix_upper)
  
  #####################################################################
  
  # A) Extracting a matrix from a list and converting it into a data frame
  result_daily_response <- final_list
  object <- final_list
  
  type <- data.frame(object[[14]])
  
  result_daily_element1 <- data.frame(object[[1]])
  
  reference_window <- object[[15]]
  
  # To keep RCMD check happy:
  
  # With the following chunk, overall_maximum and overall_minimum values of
  # result_daily_element1 matrix are calculated.
  overall_max <- max(result_daily_element1, na.rm = TRUE)
  overall_min <- min(result_daily_element1, na.rm = TRUE)
  
  # absolute vales of overall_maximum and overall_minimum are compared and
  # one of the following two if functions is used
  # There are unimportant warnings produced:
  # no non-missing arguments to max; returning -Inf
  # Based on the answer on the StackOverlow site:
  # https://stackoverflow.com/questions/24282550/no-non-missing-arguments-warning-when-using-min-or-max-in-reshape2
  # Those Warnings could be easily ignored
  if ((abs(overall_max) > abs(overall_min)) == TRUE) {
    
    # maximum value is located. Row indeces are needed to query information
    # about the window width used to calculate the maximum. Column name is
    # needed to query the starting day.
    max_result <- suppressWarnings(which.max(apply(result_daily_element1,
                                                   MARGIN = 2, max, na.rm = TRUE)))
    plot_column <- max_result
    plot_column_source <- plot_column
    max_index <- which.max(result_daily_element1[, names(max_result)])
    row_index <- row.names(result_daily_element1)[max_index]
    temporal_vector <- unlist(result_daily_element1[max_index, ])
    temporal_vector <- data.frame(temporal_vector)
    calculated_metric <- round(max(temporal_vector, na.rm = TRUE), 3)
    
    lower_bound <- final_list$boot_lower[max_index, as.numeric(max_result)]
    upper_bound <- final_list$boot_upper[max_index, as.numeric(max_result)]
    
    # Here we remove missing values at the end of the temporal_vector.
    # It is important to remove missing values only at the end of the
    # temporal_vector!
    row_count <- nrow(temporal_vector)
    delete_rows <- 0
    while (is.na(temporal_vector[row_count, ] == TRUE)){
      delete_rows <- delete_rows + 1
      row_count <-  row_count - 1
    }
    # To check if the last row is a missing value
    if (is.na(temporal_vector[nrow(temporal_vector), ] == TRUE)) {
      temporal_vector <-  temporal_vector[-c((row_count + 1):(row_count +
                                                                delete_rows)), ]
    }
    temporal_vector <- data.frame(temporal_vector)
  }
  
  if ((abs(overall_max) < abs(overall_min)) == TRUE) {
    
    # minimum value is located. Row indeces are needed to query information
    # about the window width used to calculate the minimum. Column name is
    # needed to query the starting day.
    min_result <- suppressWarnings(which.min(apply(result_daily_element1,
                                                   MARGIN = 2, min, na.rm = TRUE)))
    plot_column <- min_result
    plot_column_source <- plot_column
    min_index <- which.min(result_daily_element1[, names(min_result)])
    row_index <- row.names(result_daily_element1)[min_index]
    temporal_vector <- unlist(result_daily_element1[min_index, ])
    temporal_vector <- data.frame(temporal_vector)
    calculated_metric <- round(min(temporal_vector, na.rm = TRUE), 3)
    
    lower_bound <- final_list$boot_lower[min_index, as.numeric(min_result)]
    upper_bound <- final_list$boot_upper[min_index, as.numeric(min_result)]
    # Here we remove missing values
    # We remove missing values at the end of the temporal_vector.
    # It is important to remove missing values only at the end of the
    # temporal_vector!
    
    row_count <- nrow(temporal_vector)
    delete_rows <- 0
    while (is.na(temporal_vector[row_count, ] == TRUE)){
      delete_rows <- delete_rows + 1
      row_count <-  row_count - 1
    }
    # To check if the last row is a missing value
    if (is.na(temporal_vector[nrow(temporal_vector), ] == TRUE)) {
      temporal_vector <-  temporal_vector[-c((row_count + 1):(row_count +
                                                                delete_rows)), ]
    }
    temporal_vector <- data.frame(temporal_vector)
  }
  
  # In case of previous_year == TRUE, we calculate the day of a year
  # (plot_column), considering 366 days of previous year.
  if (nrow(temporal_vector) > 366 & plot_column > 366) {
    previous_year = TRUE
    plot_column_extra <- plot_column %% 366
  } else {
    previous_year = FALSE
    plot_column_extra <- plot_column
  }
  
  if (nrow(temporal_vector) > 366) {
    previous_year <- TRUE
  } else {
    previous_year <- FALSE
  }
  
  # Here we define a data frame of dates and corresponing day of year (doi). Later
  # this dataframe will be used to describe tht optimal sequence od days
  doy <- seq(1:730)
  date <- seq(as.Date('2013-01-01'),as.Date('2014-12-31'), by = "+1 day")
  # date[366] <- as.Date('2015-12-31')
  date <- format(date, "%b %d")
  date_codes <- data.frame(doy = doy, date = date)
  
  # Here, there is a special check if optimal window width is divisible by 2 or not.
  if (as.numeric(row_index)%%2 == 0){
    adjustment_1 = 0
    adjustment_2 = 1
  } else {
    adjustment_1 = 1
    adjustment_2 = 2
  }
  
  if (reference_window == "start"){
    Optimal_string <- paste(as.character(date_codes[plot_column_extra, 2]),"-",
                            as.character(date_codes[plot_column_extra + as.numeric(row_index) - 1, 2]))
  } else if (reference_window == "end") {
    Optimal_string <- paste(as.character(date_codes[plot_column_extra - as.numeric(row_index) + 1, 2]),"-",
                            as.character(date_codes[plot_column_extra, 2]))
  } else if (reference_window == "middle") {
    Optimal_string <- paste(as.character(date_codes[(round2((plot_column_extra - as.numeric(row_index)/2)) - adjustment_1), 2]),"-",
                            as.character(date_codes[(round2((plot_column_extra + as.numeric(row_index)/2)) - adjustment_2), 2]))
  }
  
  # Here we define titles. They differ importantly among methods and arguments
  # in the final output list from daily_response() function
  if (result_daily_response[[2]] == "cor"){
    y_lab <- NA
  } else if (result_daily_response[[2]] == "pcor"){
    y_lab <- NA
  } else if (result_daily_response[[3]] == "r.squared"){
    y_lab <- "Explained Variance"
  } else if (result_daily_response[[3]] == "adj.r.squared"){
    y_lab <- "Adjusted Explained Variance"
  }
  
  if (reference_window == 'start' &&  plot_column > 366 && nrow(temporal_vector) > 366){
    reference_string <- paste0("Starting Day of Optimal Window Width: Day ",
                               plot_column_extra, " of Current Year")}
  
  if (reference_window == 'start' &&  plot_column <= 366 && nrow(temporal_vector) > 366){
    reference_string <- paste0("Starting Day of Optimal Window Width: Day ",
                               plot_column_extra, " of Previous Year")}
  
  if (reference_window == 'start' &&  plot_column <=  366 && nrow(temporal_vector) <=  366){
    reference_string <- paste0("Starting Day of Optimal Window Width: Day ",
                               plot_column_extra)}
  
  
  if (reference_window == 'end' &&  plot_column > 366 && nrow(temporal_vector) > 366){
    reference_string <- paste0("Ending Day of Optimal Window Width: Day ",
                               plot_column_extra, " of Current Year")}
  
  if (reference_window == 'end' &&  plot_column  <= 366 && nrow(temporal_vector) > 366){
    reference_string <- paste0("Ending Day of Optimal Window Width: Day ",
                               plot_column_extra, " of Previous Year")}
  
  if (reference_window == 'end' &&  plot_column  <=  366 && nrow(temporal_vector) <=  366){
    reference_string <- paste0("Ending Day of Optimal Window Width: Day ",
                               plot_column_extra)}
  
  
  if (reference_window == 'middle' &&  plot_column > 366 && nrow(temporal_vector) > 366){
    reference_string <- paste0("Middle Day of Optimal Window Width: Day ",
                               plot_column_extra, " of Current Year")}
  
  if (reference_window == 'middle' &&  plot_column  <= 366 && nrow(temporal_vector) > 366){
    reference_string <- paste0("Middle Day of Optimal Window Width: Day ",
                               plot_column_extra, " of Previous Year")}
  
  if (reference_window == 'middle' &&  plot_column  <=  366 && nrow(temporal_vector) <=  366){
    reference_string <- paste0("Middle Day of Optimal Window Width: Day ",
                               plot_column_extra)}
  
  optimal_window_string <- paste0("Optimal Window Width: ", as.numeric(row_index),
                                  " Days")
  
  optimal_calculation <- paste0("The Highest ", y_lab,": " , calculated_metric)
  
  period_string <- paste0("Analysed Period: ", result_daily_response[[4]])
  
  if (result_daily_response[[2]] == 'cor'){
    method_string <- paste0("Correlation Coefficient (", result_daily_response[[3]], ")")
    
  } else if (result_daily_response[[2]] == 'pcor'){
    method_string <- paste0("Partial Correlation Coefficient (", result_daily_response[[3]], ")")
    
  } else if (result_daily_response[[2]] == 'lm'){
    method_string <- paste0("Linear Regression")
  } else if (result_daily_response[[2]] == 'brnn'){
    method_string <- paste0("ANN with Bayesian Regularization")
  }
  
  if (type == "monthly"){
    
    # Plural or singular?
    if (as.numeric(row_index) == 1){
      month_string <- " Month"
      
    } else {
      month_string <- " Months"
    }
    
    # In case of previous_year == TRUE, we calculate the day of a year
    # (plot_column), considering 366 days of previous year.
    
    if (ncol(result_daily_response[[1]]) == 24 & plot_column > 12) {
      plot_column_extra <- plot_column %% 12
    } else {
      plot_column_extra <- plot_column
    }
    
    
    if (ncol(result_daily_response[[1]]) == 24) {
      previous_year <- TRUE
    } else {
      previous_year <- FALSE
    }
    
    if (reference_window == 'start' &&  plot_column > 12 && ncol(result_daily_response[[1]]) == 24){
      reference_string <- paste0("Starting Month of Optimal Window Width: Month ",
                                 plot_column_extra, " of Current Year")}
    
    if (reference_window == 'start' &&  plot_column <= 12 && ncol(result_daily_response[[1]]) == 24){
      reference_string <- paste0("Starting Month of Optimal Window Width: Month ",
                                 plot_column_extra, " of Previous Year")}
    
    if (reference_window == 'start' &&  plot_column <=  12 && ncol(result_daily_response[[1]]) <=  12){
      reference_string <- paste0("Starting Month of Optimal Window Width: Month ",
                                 plot_column_extra)}
    
    optimal_window_string <- paste0("Optimal Window Width: ", as.numeric(row_index),
                                    month_string)
    
    # Here we define a data frame of months. Later
    # this dataframe will be used to describe tht optimal sequence od days
    
    if (ncol(result_daily_response[[1]]) == 24){
      date_codes <- c("Jan*", "Feb*", "Mar*", "Apr*", "May*", "Jun*", "Jul*", "Aug*", "Sep*", "Oct*", "Nov*", "Dec*",
                      "Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")
      
    } else if (ncol(result_daily_response[[1]]) == 12){
      
      date_codes <- c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")
      
    }
    
    if (reference_window == "start"){
      Optimal_string <- paste(as.character(date_codes[plot_column_source]),"-",
                              as.character(date_codes[plot_column_source + as.numeric(row_index) - 1]))
    } else if (reference_window == "end") {
      Optimal_string <- paste(as.character(date_codes[plot_column_source - as.numeric(row_index) + 1]),"-",
                              as.character(date_codes[plot_column_source]))
    } else if (reference_window == "middle") {
      Optimal_string <- paste(as.character(date_codes[(round2((plot_column_source - as.numeric(row_index)/2)) - adjustment_1)]),"-",
                              as.character(date_codes[(round2((plot_column_source + as.numeric(row_index)/2)) - adjustment_2)]))
    }
    
    if (as.numeric(row_index == 1)){
      Optimal_string <- substr(Optimal_string, 1, nchar(Optimal_string)-6)
    }
    
  }
  
  output_df <-  data.frame(Variable = c("approach",
                                        "method",
                                        "metric",
                                        "analysed_years",
                                        "maximal_calculated_metric",
                                        "lower_ci",
                                        "upper_ci",
                                        "reference_window",
                                        "analysed_previous_year",
                                        "optimal_time_window",
                                        "optimal_time_window_length"),
                           
                           Value = c(result_daily_response[[14]],
                                     method_string,
                                     y_lab,
                                     result_daily_response[[4]],
                                     calculated_metric,
                                     round(lower_bound, 3),
                                     round(upper_bound, 3),
                                     reference_string,
                                     previous_year,
                                     Optimal_string,
                                     as.numeric(row_index)))
  
############################################# Extract results 
  SPEIa_monthly <- step2_SPEI(output_df)
  colnames(SPEIa_monthly) <- paste0(colnames(SPEIa_monthly), "_monthly")
  
#############################################

final_results <- cbind(SPEIa_daily, SPEIa_monthly)

final_results$variable <- "SPEI"
  
final_results <- tibble::add_column(final_results, key_code = meta_file[m,"file_code"], .before = 1)

write.table(final_results, gsub(".crn","",paste0("C:/Users/Dendro/Desktop/GDS/calculations_yesBOOT/SPEI_daily_", meta_file[m,"file_code"], ".txt")))

print(as.character(meta_file[m,"file_code"]))

}
