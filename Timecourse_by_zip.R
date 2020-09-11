 setwd("/Users/michaellevy/covid19-philadelphia/")
 library(tidyverse)
 library(lubridate)
 library(EpiEstim)

# Stack daily files within folder and keep distinct records
build_historical_dataset <- function(data_folder) {
  list.files(path = data_folder, full.names = TRUE) %>%
    map_dfr(read_csv) %>%
    distinct() # duplicates occur when there is an extract but no data update
}
cases_by_date <- build_historical_dataset("cases_by_date")

cases_by_zipcode <- build_historical_dataset("cases_by_zipcode")

d<-cases_by_zipcode


plot_zip<-function(zip){
	dzip<-d[d$zip_code==zip,]
	startdate<-(as.Date(min(dzip$etl_timestamp)))
	enddate<-(as.Date(max(dzip$etl_timestamp)))
     plot(dzip$etl_timestamp,dzip$POS, ylim=c(0,2000), type='l')
}

add_zip<-function(zip, COL){
	 dzip<-d[d$zip_code==zip,]
     lines(dzip$etl_timestamp,dzip$POS, col=COL)
}

plot_zip(19143)
zips<-unique(d$zip_code)
for(i in zips){
	dzip<-d[d$zip_code==i,]
	lines(dzip$etl_timestamp,dzip$POS, col=i)
	text(x=last(dzip$etl_timestamp), y=last(dzip$POS), label = i, cex=.5, col=i, ltp=3)
}




###below was written by someone else--we need to stratify by zipcode
# Daily case count by test result date
incidence_data <- list.files(path = "cases_by_date", full.names = TRUE) %>% 
  last %>% 
  read_csv(col_types = cols(result_date = col_character(),
                            etl_timestamp = col_skip(),
                            negative = col_integer(),
                            positive = col_integer())) %>% 
  filter(!is.na(positive) & (date(result_date) >= date("2020-03-16"))) %>%
  mutate(dates = date(result_date)) %>%
  arrange(dates) %>%
  mutate(positivity_rate = positive / (positive + negative)) %>%
  mutate(positivity_rate_last = lag(positivity_rate, default = first(positivity_rate))) %>%
  
  
  select(dates, positive, negative, positivity_rate, positivity_rate_last) %>%
  filter(dates <= last(dates) - 3) # remove last 3 days considering lag in test results

# Plot incidence and effective reproductive number over time
# Serial interval mean and std estimates from: https://www.dhs.gov/publication/st-master-question-list-covid-19
res_parametric_si <- estimate_R(incidence_data %>% 
                                  select(dates, I = positive),
                                method="parametric_si",
                                config = make_config(list(mean_si = 5.29, std_si = 5.32)))
plot(res_parametric_si, legend = FALSE)
