library(tidyverse)
library(RColorBrewer)
library(fmsb)

# Change working directory and load data
setwd("XXX")
phl <- read.csv("PhillyCovid19_4.2.2020.csv")
zip <- readRDS("Data/zipcodes_phila_df.rds")
philly <- readRDS("Data/phillybd.rds")
adi <- readRDS("Data/ADI_zcta.rds")
adi$ZIP <- adi$GEOID

# Prep zip-code level data frame for covid test data
zip_cases <- phl %>%
  mutate(ZIP = as.character(ZIP)) %>%
  group_by(ZIP) %>%
  summarise(cases = sum(Result == "POS"), total = n(),
            percent.pos = cases/total*100)

# Merge covid test and ADI data to zipcode shapefile
zip <- left_join(zip, zip_cases)
zip <- left_join(zip, adi[,c("ZIP", "ADI")])
zip$ADI_percentile <- zip$ADI
zip$ADI_percentile[! is.na(zip$ADI)] <- percentile(na.omit(zip$ADI))

my_theme <- function(){
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.title =element_blank(),
        axis.text =element_blank(),
        axis.ticks =element_blank(),
        legend.title=element_text(size=24),
        legend.text = element_text(size = 22),
        plot.title = element_text(size = 26, face = "bold"))
}  

# Make maps
ggplot(data=zip, aes(x=long, y=lat, group=group)) +
  geom_polygon(aes(fill=cases)) +
  ggtitle("Number of COVID-19 cases by zipcode") +
  scale_fill_gradientn(name = "Positive \ntests", colors = brewer.pal(n = 9, name = "Reds")) +
  geom_path(data = philly, aes(x=long, y=lat, group=group), 
            col = "black", lwd=0.9) +
  my_theme()

ggplot(data=zip, aes(x=long, y=lat, group=group)) +
  geom_polygon(aes(fill=percent.pos)) +
  ggtitle("Percent positive by zipcode") +
  scale_fill_gradientn(name = "Percent \npositive (%)", colors = brewer.pal(n = 9, name = "BuPu")) +
  geom_path(data = philly, aes(x=long, y=lat, group=group), 
            col = "black", lwd=0.9) +
  my_theme()

ggplot(data=zip, aes(x=long, y=lat, group=group)) +
  geom_polygon(aes(fill=total)) +
  ggtitle("Tests administered by zipcode") +
  scale_fill_gradientn(name = "Total \ntests", colors = brewer.pal(n = 9, name = "Blues")) +
  geom_path(data = philly, aes(x=long, y=lat, group=group), 
            col = "black", lwd=0.9) +
  my_theme()

ggplot(data=zip, aes(x=long, y=lat, group=group)) +
  geom_polygon(aes(fill=ADI_percentile)) +
  ggtitle("Neighborhood disadvantage of Philadelphia zip codes") +
  scale_fill_gradientn(name = "Neighborhood \ndisadvantage", colors = brewer.pal(n = 9, name = "Oranges")) +
  geom_path(data = philly, aes(x=long, y=lat, group=group), 
            col = "black", lwd=0.9) +
  my_theme()

