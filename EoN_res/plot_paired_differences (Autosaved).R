# Creation of Table 6 trajectories and quantiles -------------------------------

# Part 1: Get quantile information and histogram of paired differences ---------

# Get output folder name -------------------------------------------------------
output_folder <- "~/Philly_Covid/EoN_res/output_paired_differences/"

# Parameters of Table 6 --------------------------------------------------------
Res <- c("1_5")
fusings <- c(100, 200, 50)
for (Re in Res) {
  for (fusing in fusings) {
    # Get input folder name ----------------------------------------------------
    input_folder <- paste0("~/Philly_Covid/EoN_res/Re_", Re, "_fuse_", fusing, "/")
    print(paste0("Re: ", Re))
    print(paste0("fusing: ", fusing))
    paired_differences <- vector()
    nsims <- length(list.files(input_folder))
    for (nsim in 1:nsims) {
      # Plot paired difference -------------------------------------------------
      paired_difference_name <- paste0(input_folder, "batch", nsim, "/plateau_paired_differences_", gsub("_", ".", Re), "_", fusing, ".csv")
      paired_difference <- read.csv(paired_difference_name, stringsAsFactors = F, header=F)
      paired_differences <- c(paired_differences, unname(unlist(paired_difference)))
    }
    # Write paired difference results ------------------------------------------
    o1 <- paste0("Re_", Re, "_fusing_", fusing)
    o2 <- paste0("nsims: ", nsims)
    o3 <- paste0("quantiles 0,.025,.25,.5,.75,.975, 1 :")
    o4 <- paste0(quantile(paired_differences, c(0,.025,.25,.5,.75,.975,1)))
    outputResult<-list(o1, o2, o3, o4)
    filename <- file.path(paste0(output_folder, "Re_", Re, "_fusing_", fusing, ".txt"))
    capture.output(outputResult, file = filename)
    
    # Plot paired difference histogram -----------------------------------------
    pdf(paste0(output_folder, "Re_", Re, "_fusing_", fusing, ".pdf"), width=7, height=5)
    hist(paired_differences, breaks=10)
    dev.off()
  }
} 

# Part 2: Get plots of trajectores ---------------------------------------------

# Get output folder name -------------------------------------------------------
output_folder <- "~/Philly_Covid/EoN_res/output_trajectories/"

N <- 100000
sd_date <- "2020/03/23"
easing_date <- "2020/06/05"
second_easing_date <- "2020/08/01"
evictions_date <- "2020/09/01"
second_sd_date <- "2020/11/1"
sd_to_easing <- as.numeric(as.Date(easing_date, format="%Y/%m/%d") - as.Date(sd_date, format="%Y/%m/%d"))
easing_to_second_easing <- as.numeric(as.Date(second_easing_date, format="%Y/%m/%d") - as.Date(easing_date, format="%Y/%m/%d"))
second_easing_to_evictions <- as.numeric(as.Date(evictions_date, format="%Y/%m/%d") - as.Date(second_easing_date, format="%Y/%m/%d"))
evictions_to_second_sd <- as.numeric(as.Date(second_sd_date, format="%Y/%m/%d") - as.Date(evictions_date, format="%Y/%m/%d"))
for (Re in Res) {
  for (fusing in fusings) {
    # Get input folder name ----------------------------------------------------
    input_folder <- paste0("~/Philly_Covid/EoN_res/Re_", Re, "_fuse_", fusing, "/")
    print(paste0("Re: ", Re))
    print(paste0("fusing: ", fusing))
    align_point_FSD <- NA
    align_point_EQ <- NA
    start_of_evictions_day <- NA
    nsims <- length(list.files(input_folder))
    for (nsim in 1:nsims) {
      # Plot aligned trajectory ------------------------------------------------
      FSD_trajectory <- read.csv(paste0(input_folder, "batch", nsim, "/csvs/0_FSD.csv"), stringsAsFactors = F, header = F)
      EQ_trajectory <- read.csv(paste0(input_folder, "batch", nsim, "/csvs/0_EQ.csv"), stringsAsFactors = F, header = F)
      
      if (nsim == 1) {
        align_point_FSD <- FSD_trajectory$V1[min(which(((FSD_trajectory$V4 + FSD_trajectory$V3) / N) >= 0.02))]
        align_point_EQ <- EQ_trajectory$V1[min(which(((EQ_trajectory$V4 + EQ_trajectory$V3) / N) >= 0.02))]
        png(paste0(output_folder, "Re_", Re, "_fusing_", fusing, ".png"), width=1200, height=800)
        plot(FSD_trajectory$V1, FSD_trajectory$V4 / N, col='black', pch='.', ylim=c(0, 0.02), xlab="time (days)", ylab="% Infectious")
        
        start_of_evictions_day <- align_point_FSD + sd_to_easing + easing_to_second_easing + second_easing_to_evictions
        
        min_dex <- min(which(EQ_trajectory$V1 >= start_of_evictions_day))
        lines(EQ_trajectory$V1[min_dex:nrow(EQ_trajectory)], EQ_trajectory$V4[min_dex:nrow(EQ_trajectory)] / N, col='red', pch='.')
      } else {
        align_point_FSD_focal <- FSD_trajectory$V1[min(which(((FSD_trajectory$V4 + FSD_trajectory$V3) / N) >= 0.02))]
        align_point_EQ_focal <- EQ_trajectory$V1[min(which(((EQ_trajectory$V4 + EQ_trajectory$V3) / N) >= 0.02))]
        to_move_FSD <- align_point_FSD - align_point_FSD_focal
        to_move_EQ <- align_point_EQ - align_point_EQ_focal
        stopifnot(to_move_FSD == to_move_EQ)
        new_time_FSD <- FSD_trajectory$V1 + to_move_FSD
        new_time_EQ <- EQ_trajectory$V1 + to_move_EQ
        
        min_dex <- min(which(new_time_EQ >= start_of_evictions_day))
        lines(new_time_FSD, FSD_trajectory$V4 / N, col='black', lwd=.1)
        lines(new_time_EQ[min_dex:nrow(EQ_trajectory)], EQ_trajectory$V4[min_dex:nrow(EQ_trajectory)] / N, col='red', lwd=.1)
      }
    }
    abline(v=align_point_FSD, col='blue')
    abline(v=align_point_FSD + sd_to_easing, col='blue', lty="dashed")
    abline(v=align_point_FSD + sd_to_easing + easing_to_second_easing, col='blue', lty="dashed")
    abline(v=align_point_FSD + sd_to_easing + easing_to_second_easing + second_easing_to_evictions, col='blue', lty="dashed")
    abline(v=align_point_FSD + sd_to_easing + easing_to_second_easing + second_easing_to_evictions + evictions_to_second_sd, col='blue')
    dev.off()
  }
}
