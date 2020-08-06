# Creation of Table 6 trajectories and quantiles -------------------------------

# Get output folder name -------------------------------------------------------
output_folder <- "~/Philly_Covid/EoN_res/output_paired_differences/"

# Parameters of Table 6 --------------------------------------------------------
Res <- c("1_5")
fusings <- c(100, 200, 50)
nsims <- 30

for (Re in Res) {
  for (fusing in fusings) {
    # Get input folder name ----------------------------------------------------
    input_folder <- paste0("~/Philly_Covid/EoN_res/Re_", Re, "_fuse_", fusing, "/")
    print(paste0("Re: ", Re))
    print(paste0("fusing: ", fusing))
    align_point_FSD <- NA
    align_point_EQ <- NA
    paired_differences <- vector()
    for (nsim in 1:nsims) {
      # Plot paired difference -------------------------------------------------
      paired_difference_name <- paste0(input_folder, "batch", nsim, "/plateau_paired_differences_", gsub("_", ".", Re), "_", fusing, ".csv")
      paired_difference <- read.csv(paired_difference_name, stringsAsFactors = F, header=F)
      paired_differences <- c(paired_differences, unname(unlist(paired_difference)))
    
    # Write paired difference results ------------------------------------------
    o1 <- paste0("Re_", Re, "_fusing_", fusing)
    o2 <- paste0("nsims: ", nsims)
    o3 <- paste0("quantiles:")
    o4 <- paste0(quantile(paired_differences))
    outputResult<-list(o1, o2, o3, o4)
    filename <- file.path(paste0(output_folder, "Re_", Re, "_fusing_", fusing, ".txt"))
    capture.output(outputResult, file = filename)
    
    # Plot paired difference histogram -----------------------------------------
    pdf(paste0(output_folder, "Re_", Re, "_fusing_", fusing, ".pdf"), width=5, height=5)
    hist(paired_differences, breaks=10)
    dev.off()
    }
  }
} 


for (Re in Res) {
  for (fusing in fusings) {
    # Plot aligned trajectory ------------------------------------------------
    FSD_trajectory <- read.csv(paste0(input_folder, "batch", nsim, "/csvs/0_FSD.csv"), stringsAsFactors = F, header = F)
    EQ_trajectory <- read.csv(paste0(input_folder, "batch", nsim, "/csvs/0_EQ.csv"), stringsAsFactors = F, header = F)
    
    if (nsim == 1) {
      align_point_FSD <- min(which(FSD_trajectory) > 0.01)
      align_point_EQ <- min(which(EQ_trajectory) > 0.01)
      pdf(paste0(output_folder, "/Re_", Re, "_fusing_", fusing, ".pdf"))
      plot(1:length(FSD_trajectory), FSD_trajectory, col='black')
      lines(1:length(EQ_trajectory), EQ_trajectory, col='red')
    } else {
      align_point_FSD_focal <- min(which(FSD_trajectory) > 0.01)
      align_point_EQ_focal <- min(which(EQ_trajectory) > 0.01)
      to_move_FSD <- align_point_FSD_focal - align_point_FSD
      to_move_EQ <- align_point_EQ_focal - align_point_EQ
      FSD_trajectory$time <- FSD_trajectory$time + to_move_FSD
      EQ_trajectory$time <- EQ_trajectory$time + to_move_EQ
      
      lines(FSD_trajectory$time, FSD_trajectory, col='black')
      lines(EQ_trajectory$time, EQ_trajectory, col='red')
    }
    dev.off()
  }
}
