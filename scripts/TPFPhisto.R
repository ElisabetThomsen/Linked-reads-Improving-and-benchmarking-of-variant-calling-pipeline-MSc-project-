# Making density plots of TP and FP
# http://www.sthda.com/english/wiki/ggplot2-density-plot-quick-start-guide-r-software-and-data-visualization

### Import data
TPdata <- read.table("giab_ann.txt.TP", sep="\t", header=TRUE)
FPdata <- read.table("giab_ann.txt.FP", sep="\t", header=TRUE)


### Prepare the data
set.seed(1234)
df <- data.frame(
                   Type=factor(rep(c("FP", "TP"), c(length(FPdata$QUAL),length(TPdata$QUAL)))),
                     Quality=round(c(FPdata$QUAL,TPdata$QUAL)),
                     Depth=round(c(FPdata$DP,TPdata$DP)),
                     QD=round(c(FPdata$QD,TPdata$QD)),
                     FS=round(c(FPdata$FS,TPdata$FS)),
                     SOR=c(FPdata$SOR,TPdata$SOR),
                     MQ=round(c(FPdata$MQ,TPdata$MQ)),
                     MQRankSum=round(c(FPdata$MQRankSum,TPdata$MQRankSum)),
                     ReadPosRankSum=round(c(FPdata$ReadPosRankSum,TPdata$ReadPosRankSum))
                     )

### Calculate mean of each group
library(plyr)
qual_mu <- ddply(df, "Type", summarise, grp.mean=mean(Quality))

### Basic density plots
library(ggplot2)
# Basic density
# Color and fill by groups
qual_p <- ggplot(df, aes(x=Quality, fill=Type)) +
      # Make fill transparent (alpha)
      geom_histogram(alpha=0.4, binwidth = 1, position = "identity") +
        # Add mean line
        geom_vline(xintercept=30, color="purple", linetype="dashed") +
  # Add titles
  labs(title="Quality Score",x="Quality Score (QS)", y = "# variants with x QS") +
    # Center the title
    theme(plot.title = element_text(hjust = 0.5)) +
      scale_x_continuous(
                           labels = scales::number_format(accuracy = 1)) +
  scale_y_continuous(
                         labels = scales::number_format(accuracy = 0.001)) +
  xlim(0, 1000)


  ### Plot depth
  db_mu <- ddply(df, "Type", summarise, grp.mean=mean(Depth))

  db_p <- ggplot(df, aes(x=Depth, fill=Type)) +
        # Make fill transparent (alpha)
        geom_histogram(alpha=0.4, binwidth = 1, position = "identity") +
          # Add mean line
          geom_vline(xintercept=10, color="purple", linetype="dashed") +
  # Add titles
  labs(title="Depth of Coverage",x="Depth of Coverage", y = "# variants with x Coverage") +
    # Center the title
    theme(plot.title = element_text(hjust = 0.5)) +
      scale_x_continuous(
                           labels = scales::number_format(accuracy = 1)) +
  scale_y_continuous(
                         labels = scales::number_format(accuracy = 0.001)) +
  xlim(0, 400)


  ### Plot QD
  qd_mu <- ddply(df, "Type", summarise, grp.mean=mean(QD))

  qd_p <- ggplot(df, aes(x=QD, fill=Type)) +
        # Make fill transparent (alpha)
        geom_histogram(alpha=0.4, binwidth = 1, position = "identity") +
          # Add mean line
          geom_vline(xintercept=2, color="purple", linetype="dashed") +
  # Add titles
  labs(title="Quality by depth",x="Quality by depth (QD)", y = "# variants with x QD") +
    # Center the title
    theme(plot.title = element_text(hjust = 0.5)) +
      scale_x_continuous(
  labels = scales::number_format(accuracy = 1)) +
  scale_y_continuous(
                         labels = scales::number_format(accuracy = 0.001)) +
  xlim(0, 40)


  ### Plot FS
  fs_mu <- ddply(df, "Type", summarise, grp.mean=mean(FS))

  fs_p <- ggplot(df, aes(x=FS, fill=Type)) +
        # Make fill transparent (alpha)
        geom_histogram(alpha=0.4, binwidth = 1, position = "identity") +
          # Add mean line
          geom_vline(xintercept=60, color="purple", linetype="dashed") +
  # Add titles
  labs(title="FisherStrand",x="FisherStrand (FS)", y = "# variants with x FS") +
    # Center the title
    theme(plot.title = element_text(hjust = 0.5)) +
      scale_x_continuous(
  labels = scales::number_format(accuracy = 1)) +
  scale_y_continuous(
                         labels = scales::number_format(accuracy = 0.001)) +
  xlim(0, 65)


### Plot SOR
  sor_mu <- ddply(df, "Type", summarise, grp.mean=mean(SOR))

  sor_p <- ggplot(df, aes(x=SOR, fill=Type)) +
        # Make fill transparent (alpha)
      geom_histogram(alpha=0.4, binwidth = 0.1, position = "identity") +
          # Add mean line
      geom_vline(xintercept=3, color="purple", linetype="dashed") +
  # Add titles
      labs(title="StrandOddsRatio",x="StrandOddsRatio (SOR)", y = "# variants with x SOR") +
    # Center the title
    theme(plot.title = element_text(hjust = 0.5)) +
  scale_x_continuous(labels = scales::number_format(accuracy = 1)) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.001)) +
  xlim(0, 5)


### Plot MQ
  mq_mu <- ddply(df, "Type", summarise, grp.mean=mean(MQ))

  mq_p <- ggplot(df, aes(x=MQ, fill=Type)) +
      # Make fill transparent (alpha)
      geom_histogram(alpha=0.4, binwidth = 1, position = "identity") +
        # Add mean line
        geom_vline(xintercept=40, color="purple", linetype="dashed") +
  # Add titles
  labs(title="RMS Mapping Quality",x="RMS Mapping Quality (MQ)", y = "# variants with x MQ") +
    # Center the title
    theme(plot.title = element_text(hjust = 0.5)) +
      scale_x_continuous(
                           labels = scales::number_format(accuracy = 1)) +
  scale_y_continuous(
                         labels = scales::number_format(accuracy = 0.001)) +
  xlim(30, 90)

### Plot MQRankSum
  mqrs_mu <- ddply(df, "Type", summarise, grp.mean=mean(MQRankSum))

  mqrs_p <- ggplot(df, aes(x=MQRankSum, fill=Type)) +
      # Make fill transparent (alpha)
      geom_histogram(alpha=0.4, binwidth = 1, position = "identity") +
        # Add mean line
        geom_vline(xintercept=-12.5, color="purple", linetype="dashed") +
  # Add titles
  labs(title="MappingQualityRankSumTest",x="MQRankSum", y = "# variants with x MQRankSum") +
    # Center the title
    theme(plot.title = element_text(hjust = 0.5)) +
      scale_x_continuous(
  labels = scales::number_format(accuracy = 1)) +
  scale_y_continuous(
                         labels = scales::number_format(accuracy = 0.001)) +
  xlim(-20, 20)

### Plot ReadPosRankSum
  rprs_mu <- ddply(df, "Type", summarise, grp.mean=mean(ReadPosRankSum))

  rprs_p <- ggplot(df, aes(x=ReadPosRankSum, fill=Type)) +
      # Make fill transparent (alpha)
      geom_histogram(alpha=0.4, binwidth = 1, position = "identity") +
        # Add mean line
        geom_vline(xintercept=-8, color="purple", linetype="dashed") +
  # Add titles
  labs(title="ReadPosRankSumTest",x="ReadPosRankSum", y = "# variants with x ReadPosRankSum") +
    # Center the title
    theme(plot.title = element_text(hjust = 0.5)) +
      scale_x_continuous(
  labels = scales::number_format(accuracy = 1)) +
  scale_y_continuous(
                         labels = scales::number_format(accuracy = 0.001)) +
  xlim(-20, 20)


pdf("TPFPhisto.pdf")
qual_p
db_p
qd_p
fs_p
sor_p
mq_p
mqrs_p
rprs_p

  # Print plots nicely
  #library(gridExtra)
  # Arranges the plots
  #arr_plots <- grid.arrange(qual_p, db_p, nrow = 1)
  # Save the plots
  #dev.copy(pdf,"TPFPplots.pdf", width = 8, height = 4)
  dev.off()


  ### How to customize
  # Use custom color palettes
  #p+scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9"))
  #+scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9"))
  # Change the legend position
  #p + theme(legend.position="top")
  #p + theme(legend.position="bottom")
  #p + theme(legend.position="none") # Remove legend
