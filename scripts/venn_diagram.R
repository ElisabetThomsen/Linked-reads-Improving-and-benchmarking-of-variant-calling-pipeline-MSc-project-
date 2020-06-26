# https://rstudio-pubs-static.s3.amazonaws.com/13301_6641d73cfac741a59c0a851feb99e98b.html

#install.packages('VennDiagram')
library(VennDiagram)

# How to get the number to put in see https://github.com/ElisabetThomsen/MSc_thesis#venn-diagram

# Only target INDELs - NA12878
grid.newpage()
draw.triple.venn(area1 = 22902, area2 = 16481, area3 = 10040,
                 n12 = 11290, n23 = 7803, n13 = 8766, 
                 n123 = 7450,
                 category = c("LinkSeq", "Long Ranger", "Giab Truth"),
                 lty = "blank", 
                 fill = c("blue", "green", "red"),
                 col = "transparent",
                 alpha = 0.30,
                 cex = 0.7,
                 cat.cex = 0.7,
                 print.mode=c("raw","percent"))

# Only target SNPs - NA12878
grid.newpage()
draw.triple.venn(area1 = 104830, area2 = 118532, area3 = 73881,
                 n12 = 96448, n23 = 73533, n13 = 72900, 
                 n123 = 72829,
                 category = c("LinkSeq", "Long Ranger", "Giab Truth"),
                 lty = "blank", 
                 fill = c("blue", "green", "red"),
                 col = "transparent",
                 alpha = 0.30,
                 cex = 0.7,
                 cat.cex = 0.7,
                 print.mode=c("raw","percent"))