library(ggplot2)
library(tidyr)
#library(grid)
#library(gridExtra)

day <- format(as.Date(date(), format="%a %b %d %H:%M:%S %Y"), format='%y%m%d')

dat <- read.csv("Data/niche_PC_4axes190813.csv", stringsAsFactors = FALSE)
dat$bin <- - dat$bin
pc1 <- dat[,1:4]
pc_wide <- spread(pc1, stat, PC1)

# 3 species have no more than 2 consecutive observations.
spp <- sort(unique(dat$sp))
rare <- c("Globoquadrina conglomerata","Truncorotalia crassula","Truncorotalia truncatulinoides")
spp <- setdiff(spp, rare)

# There are too many species to plot at once; put them on 2 figures.
set1 <- spp[1:12]
set2 <- spp[13:24]
set1_rows <- pc_wide$sp %in% set1
p1 <- pc_wide[set1_rows,]
set2_rows <- pc_wide$sp %in% set2
p2 <- pc_wide[set2_rows,]

lw <- 0.5
min_y <- -3
max_y <- 8 # max(dat$PC1, na.rm=TRUE)*1.1

# Names are too long to print. Omit genus name.
get_sp_nm <- function(txt){
  splt <- strsplit(txt, ' ')[[1]]
  tail(splt,1)
}
spp_short <- sapply(spp, get_sp_nm)

empty <- ggplot() + 
  theme_bw() + 
  # labs(title = '') + theme(plot.title = element_blank()) +
  scale_x_continuous(name = 'ka', expand=c(0,0), limits=c(-800,0),
                     breaks = seq(-800, 0, by= 32),
                     labels = paste(-seq(-800, 0, by=32)) ) +
  scale_y_continuous(name= 'PC1', expand=c(0,0), limits=c(min_y,max_y),
                     breaks=c(0,4,8)) +
  theme(axis.text.x = element_text(angle=90, vjust=0.5),
        legend.text = element_text(size=8)
  )

plot1 <- empty +
    geom_ribbon(data=p1, aes(x=bin, ymin=min, ymax=max), fill="blue", alpha="0.5") +
    geom_line(data=p1, aes(x=bin, y=min), lwd=lw) +
    geom_line(data=p1, aes(x=bin, y=med), lwd=lw) +
    geom_line(data=p1, aes(x=bin, y=max), lwd=lw) +
    geom_point(data=p1, aes(x=bin, y=med), size = 1, stroke = 0, shape = 16) +
    facet_grid(sp ~ ., labeller = labeller(sp = spp_short[1:12]))

plot2 <- empty +
    geom_ribbon(data=p2, aes(x=bin, ymin=min, ymax=max), fill="blue", alpha="0.5") +
    geom_line(data=p2, aes(x=bin, y=min), lwd=lw) +
    geom_line(data=p2, aes(x=bin, y=med), lwd=lw) +
    geom_line(data=p2, aes(x=bin, y=max), lwd=lw) +
    geom_point(data=p2, aes(x=bin, y=med), size = 1, stroke = 0, shape = 16) +
    facet_grid(sp ~ ., labeller = labeller(sp = spp_short[13:24]))
  
  # add vertical lines at ice ages

plot_nm <- paste('Figs/tseries_niche_pc1_pg1_', day, '.pdf', sep='')
pdf(plot_nm, width=8, height=11)
print(plot1)
dev.off()

plot_nm <- paste('Figs/tseries_niche_pc1_pg2_', day, '.pdf', sep='')
pdf(plot_nm, width=8, height=11)
print(plot2)
dev.off()
