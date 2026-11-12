source("R/compare_gold.R")

png('ethiopia_chpt_vs_onset.png', height=600,width=600)
par(mfrow=c(2,2))

par(mgp=c(2,.4,0))
d <- gold_chpt_pairs(run)
plot_gold_vs_chpt(d, main = "all", xlim = c(100,260))
d <- gold_chpt_pairs(run, lon_range=c(0,40), lat_range=c(8,100))
plot_gold_vs_chpt(d, main = "northwest (>8N,<40E)", xlim = c(100,260))
d <- gold_chpt_pairs(run, dev_range=c(.3,Inf))
plot_gold_vs_chpt(d, main = "dev > 0.3", xlim = c(100,260))
d <- gold_chpt_pairs(run, lon_range=c(0,40), lat_range=c(8,100), dev_range=c(.3,Inf))
plot_gold_vs_chpt(d, main = "northwest (>8N,<40E); dev > 0.3", xlim = c(100,260))

dev.off()
