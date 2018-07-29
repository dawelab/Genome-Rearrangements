test<-read.table(lambda_aligned_filtered_sum.txt",sep = "\t")
head(test)
size<-test[['V5']]
nrow(test)
head(size)
library(ggplot2)
img = ggplot(data.frame(size), aes(x=size)) + 
  geom_histogram(color="grey37", fill="lightblue",bins = 50, alpha=0.8,show.legend = FALSE) + 
  geom_rug(color="tomato2") + 
  geom_vline(aes(xintercept=mean(size), color="mean"), linetype="dashed",show.legend = FALSE) + 
  geom_vline(aes(xintercept=median(size), color="median"), linetype="dashed") + 
  scale_x_continuous(breaks = pretty(size, n = 12)) + 
  labs(x="Size (bp)", y = "Number of fragments", color = "") + geom_density(alpha=0.6) + 
  scale_fill_manual(name = "statistics", values = c(mean = "red",median = "royalblue")) + 
  theme_classic() +
  ylim(0,20000) + 
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold"))

ggsave(filename="distribution.svg", plot=img, device = "svg")
