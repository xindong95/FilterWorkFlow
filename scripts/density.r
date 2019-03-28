library("ggplot2")
args=commandArgs(T)
data <- read.table(args[1])
#data <- read.table("test.txt")
data <- abs(data[1])
max(data)
data <- data[data<1000]
max(data)
data<- data.frame(data)
colnames(data) <- "FragmentSize"
ggplot(data, aes(x = FragmentSize))+
geom_density(alpha=0.3)+
labs(title = "Fragment Size Distribution",
       x = "Fragment Size",
       y = "Density")+
scale_x_continuous(limits=c(0,1000),breaks = seq(0,1000,200))
ggsave(args[2])
