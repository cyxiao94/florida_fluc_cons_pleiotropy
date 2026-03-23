
library(ggplot2)
library(data.table)
library(ggsignif)
library(ggpubr)

#PC1 only
euclidean_distance=NULL
for (model in 1:4){
    temp_f=fread(paste0("../model",model,"/04analysis/weighted_Euclidean_distant_PC1.txt"))
    temp_f$model=paste0("model",model)
    euclidean_distance=rbind(euclidean_distance,temp_f)
}

euclidean_distance1=melt(euclidean_distance,id.vars=c("set","model"), measure.var=c("anc_fluc","anc_cons","fluc_cons"))
colnames(euclidean_distance1)=c("set","model","comparison","distance")
euclidean_distance1$model=factor(euclidean_distance1$model)
euclidean_distance1$comparison=factor(euclidean_distance1$comparison)


# Step 1: Identify the largest distance for each replicate
euclidean_distance[, Largest := apply(.SD, 1, function(row) names(row)[which.max(row)]), 
   by = model, .SDcols = c("anc_fluc", "anc_cons", "fluc_cons")]

# Step 2: Calculate proportions within each model
proportion_data <- euclidean_distance[, .N, by = .(model, Largest)]
proportion_data[, Proportion := N / 100]  # 100 replicates per model

# Step 3: Plot the proportions
p1=ggplot(proportion_data, aes(x = model, y = Proportion, fill = Largest)) +
  geom_bar(stat = "identity")+
  geom_text(aes(label=Proportion*100),size=3, hjust=0.5,position=position_stack(vjust=0.5))+
  labs(x = "Model", y = "Proportion of Replicates with Largest Distance", fill = "Largest Distance") +
  scale_y_continuous(labels = scales::percent) +
  theme_minimal()
ggsave("04weighted_proportion_seperation_PC1.pdf",plot=p1, width = 7, height = 7)

#PC1 and PC2
euclidean_distance=NULL
for (model in 1:4){
    temp_f=fread(paste0("../model",model,"/04analysis/weighted_Euclidean_distant_PC1_PC2.txt"))
    temp_f$model=paste0("model",model)
    euclidean_distance=rbind(euclidean_distance,temp_f)
}

euclidean_distance1=melt(euclidean_distance,id.vars=c("set","model"), measure.var=c("anc_fluc","anc_cons","fluc_cons"))
colnames(euclidean_distance1)=c("set","model","comparison","distance")
euclidean_distance1$model=factor(euclidean_distance1$model)
euclidean_distance1$comparison=factor(euclidean_distance1$comparison)


# Step 1: Identify the largest distance for each replicate
euclidean_distance[, Largest := apply(.SD, 1, function(row) names(row)[which.max(row)]), 
   by = model, .SDcols = c("anc_fluc", "anc_cons", "fluc_cons")]

# Step 2: Calculate proportions within each model
proportion_data <- euclidean_distance[, .N, by = .(model, Largest)]
proportion_data[, Proportion := N / 100]  # 100 replicates per model

# Step 3: Plot the proportions
p12=ggplot(proportion_data, aes(x = model, y = Proportion, fill = Largest)) +
  geom_bar(stat = "identity")+
  geom_text(aes(label=Proportion*100),size=3, hjust=0.5,position=position_stack(vjust=0.5))+
  labs(x = "Model", y = "Proportion of Replicates with Largest Distance", fill = "Largest Distance") +
  scale_y_continuous(labels = scales::percent) +
  theme_minimal()
ggsave("04weighted_proportion_seperation_PC1_PC2.pdf",plot=p12, width = 7, height = 7)

#PC1 to PC20 
euclidean_distance=NULL
for (model in 1:4){
    temp_f=fread(paste0("../model",model,"/04analysis/weighted_Euclidean_distant_PC1_to_PC20.txt"))
    temp_f$model=paste0("model",model)
    euclidean_distance=rbind(euclidean_distance,temp_f)
}

euclidean_distance1=melt(euclidean_distance,id.vars=c("set","model"), measure.var=c("anc_fluc","anc_cons","fluc_cons"))
colnames(euclidean_distance1)=c("set","model","comparison","distance")
euclidean_distance1$model=factor(euclidean_distance1$model)
euclidean_distance1$comparison=factor(euclidean_distance1$comparison)


# Step 1: Identify the largest distance for each replicate
euclidean_distance[, Largest := apply(.SD, 1, function(row) names(row)[which.max(row)]), 
   by = model, .SDcols = c("anc_fluc", "anc_cons", "fluc_cons")]

# Step 2: Calculate proportions within each model
proportion_data <- euclidean_distance[, .N, by = .(model, Largest)]
proportion_data[, Proportion := N / 100]  # 100 replicates per model

# Step 3: Plot the proportions
p120=ggplot(proportion_data, aes(x = model, y = Proportion, fill = Largest)) +
  geom_bar(stat = "identity")+
  geom_text(aes(label=Proportion*100),size=3, hjust=0.5,position=position_stack(vjust=0.5))+
  labs(x = "Model", y = "Proportion of Replicates with Largest Distance", fill = "Largest Distance") +
  scale_y_continuous(labels = scales::percent) +
  theme_minimal()
ggsave("04weighted_proportion_seperation_PC1_PC20.pdf",plot=p120, width = 7, height = 7)



