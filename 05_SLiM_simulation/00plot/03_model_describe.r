
#library packages

library(scales)
library(ggplot2)
library(ggpubr)

#polygenicity
# Sample data frame with proportions
data <- data.frame(
  category = c("Lab", "23", "2818"),
  proportion = c(0.3333, 0.3333, 0.3333)  # Each proportion represents a part of the bar
)

# Plot a single stacked bar with three colors
p_poly_equal=ggplot(data, aes(x = 1, y = proportion, fill = category)) +
  geom_bar(stat = "identity", width = 0.5) +
  scale_fill_manual(values = c("Lab" = "#8c0000", "23" = "#4fc0c6", "2818" = "#e8a11c")) +
  labs(x = NULL, y = NULL) +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        legend.position = "none",
        axis.ticks = element_blank(), axis.text = element_blank()) +
  coord_flip()
  # Sample data frame with proportions
data <- data.frame(
  category = c("Lab", "23", "2818"),
  proportion = c(0.1, 0.45, 0.45)  # Each proportion represents a part of the bar
)

# Plot a single stacked bar with three colors
p_poly_unequal=ggplot(data, aes(x = 1, y = proportion, fill = category)) +
  geom_bar(stat = "identity", width = 0.5) +
  scale_fill_manual(values = c("Lab" = "#8c0000", "23" = "#4fc0c6", "2818" = "#e8a11c")) +
  labs(x = NULL, y = NULL) +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        legend.position = "none",
        axis.ticks = element_blank(), axis.text = element_blank()) +
  coord_flip()

#pleiotropy
data <- data.frame(
  category = c("Lab", "23", "2818"),
  proportion = c(1, 1, 1)  # Each proportion represents a part of the bar
)


# Plot a single stacked bar with three colors
p_poly_pleio=ggplot(data, aes(x = 1, y = proportion, fill = category)) +
  geom_bar(stat = "identity", width = 0.2, position = "dodge") +
  scale_fill_manual(values = c("Lab" = "#8c0000", "23" = "#4fc0c6", "2818" = "#e8a11c")) +
  labs(x = NULL, y = NULL) +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        legend.position = "none",
        axis.ticks = element_blank(), axis.text = element_blank()) +
  coord_flip()





#model1
#generate data for distribtion-model1
anc_mean=0.0
anc_sd=20
fit1_mean=anc_mean+3*anc_sd
fit1_sd=anc_sd
fit2_mean=anc_mean+anc_sd
fit2_sd=3.6*anc_sd

#anc_distrbution<-data.frame()

traits=seq(-120,120, length=1001)
anc=dnorm(traits, mean=anc_mean, sd=anc_sd)
fit1=dnorm(traits, mean=fit1_mean, sd =fit1_sd)
fit2=dnorm(traits, mean=fit2_mean, sd =fit2_sd)
anc_distrbution<-data.frame(phenotypic_value=traits,
                            density=anc)
fit_distrbution<-data.frame(phenotypic_value=traits,
                            fitness1=fit1,
                            fitness2=fit2
                            )


#anc_fit_distrbution$distribution<-as.factor(anc_fit_distrbution$distribution)
plot_fitness=function(traits,den_anc,fit,color_anc,color_fit){
    ymax1=dnorm(0,0,20)
    anc_distrbution=data.frame(phenotypic_value=traits,
                            density=anc)
    fit_distrbution=data.frame(phenotypic_value=traits,
                            fitness=fit)
                            #fitness2=fit2
                            
    #p1=ggplot()+
    p_fit=ggplot()+
        geom_line(data=anc_distrbution, aes(x=phenotypic_value, y=density), color=color_anc)+
        annotate("segment",x=0.0, xend = 60,y=0.02,yend=0.02,color="#219ebc",linewidth=0.5, 
        arrow=arrow(ends = "both", angle = 40, length = unit(0.2,"cm")))+
        geom_line(data=fit_distrbution, aes(x=phenotypic_value,y=rescale(fitness, c(0,ymax1))), color=color_fit)+
        annotate("text",x=-20,y=0.005, label="Ancestral Phenotypic Distribution",color=color_anc,size=2)+
        annotate("text",x=60,y=0.015, label="Fitness Function", color=color_fit,size=2)+
        xlab("Phenotypic value")+
        scale_y_continuous(breaks = seq(0,ymax1, 0.01), limits = c(0,ymax1), name="Density",
                     sec.axis = sec_axis(~./0.02+0.5, name = "Fitness"))+
        theme_light()+
        theme(panel.grid = element_blank(),
            axis.title=element_text(size=15),
            axis.text = element_text(size=12),
            axis.title.y.right = element_text(color="#AF1B3F",size=15),
            axis.text.y.right = element_text(color = "#AF1B3F",size=12))
return(p_fit)
}








#model1
p_lab=plot_fitness(traits,anc,fit1,color_anc="#6e326e",color_fit="#8c0000")
p_23=plot_fitness(traits,anc,fit1,color_anc="#6e326e",color_fit="#4fc0c6")
p_2818=plot_fitness(traits,anc,fit1,color_anc="#6e326e",color_fit="#e8a11c")

p_fitness <- ggarrange(p_lab+ggtitle("Lab"),
            p_2818+ggtitle("2818"),
            p_23+ggtitle("23"),nrow=3, heights=c(3,3,3))

p_fitness<- annotate_figure(p_fitness, 
                            top = text_grob("Fitness function", size = 15, face = "bold", hjust = 1))
                      
p=ggarrange(p_poly_equal,
            p_fitness ,nrow=2, heights=c(1,9))
ggsave("03fitness_function_model1.pdf",plot=p, width=4, height=7)


#model2
p_lab=plot_fitness(traits,anc,fit1,color_anc="#6e326e",color_fit="#8c0000")
p_23=plot_fitness(traits,anc,fit2,color_anc="#6e326e",color_fit="#4fc0c6")
p_2818=plot_fitness(traits,anc,fit2,color_anc="#6e326e",color_fit="#e8a11c")

p_fitness <- ggarrange(p_lab+ggtitle("Lab"),
            p_2818+ggtitle("2818"),
            p_23+ggtitle("23"),nrow=3, heights=c(3,3,3))

p_fitness<- annotate_figure(p_fitness, 
                            top = text_grob("Fitness function", size = 15, face = "bold", hjust = 1))
                      
p=ggarrange(p_poly_equal,
            p_fitness ,nrow=2, heights=c(1,9))
ggsave("03fitness_function_model2.pdf",plot=p, width=4, height=7)

#model3
p_lab=plot_fitness(traits,anc,fit1,color_anc="#6e326e",color_fit="#8c0000")
p_23=plot_fitness(traits,anc,fit2,color_anc="#6e326e",color_fit="#4fc0c6")
p_2818=plot_fitness(traits,anc,fit2,color_anc="#6e326e",color_fit="#e8a11c")

p_fitness <- ggarrange(p_lab+ggtitle("Lab"),
            p_2818+ggtitle("2818"),
            p_23+ggtitle("23"),nrow=3, heights=c(3,3,3))

p_fitness<- annotate_figure(p_fitness, 
                            top = text_grob("Fitness function", size = 15, face = "bold", hjust = 1))
                      
p=ggarrange(p_poly_unequal,
            p_fitness ,nrow=2, heights=c(1,9))
ggsave("03fitness_function_model3.pdf",plot=p, width=4, height=7)


#model4
p_lab=plot_fitness(traits,anc,fit1,color_anc="#6e326e",color_fit="#8c0000")
p_23=plot_fitness(traits,anc,fit2,color_anc="#6e326e",color_fit="#4fc0c6")
p_2818=plot_fitness(traits,anc,fit2,color_anc="#6e326e",color_fit="#e8a11c")

p_fitness <- ggarrange(p_lab+ggtitle("Lab"),
            p_2818+ggtitle("2818"),
            p_23+ggtitle("23"),nrow=3, heights=c(3,3,3))

p_fitness<- annotate_figure(p_fitness, 
                            top = text_grob("Fitness function", size = 15, face = "bold", hjust = 1))
                      
p=ggarrange(p_poly_pleio,
            p_fitness ,nrow=2, heights=c(1,9))
ggsave("03fitness_function_model4.pdf",plot=p, width=4, height=7)



