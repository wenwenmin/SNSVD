source('F1_Sim_Test_50Times.R')
# packages
library(ggplot2)
library(gridExtra)
library(cowplot)

# --------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------
getFigDat = function(Res50,lambda){
   specificity.u = Sensitivity.u = accuracy.u = matrix(0,nrow = 4, ncol = length(lambda))
   specificity.v = Sensitivity.v = accuracy.v = matrix(0,nrow = 4, ncol = length(lambda))
   
   for(kk in 1:50){
      Res = Res50[[kk]]
      
      Sensitivity.u =  Res$Sensitivity.u + Sensitivity.u
      specificity.u = Res$specificity.u + specificity.u
      accuracy.u = Res$accuracy.u + accuracy.u
      
      Sensitivity.v = Res$Sensitivity.v + Sensitivity.v
      specificity.v = Res$specificity.v + specificity.v
      accuracy.v = Res$accuracy.v + accuracy.v
   }
   
   Sensitivity.u = Sensitivity.u/50 
   specificity.u = specificity.u/50
   accuracy.u = accuracy.u/50
   
   Sensitivity.v = Sensitivity.v/50 
   specificity.v = specificity.v/50
   accuracy.v = accuracy.v/50
   
   a1 = data.frame(lambda = lambda, method = rep("SNSVD", length(lambda)), spe.u = specificity.u[1,], sen.u = Sensitivity.u[1,], acc.u = accuracy.u[1,])
   a2 = data.frame(lambda = lambda, method = rep("L0SVD", length(lambda)), spe.u = specificity.u[2,], sen.u = Sensitivity.u[2,], acc.u = accuracy.u[2,])
   a3 = data.frame(lambda = lambda, method = rep("ALSVD", length(lambda)), spe.u = specificity.u[3,], sen.u = Sensitivity.u[3,], acc.u = accuracy.u[3,])
   a4 = data.frame(lambda = lambda, method = rep("SCADSVD", length(lambda)), spe.u = specificity.u[4,], sen.u = Sensitivity.u[4,], acc.u = accuracy.u[4,])
   
   b1 = data.frame(lambda = lambda, method = rep("SNSVD", length(lambda)), spe.v = specificity.v[1,], sen.v = Sensitivity.v[1,], acc.v = accuracy.v[1,])
   b2 = data.frame(lambda = lambda, method = rep("L0SVD", length(lambda)), spe.v = specificity.v[2,], sen.v = Sensitivity.v[2,], acc.v = accuracy.v[2,])
   b3 = data.frame(lambda = lambda, method = rep("ALSVD", length(lambda)), spe.v = specificity.v[3,], sen.v = Sensitivity.v[3,], acc.v = accuracy.v[3,])
   b4 = data.frame(lambda = lambda, method = rep("SCADSVD", length(lambda)), spe.v = specificity.v[4,], sen.v = Sensitivity.v[4,], acc.v = accuracy.v[4,])
   
   data.u = rbind(a1,rbind(a2,rbind(a3,a4)))
   data.v = rbind(b1,rbind(b2,rbind(b3,b4)))
   
   return(list(data.u=data.u,data.v=data.v))
}
# --------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------

FigDat = getFigDat(Res50,lambda) 
data.u = FigDat$data.u
data.v = FigDat$data.v

# network-regularized norm for v

pp1 = ggplot(data.u, aes(x=lambda, y = sen.u, color = method, shape = method)) + geom_point(size=3) +
      xlab(expression(gamma))+ylab("Sensitivity") + ggtitle("u")
pp1 = pp1 + scale_shape(solid = FALSE) + geom_line(size=0.5,linetype = 1)

pp2 = ggplot(data.u, aes(x=lambda, y = spe.u, color = method, shape = method)) + geom_point(size=3) + 
      xlab(expression(gamma))+ylab("Specificity") + ggtitle("u") #+ theme(legend.position = "none")
pp2 = pp2 + scale_shape(solid = FALSE) + geom_line(size=0.5,linetype = 1) 

pp3 = ggplot(data.v, aes(x=lambda, y = sen.v, color = method, shape = method)) + geom_point(size=3) +
      xlab(expression(gamma)) + ylab("Sensitivity")  + ggtitle("v")
pp3 = pp3 + scale_shape(solid = FALSE) + geom_line(size=0.5,linetype = 1) 

pp4 = ggplot(data.v, aes(x=lambda, y = spe.v, color = method, shape = method)) + geom_point(size=3) + 
      xlab(expression(gamma)) + ylab("Specificity") + ggtitle("v") #+ theme(legend.position = "none") 

pp4 = pp4 + scale_shape(solid = FALSE) + geom_line(size=0.5,linetype = 1) 

# ACC
pp5 = ggplot(data.u, aes(x=lambda, y = acc.u, color = method, shape = method)) + geom_point(size=3) +
   xlab(expression(gamma))+ylab("Accuracy") + ggtitle("u")
pp5 = pp5 + scale_shape(solid = FALSE) + geom_line(size=0.5,linetype = 1) 

pp6 = ggplot(data.v, aes(x=lambda, y = acc.v, color = method, shape = method)) + geom_point(size=3) +
   xlab(expression(gamma)) + ylab("Accuracy")  + ggtitle("v")
pp6 = pp6 + scale_shape(solid = FALSE) + geom_line(size=0.5,linetype = 1) 


##------------------------------------------------------------------
pp1 = pp1 + theme(axis.title = element_text(face="plain", colour="black")) + theme_classic() + theme(legend.position = c(0.23, 0.26)) +
   theme(legend.text = element_text(colour = "black", angle = 0), legend.title = element_text(size = 0, colour = "white"))

pp2 = pp2 + theme(axis.title = element_text(face="plain", colour="black")) + theme_classic() + theme(legend.position = c(0.23, 0.26)) +
   theme(legend.text = element_text(colour = "black", angle = 0), legend.title = element_text(size = 0, colour = "white"))

pp3 = pp3 + theme(axis.title = element_text(face="plain", colour="black")) + theme_classic() + theme(legend.position = c(0.23, 0.26)) +
   theme(legend.text = element_text(colour = "black", angle = 0), legend.title = element_text(size = 0, colour = "white"))

pp4 = pp4 + theme(axis.title = element_text(face="plain", colour="black")) + theme_classic() + theme(legend.position = c(0.23, 0.26)) +
   theme(legend.text = element_text(colour = "black", angle = 0), legend.title = element_text(size = 0, colour = "white"))

pp5 = pp5 + theme(axis.title = element_text(face="plain", colour="black")) + theme_classic() + theme(legend.position = c(0.23, 0.26)) +
   theme(legend.text = element_text(colour = "black", angle = 0), legend.title = element_text(size = 0, colour = "white"))

pp6 = pp6 + theme(axis.title = element_text(face="plain", colour="black")) + theme_classic() + theme(legend.position = c(0.23, 0.26)) +
   theme(legend.text = element_text(colour = "black", angle = 0), legend.title = element_text(size = 0, colour = "white"))
##------------------------------------------------------------------

colors = c("red","blue", "seagreen","purple") 
II = c(1,2,6,0)

pp1 = pp1 + scale_shape_manual(values=II) + scale_color_manual(values = colors)
pp2 = pp2 + scale_shape_manual(values=II) + scale_color_manual(values = colors)
pp3 = pp3 + scale_shape_manual(values=II) + scale_color_manual(values = colors)
pp4 = pp4 + scale_shape_manual(values=II) + scale_color_manual(values = colors)
pp5 = pp5 + scale_shape_manual(values=II) + scale_color_manual(values = colors)
pp6 = pp6 + scale_shape_manual(values=II) + scale_color_manual(values = colors)

plot_grid(pp1,pp2, pp5, pp3,pp4, pp6, 
          ncol=3, nrow =2, align="hv", labels =c("","","",""),
          label_size = 20)
ggsave(file="SNSVD_Figure2+++.pdf",width=9, height=6)

# plot_grid(pp1,pp2,pp3,pp4,
#           ncol=2, nrow =2, align="hv", labels =c("","","",""),
#           label_size = 20)
# ggsave(file="SNSVD_Figure2--.pdf",width=7, height=6)