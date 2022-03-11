#Written by Jeffrey Malins

library(dplyr); library(tidyr); library(ggplot2)

#Multiple Regression Analyses - Mean and Trial-by-Trial Activation Variability

#Print
nuisancevars <- c("scale(Age)",
                  "scale(average_motion_TR)",
                  "scale(NumTrials_V)")

for (ROI in c("LIFGop","LIFGpt","LSPL","LMTG","LITG","LSTG","LThal")){
  NeuroimagingSample$ROI.meanV=NeuroimagingSample[[paste(ROI,"_meanV",sep="")]]
  NeuroimagingSample$ROI.sdV=NeuroimagingSample[[paste(ROI,"_sdV",sep="")]]
  
  print(ROI)
  
  m.full <- lm(scale(WJ3LWID_RS) ~ 
                 scale(Age) + 
                 scale(average_motion_TR) +
                 scale(NumTrials_V) +
                 scale(ROI.meanV) +
                 scale(ROI.sdV), 
               data = NeuroimagingSample); summary(m.full)
  
  m.dropfirst <- MASS::dropterm(m.full,test="F",scope=nuisancevars)
  P1 <- max(m.dropfirst$`Pr(F)`, na.rm = TRUE)
  
  if (P1 >= .05) {
    varDrop1 <- noquote(row.names(m.dropfirst)[which(m.dropfirst$`Pr(F)`==P1)]); print(varDrop1)
    m.drop1 = update(m.full, paste(". ~ . -",varDrop1))
    
    m.dropsecond <- MASS::dropterm(m.drop1,test="F",scope=nuisancevars[-which(nuisancevars==varDrop1)])
    P2 <- max(m.dropsecond$`Pr(F)`, na.rm = TRUE)
    
    if (P2 >= .05){
      varDrop2 <- noquote(row.names(m.dropsecond)[which(m.dropsecond$`Pr(F)`==P2)]); print(varDrop2)
      m.drop2 <- update(m.drop1, paste(". ~ . -",varDrop2))
      
      m.dropthird <- MASS::dropterm(m.drop2,test="F",scope=nuisancevars[-which(nuisancevars==varDrop1 | nuisancevars==varDrop2)])
      P3 <- max(m.dropthird$`Pr(F)`, na.rm = TRUE)
      
      if (P3 >= .05){
        varDrop3 <- noquote(row.names(m.dropthird)[which(m.dropthird$`Pr(F)`==P3)]); print(varDrop3)
        m.final <- update(m.drop2, paste(". ~ . -",varDrop3))
        
      } else {
        m.final <- m.drop2
      } 
    } else {
      m.final <- m.drop1
    }
  } else {
    m.final <- m.full
  }
  
  mfinal.f <- MASS::dropterm(m.final,test="F"); mfinal.f
  pvals.V <- c(pvals.V,mfinal.f$`Pr(F)`[-1])
  
  m.noSD <- update(m.final, . ~ . -scale(ROI.sdV))
  adj.R2.change.SD = summary(m.final)$adj.r.squared - summary(m.noSD)$adj.r.squared
  comp.AIC.SD <- AIC(m.final,m.noSD)
  AIC.change.SD = comp.AIC.SD$AIC[1]-comp.AIC.SD$AIC[2]
  comp.BIC.SD <- BIC(m.final,m.noSD)
  BIC.change.SD = comp.BIC.SD$BIC[1]-comp.BIC.SD$BIC[2]
  
  m.noMean <- update(m.final, . ~ . -scale(ROI.meanV))
  adj.R2.change.mean = summary(m.final)$adj.r.squared - summary(m.noMean)$adj.r.squared
  comp.AIC.mean <- AIC(m.final,m.noMean)
  AIC.change.mean = comp.AIC.mean$AIC[1]-comp.AIC.mean$AIC[2]
  comp.BIC.mean <- BIC(m.final,m.noMean)
  BIC.change.mean = comp.BIC.mean$BIC[1]-comp.BIC.mean$BIC[2]
  
  #Values for table
  nrow.SD <- which(row.names(mfinal.f)=="scale(ROI.sdV)")
  table.Beta.SD=round(summary(m.final)$coefficients[nrow.SD,1],3) #Beta
  table.SE.SD=round(summary(m.final)$coefficients[nrow.SD,2],3) #SE
  table.AIC.change.SD=round(AIC.change.SD,3) #AIC
  table.BIC.change.SD=round(BIC.change.SD,3) #BIC
  table.F.SD=round(mfinal.f$`F Value`[nrow.SD],3) #F
  table.p.SD=round(mfinal.f$`Pr(F)`[nrow.SD],3) #p-value
  table.R2adj.SD=round(adj.R2.change.SD,3) #Adjusted R-squared
  
  nrow.mean <- which(row.names(mfinal.f)=="scale(ROI.meanV)")
  table.Beta.mean=round(summary(m.final)$coefficients[nrow.mean,1],3) #Beta
  table.SE.mean=round(summary(m.final)$coefficients[nrow.mean,2],3) #SE
  table.AIC.change.mean=round(AIC.change.mean,3) #AIC
  table.BIC.change.mean=round(BIC.change.mean,3) #BIC
  table.F.mean=round(mfinal.f$`F Value`[nrow.mean],3) #F
  table.p.mean=round(mfinal.f$`Pr(F)`[nrow.mean],3) #p-value
  table.R2adj.mean=round(adj.R2.change.mean,3) #Adjusted R-squared
  
  print(c(table.Beta.mean,table.SE.mean,table.AIC.change.mean,table.BIC.change.mean,table.R2adj.mean,table.F.mean,table.p.mean))
  print(c(table.Beta.SD,table.SE.SD,table.AIC.change.SD,table.BIC.change.SD,table.R2adj.SD,table.F.SD,table.p.SD))
}



#Speech
nuisancevars <- c("scale(Age)",
                  "scale(average_motion_TR)",
                  "scale(NumTrials_A)")

for (ROI in c("LIFGop","LIFGpt","LSPL","LMTG","LITG","LSTG","LThal")){
  NeuroimagingSample$ROI.meanA=NeuroimagingSample[[paste(ROI,"_meanA",sep="")]]
  NeuroimagingSample$ROI.sdA=NeuroimagingSample[[paste(ROI,"_sdA",sep="")]]
  
  print(ROI)
  
  m.full <- lm(scale(WJ3LWID_RS) ~ 
                 scale(Age) + 
                 scale(average_motion_TR) +
                 scale(NumTrials_A) +
                 scale(ROI.meanA) +
                 scale(ROI.sdA), 
               data = NeuroimagingSample); summary(m.full)
  
  m.dropfirst <- MASS::dropterm(m.full,test="F",scope=nuisancevars)
  P1 <- max(m.dropfirst$`Pr(F)`, na.rm = TRUE)
  
  if (P1 >= .05) {
    varDrop1 <- noquote(row.names(m.dropfirst)[which(m.dropfirst$`Pr(F)`==P1)]); print(varDrop1)
    m.drop1 = update(m.full, paste(". ~ . -",varDrop1))
    
    m.dropsecond <- MASS::dropterm(m.drop1,test="F",scope=nuisancevars[-which(nuisancevars==varDrop1)])
    P2 <- max(m.dropsecond$`Pr(F)`, na.rm = TRUE)
    
    if (P2 >= .05){
      varDrop2 <- noquote(row.names(m.dropsecond)[which(m.dropsecond$`Pr(F)`==P2)]); print(varDrop2)
      m.drop2 <- update(m.drop1, paste(". ~ . -",varDrop2))
      
      m.dropthird <- MASS::dropterm(m.drop2,test="F",scope=nuisancevars[-which(nuisancevars==varDrop1 | nuisancevars==varDrop2)])
      P3 <- max(m.dropthird$`Pr(F)`, na.rm = TRUE)
      
      if (P3 >= .05){
        varDrop3 <- noquote(row.names(m.dropthird)[which(m.dropthird$`Pr(F)`==P3)]); print(varDrop3)
        m.final <- update(m.drop2, paste(". ~ . -",varDrop3))
        
      } else {
        m.final <- m.drop2
      } 
    } else {
      m.final <- m.drop1
    }
  } else {
    m.final <- m.full
  }
  
  mfinal.f <- MASS::dropterm(m.final,test="F"); mfinal.f
  pvals.A <- c(pvals.A,mfinal.f$`Pr(F)`[-1])
  
  m.noSD <- update(m.final, . ~ . -scale(ROI.sdA))
  adj.R2.change.SD = summary(m.final)$adj.r.squared - summary(m.noSD)$adj.r.squared
  comp.AIC.SD <- AIC(m.final,m.noSD)
  AIC.change.SD = comp.AIC.SD$AIC[1]-comp.AIC.SD$AIC[2]
  comp.BIC.SD <- BIC(m.final,m.noSD)
  BIC.change.SD = comp.BIC.SD$BIC[1]-comp.BIC.SD$BIC[2]
  
  m.noMean <- update(m.final, . ~ . -scale(ROI.meanA))
  adj.R2.change.mean = summary(m.final)$adj.r.squared - summary(m.noMean)$adj.r.squared
  comp.AIC.mean <- AIC(m.final,m.noMean)
  AIC.change.mean = comp.AIC.mean$AIC[1]-comp.AIC.mean$AIC[2]
  comp.BIC.mean <- BIC(m.final,m.noMean)
  BIC.change.mean = comp.BIC.mean$BIC[1]-comp.BIC.mean$BIC[2]
  
  #Values for table
  nrow.SD <- which(row.names(mfinal.f)=="scale(ROI.sdA)")
  table.Beta.SD=round(summary(m.final)$coefficients[nrow.SD,1],3) #Beta
  table.SE.SD=round(summary(m.final)$coefficients[nrow.SD,2],3) #SE
  table.AIC.change.SD=round(AIC.change.SD,3) #AIC
  table.BIC.change.SD=round(BIC.change.SD,3) #BIC
  table.F.SD=round(mfinal.f$`F Value`[nrow.SD],3) #F
  table.p.SD=round(mfinal.f$`Pr(F)`[nrow.SD],3) #p-value
  table.R2adj.SD=round(adj.R2.change.SD,3) #Adjusted R-squared
  
  nrow.mean <- which(row.names(mfinal.f)=="scale(ROI.meanA)")
  table.Beta.mean=round(summary(m.final)$coefficients[nrow.mean,1],3) #Beta
  table.SE.mean=round(summary(m.final)$coefficients[nrow.mean,2],3) #SE
  table.AIC.change.mean=round(AIC.change.mean,3) #AIC
  table.BIC.change.mean=round(BIC.change.mean,3) #BIC
  table.F.mean=round(mfinal.f$`F Value`[nrow.mean],3) #F
  table.p.mean=round(mfinal.f$`Pr(F)`[nrow.mean],3) #p-value
  table.R2adj.mean=round(adj.R2.change.mean,3) #Adjusted R-squared
  
  print(c(table.Beta.mean,table.SE.mean,table.AIC.change.mean,table.BIC.change.mean,table.R2adj.mean,table.F.mean,table.p.mean))
  print(c(table.Beta.SD,table.SE.SD,table.AIC.change.SD,table.BIC.change.SD,table.R2adj.SD,table.F.SD,table.p.SD))
}



#Assessing whether removal of the activation variability term from the model changes signifiance of mean activation term

#Print
nuisancevars <- c("scale(Age)",
                  "scale(average_motion_TR)",
                  "scale(NumTrials_V)")

for (ROI in c("LIFGop","LIFGpt","LSPL","LMTG","LITG","LSTG","LThal")){
  NeuroimagingSample$ROI.meanV=NeuroimagingSample[[paste(ROI,"_meanV",sep="")]]
  
  print(ROI)
  
  m.full <- lm(scale(WJ3LWID_RS) ~ 
                 scale(Age) + 
                 scale(average_motion_TR) +
                 scale(NumTrials_V) +
                 scale(ROI.meanV), 
               data = NeuroimagingSample); summary(m.full)
  
  m.dropfirst <- MASS::dropterm(m.full,test="F",scope=nuisancevars)
  P1 <- max(m.dropfirst$`Pr(F)`, na.rm = TRUE)
  
  if (P1 >= .05) {
    varDrop1 <- noquote(row.names(m.dropfirst)[which(m.dropfirst$`Pr(F)`==P1)]); print(varDrop1)
    m.drop1 = update(m.full, paste(". ~ . -",varDrop1))
    
    m.dropsecond <- MASS::dropterm(m.drop1,test="F",scope=nuisancevars[-which(nuisancevars==varDrop1)])
    P2 <- max(m.dropsecond$`Pr(F)`, na.rm = TRUE)
    
    if (P2 >= .05){
      varDrop2 <- noquote(row.names(m.dropsecond)[which(m.dropsecond$`Pr(F)`==P2)]); print(varDrop2)
      m.drop2 <- update(m.drop1, paste(". ~ . -",varDrop2))
      
      m.dropthird <- MASS::dropterm(m.drop2,test="F",scope=nuisancevars[-which(nuisancevars==varDrop1 | nuisancevars==varDrop2)])
      P3 <- max(m.dropthird$`Pr(F)`, na.rm = TRUE)
      
      if (P3 >= .05){
        varDrop3 <- noquote(row.names(m.dropthird)[which(m.dropthird$`Pr(F)`==P3)]); print(varDrop3)
        m.final <- update(m.drop2, paste(". ~ . -",varDrop3))
        
      } else {
        m.final <- m.drop2
      } 
    } else {
      m.final <- m.drop1
    }
  } else {
    m.final <- m.full
  }
  
  mfinal.f <- MASS::dropterm(m.final,test="F"); mfinal.f
  pvals.V <- c(pvals.V,mfinal.f$`Pr(F)`[-1])
  
  m.noMean <- update(m.final, . ~ . -scale(ROI.meanV))
  
  adj.R2.change = summary(m.final)$adj.r.squared - summary(m.noMean)$adj.r.squared
  comp.AIC <- AIC(m.final,m.noMean)
  AIC.change = comp.AIC$AIC[1]-comp.AIC$AIC[2]
  comp.BIC <- BIC(m.final,m.noMean)
  BIC.change = comp.BIC$BIC[1]-comp.BIC$BIC[2]
  
  #Values for table
  nrow <- which(row.names(mfinal.f)=="scale(ROI.meanV)")
  table.Beta=round(summary(m.final)$coefficients[nrow,1],3) #Beta
  table.SE=round(summary(m.final)$coefficients[nrow,2],3) #SE
  table.AIC.change=round(AIC.change,3) #AIC
  table.BIC.change=round(BIC.change,3) #BIC
  table.F=round(mfinal.f$`F Value`[nrow],3) #F
  table.p=round(mfinal.f$`Pr(F)`[nrow],3) #p-value
  table.R2adj=round(adj.R2.change,3) #Adjusted R-squared
  print(c(table.AIC.change,table.BIC.change,table.R2adj,table.F,table.p))
}



#Speech
nuisancevars <- c("scale(Age)",
                  "scale(average_motion_TR)",
                  "scale(NumTrials_A)")

for (ROI in c("LIFGop","LIFGpt","LSPL","LMTG","LITG","LSTG","LThal")){
  NeuroimagingSample$ROI.meanA=NeuroimagingSample[[paste(ROI,"_meanA",sep="")]]
  
  print(ROI)
  
  m.full <- lm(scale(WJ3LWID_RS) ~ 
                 scale(Age) + 
                 scale(average_motion_TR) +
                 scale(NumTrials_A) +
                 scale(ROI.meanA), 
               data = NeuroimagingSample); summary(m.full)
  
  m.dropfirst <- MASS::dropterm(m.full,test="F",scope=nuisancevars)
  P1 <- max(m.dropfirst$`Pr(F)`, na.rm = TRUE)
  
  if (P1 >= .05) {
    varDrop1 <- noquote(row.names(m.dropfirst)[which(m.dropfirst$`Pr(F)`==P1)]); print(varDrop1)
    m.drop1 = update(m.full, paste(". ~ . -",varDrop1))
    
    m.dropsecond <- MASS::dropterm(m.drop1,test="F",scope=nuisancevars[-which(nuisancevars==varDrop1)])
    P2 <- max(m.dropsecond$`Pr(F)`, na.rm = TRUE)
    
    if (P2 >= .05){
      varDrop2 <- noquote(row.names(m.dropsecond)[which(m.dropsecond$`Pr(F)`==P2)]); print(varDrop2)
      m.drop2 <- update(m.drop1, paste(". ~ . -",varDrop2))
      
      m.dropthird <- MASS::dropterm(m.drop2,test="F",scope=nuisancevars[-which(nuisancevars==varDrop1 | nuisancevars==varDrop2)])
      P3 <- max(m.dropthird$`Pr(F)`, na.rm = TRUE)
      
      if (P3 >= .05){
        varDrop3 <- noquote(row.names(m.dropthird)[which(m.dropthird$`Pr(F)`==P3)]); print(varDrop3)
        m.final <- update(m.drop2, paste(". ~ . -",varDrop3))
        
      } else {
        m.final <- m.drop2
      } 
    } else {
      m.final <- m.drop1
    }
  } else {
    m.final <- m.full
  }
  
  mfinal.f <- MASS::dropterm(m.final,test="F"); mfinal.f
  pvals.A <- c(pvals.A,mfinal.f$`Pr(F)`[-1])
  
  m.noMean <- update(m.final, . ~ . -scale(ROI.meanA))
  
  adj.R2.change = summary(m.final)$adj.r.squared - summary(m.noMean)$adj.r.squared
  comp.AIC <- AIC(m.final,m.noMean)
  AIC.change = comp.AIC$AIC[1]-comp.AIC$AIC[2]
  comp.BIC <- BIC(m.final,m.noMean)
  BIC.change = comp.BIC$BIC[1]-comp.BIC$BIC[2]
  
  #Values for table
  nrow <- which(row.names(mfinal.f)=="scale(ROI.meanA)")
  table.Beta=round(summary(m.final)$coefficients[nrow,1],3) #Beta
  table.SE=round(summary(m.final)$coefficients[nrow,2],3) #SE
  table.AIC.change=round(AIC.change,3) #AIC
  table.BIC.change=round(BIC.change,3) #BIC
  table.F=round(mfinal.f$`F Value`[nrow],3) #F
  table.p=round(mfinal.f$`Pr(F)`[nrow],3) #p-value
  table.R2adj=round(adj.R2.change,3) #Adjusted R-squared
  print(c(table.AIC.change,table.BIC.change,table.R2adj,table.F,table.p))
}





#Plotting the partial correlation results
m.lw <- lm(scale(WJ3LWID_RS) ~ 
             scale(NumTrials_V) +
             scale(LIFGpt_meanV), 
           data = NeuroimagingSample)
m.sd <- lm(scale(LIFGpt_sdV) ~ 
             scale(NumTrials_V) +
             scale(LIFGpt_meanV), 
           data = NeuroimagingSample)

NeuroimagingSample$res.lw <- residuals.lm(m.lw)
NeuroimagingSample$res.sd <- residuals.lm(m.sd)

ggplot(NeuroimagingSample,aes(res.sd,res.lw)) + 
  geom_point() +
  geom_smooth(method="lm", color="red") +
  labs(x="",y="") +
  coord_cartesian(ylim=c(-2.1,2.1)) +
  theme_bw(15) +
  #eliminates background, gridlines, and chart border
  theme(
    plot.background = element_blank()
    ,panel.grid.major = element_blank()
    ,panel.grid.minor = element_blank()
    ,panel.border = element_blank()
  ) +
  #draws x and y axis line
  theme(axis.line.x = element_line(color="black", size = .3),
        axis.line.y = element_line(color="black", size = .3))
ggsave("LIFGpt_Plot.pdf",width=7,height=6, useDingbats=FALSE)


