#### Initial ####
rm(list=ls())
require(plyr)
require(ggplot2)
require(cowplot)

#### Load data ####

load('commDFA.RData')
load('commDFB.RData')
load('commDFC.RData')
load('commDFD.RData')

commDF <- rbind(commDFA,commDFB,commDFC,commDFD)
rm(commDFA,commDFB,commDFC,commDFD)

#### Figure 2 ####

commSummDF<-ddply(commDF,c("landscape","treatment","effort"),summarize,N=length(community),
            gS=mean(gammaSimps2,na.rm=T),aS=mean(alphaSimps2,na.rm=T),
            gR=mean(gammaRich2,na.rm=T),aR=mean(alphaRich2,na.rm=T),
            tN=mean(totN2,na.rm=T),rM=mean(rangeMean2,na.rm=T),
            dN=mean(dispMean2,na.rm=T),tolM=mean(tolerMean2,na.rm=T),
            ex=mean(extinctions2,na.rm=T),g0=mean(gammaRich1,na.rm=T),
            d0=mean(dispMean1,na.rm=T),a0=mean(alphaSimps1,na.rm=T),
            exP=mean(extinctions2/gammaRich1,na.rm=T))

e1<-ggplot(commSummDF[commSummDF$treatment %in% c(1:3,5) & commSummDF$landscape==3,],aes(x=effort,y=ex,color=factor(treatment),shape=factor(treatment)))+
  geom_line(size=1)+
  geom_point(size=2)+
  scale_colour_manual(name="Management\nstrategy",values=c("red","blue","magenta","gray"),labels=c("Corridor establishment","Stepping-stone reserves","Reinforce reserves","Restore all locations"))+
  scale_shape_manual(name="Management\nstrategy",values=c("triangle","triangle","triangle","triangle"),labels=c("Corridor establishment","Stepping-stone reserves","Reinforce reserves","Restore everywhere"))+
  theme_classic()+
  scale_y_continuous(name="Average number of extinctions")+
  scale_x_continuous(name=expression(paste("Total area restored (\u00D7",2^12,")")))+
  coord_cartesian(ylim=c(0,2.5))+
  theme(legend.position = "none")

e2<-ggplot(commSummDF[commSummDF$treatment %in% c(1:3,5) & commSummDF$landscape==5,],aes(x=effort,y=ex,color=factor(treatment),shape=factor(treatment)))+
  geom_line(size=1)+
  geom_point(size=2)+
  scale_colour_manual(name="Management\nstrategy",values=c("red","blue","magenta","gray"),labels=c("Corridor establishment","Stepping-stone reserves","Reinforce reserves","Restore all locations"))+
  scale_shape_manual(name="Management\nstrategy",values=c("triangle","triangle","triangle","triangle"),labels=c("Corridor establishment","Stepping-stone reserves","Reinforce reserves","Restore everywhere"))+
  theme_classic()+
  scale_y_continuous(name="Average number of extinctions")+
  scale_x_continuous(name=expression(paste("Total area restored (\u00D7",2^12,")")))+
  coord_cartesian(ylim=c(0,2.5))+
  theme(legend.position = "none")

e3<-ggplot(commSummDF[commSummDF$treatment %in% 7:8 & commSummDF$landscape==3,],aes(x=effort,y=ex,color=factor(treatment),shape=factor(treatment)))+
  geom_line(size=1)+
  geom_point(size=2)+
  scale_colour_manual(name="Management\nstrategy",values=c("green","cyan"),labels=c("Assisted migration (threshold 50)","Assisted migration (threshold 75)"))+
  scale_shape_manual(name="Management\nstrategy",values=c("square","square"),labels=c("Assisted migration (threshold 50)","Assisted migration (threshold 75)"))+
  theme_classic()+
  scale_y_continuous(name="Average number of extinctions")+
  scale_x_continuous(name="Maximum relocations allowed")+
  coord_cartesian(ylim=c(0,2.5))+
  theme(legend.position = "none")

e4<-ggplot(commSummDF[commSummDF$treatment %in% 7:8 & commSummDF$landscape==5,],aes(x=effort,y=ex,color=factor(treatment),shape=factor(treatment)))+
  geom_line(size=1)+
  geom_point(size=2)+
  scale_colour_manual(name="Management\nstrategy",values=c("green","cyan"),labels=c("Assisted migration (threshold 50)","Assisted migration (threshold 75)"))+
  scale_shape_manual(name="Management\nstrategy",values=c("square","square"),labels=c("Assisted migration (threshold 50)","Assisted migration (threshold 75)"))+
  theme_classic()+
  scale_y_continuous(name="Average number of extinctions")+
  scale_x_continuous(name="Maximum relocations allowed")+
  coord_cartesian(ylim=c(0,2.5))+
  theme(legend.position = "none")

legs <- ggplot(commSummDF[commSummDF$treatment %in% c(1:3,5,7:8) & commSummDF$landscape==3 & commSummDF$effort,],aes(x=effort,y=ex,color=factor(treatment),shape=factor(treatment)))+
  geom_line(size=1)+
  geom_point(size=2)+
  scale_colour_manual(name="Management\nstrategy",values=c("red","blue","magenta","gray","green","cyan"),
                      labels=c("Corridor establishment","Stepping-stone reserves","Reinforce reserves",
                               "Restore everywhere","Assisted migration\n(threshold 50)","Assisted migration\n(threshold 75)"))+
  scale_shape_manual(name="Management\nstrategy",values=c("triangle","triangle","triangle","triangle","square","square"),
                      labels=c("Corridor establishment","Stepping-stone reserves","Reinforce reserves",
                               "Restore everywhere","Assisted migration\n(threshold 50)","Assisted migration\n(threshold 75)"))+
  theme_classic()+
  scale_y_continuous(name="Average number of extinctions")+
  scale_x_continuous(name=expression(paste("Total area restored (\u00D7",2^12,")")))+
  coord_cartesian(ylim=c(0,2.5))+
  theme(legend.title = element_text(size = 16),legend.text = element_text(size = 12),legend.box.margin = ggplot2::margin(0, 30, 30, 30))


plot_grid(e1,e2,e3,e4)


legs2 <- get_legend(legs)

t1 <- ggdraw() + 
  draw_label(
    "Few\nwide gaps",
    fontface = 'bold',
    x = 0,
    hjust = 0
  ) +
  theme(
    # add margin on the left of the drawing canvas,
    # so title is aligned with left edge of first plot
    plot.margin = margin(0, 0, 0, 7)
  )

t2 <- ggdraw() + 
  draw_label(
    "Many\nnarrow gaps",
    fontface = 'bold',
    x = 0,
    hjust = 0
  ) +
  theme(
    # add margin on the left of the drawing canvas,
    # so title is aligned with left edge of first plot
    plot.margin = margin(0, 0, 0, 7)
  )

eSquare<-plot_grid(e1,e3,e2,e4,align='hv',labels=c("a","b","c","d"),ncol=2)
eWords<-plot_grid(t1,legs2,t2,ncol=1,rel_heights = c(0.3,1,0.3))
eFull<-plot_grid(eWords,eSquare,rel_widths = c(0.6,2))

eFull

#### Figure 3 ####
# treatment is corridor, stepping stone, all locations, AM, or NOTHING
# 4 effort or 8 max relocations
# landscapes are not flat
commDF2 <- commDF[ ((commDF$treatment==1 & commDF$effort==4) | (commDF$treatment==8 & commDF$effort==8)| commDF$treatment==0) & commDF$landscape %in% c(3,5),]

# only action rows
# 2: AM
# 3: corridor
# 4: stepping stone
# 5: all locations
commDF2Action <- commDF2[c(seq(2,nrow(commDF2),by=3),seq(3,nrow(commDF2),by=3)),]

# only no action rows
commDF2NoAction <- commDF2[seq(1,nrow(commDF2),by=3),]

# How many extinctions when no action - how many extinctions with action
commDF2Action$extinctions3 <- commDF2NoAction$extinctions2-commDF2Action$extinctions2 

# Make a new set of SD categories
commDF3 <- commDF2Action
commDF3$regY <- (commDF3$tempYSD>0.125)+(commDF3$tempYSD>0.25)+(commDF3$tempYSD>0.375)+(commDF3$tempYSD>0.5)+(commDF3$tempYSD>0.625)+(commDF3$tempYSD>0.75)+(commDF3$tempYSD>0.875)
commDF3$regL <- (commDF3$tempLSD>0.25)+(commDF3$tempLSD>0.5)+(commDF3$tempLSD>0.75)+(commDF3$tempLSD>1)+(commDF3$tempLSD>1.25)+(commDF3$tempLSD>1.5)+(commDF3$tempLSD>1.75)

commSummDF2<-ddply(commDF3,c("landscape","treatment","effort","regY","regL"),summarize,N=length(community),
      gS=mean(gammaSimps2,na.rm=T),aS=mean(alphaSimps2,na.rm=T),
      gR=mean(gammaRich2,na.rm=T),aR=mean(alphaRich2,na.rm=T),
      tN=mean(totN2,na.rm=T),rM=mean(rangeMean2,na.rm=T),
      dN=mean(dispMean2,na.rm=T),tolM=mean(tolerMean2,na.rm=T),
      ex=mean(extinctions2,na.rm=T),g0=mean(gammaRich1,na.rm=T),
      d0=mean(dispMean1,na.rm=T),a0=mean(alphaSimps1,na.rm=T),
      ex3=mean(extinctions3,na.rm=T),exsd3=sd(extinctions3,na.rm=T),
      tysd=mean(tempYSD),tlsd=mean(tempLSD))


# test
cdm3 <- ggplot(commSummDF2[commSummDF2$landscape==3 & commSummDF2$treatment==1,],aes(x=regL/3.5,y=regY/7,fill=ex3))+
  geom_raster()+
  scale_fill_viridis_c(limits=c(0,2.3),option="turbo")+
  theme_classic()+
  scale_x_continuous(name="Local heterogeneity (H)",expand=c(0,0),breaks=c(0,1,2))+
  scale_y_continuous(name="Environmental stochasticity (S)",expand=c(0,0),breaks=c(0,0.5,1))+
  theme(aspect.ratio = 1)+
  coord_fixed()+
  theme(legend.position = "none")
cdm5 <- ggplot(commSummDF2[commSummDF2$landscape==5 & commSummDF2$treatment==1,],aes(x=regL/3.5,y=regY/7,fill=ex3))+
  geom_raster()+
  scale_fill_viridis_c(limits=c(0,2.3),option="turbo")+
  theme_classic()+
  scale_x_continuous(name="Local heterogeneity (H)",expand=c(0,0),breaks=c(0,1,2))+
  scale_y_continuous(name="Environmental stochasticity (S)",expand=c(0,0),breaks=c(0,0.5,1))+
  theme(aspect.ratio = 1)+
  coord_fixed()+
  theme(legend.position = "none")
amm3 <- ggplot(commSummDF2[commSummDF2$landscape==3 & commSummDF2$treatment==8,],aes(x=regL/3.5,y=regY/7,fill=ex3))+
  geom_raster()+
  scale_fill_viridis_c(limits=c(0,2.3),option="turbo")+
  theme_classic()+
  scale_x_continuous(name="Local heterogeneity (H)",expand=c(0,0),breaks=c(0,1,2))+
  scale_y_continuous(name="Environmental stochasticity (S)",expand=c(0,0),breaks=c(0,0.5,1))+
  theme(aspect.ratio = 1)+
  coord_fixed()+
  theme(legend.position = "none")
amm5 <- ggplot(commSummDF2[commSummDF2$landscape==5 & commSummDF2$treatment==8,],aes(x=regL/3.5,y=regY/7,fill=ex3))+
  geom_raster()+
  scale_fill_viridis_c(limits=c(0,2.3),option="turbo")+
  theme_classic()+
  scale_x_continuous(name="Local heterogeneity (H)",expand=c(0,0),breaks=c(0,1,2))+
  scale_y_continuous(name="Environmental stochasticity (S)",expand=c(0,0),breaks=c(0,0.5,1))+
  theme(aspect.ratio = 1)+
  coord_fixed()+
  theme(legend.position = "none")


legs <- ggplot(commSummDF2[commSummDF2$landscape==5 & commSummDF2$treatment==8,],aes(x=regL/3.5,y=regY/7,fill=ex3))+
  geom_raster()+
  scale_fill_viridis_c(name="Average number\nof extinctions\nprevented",limits=c(0,2.3),option="turbo")+
  theme_classic()+
  scale_x_continuous(name="Local heterogeneity (H)",expand=c(0,0),breaks=c(0,1,2))+
  scale_y_continuous(name="Environmental stochasticity (S)",expand=c(0,0),breaks=c(0,0.5,1))+
  theme(aspect.ratio = 1)+
  coord_fixed()

legs2 <- get_legend(legs)

t1 <- ggdraw() + 
  draw_label(
    "Few\nwide gaps",
    fontface = 'bold',
    x = 0,
    hjust = 0
  ) +
  theme(
    # add margin on the left of the drawing canvas,
    # so title is aligned with left edge of first plot
    plot.margin = margin(0, 0, 0, 7)
  )

t2 <- ggdraw() + 
  draw_label(
    "Many\nnarrow gaps",
    fontface = 'bold',
    x = 0,
    hjust = 0
  ) +
  theme(
    # add margin on the left of the drawing canvas,
    # so title is aligned with left edge of first plot
    plot.margin = margin(0, 0, 0, 7)
  )

u1 <- ggdraw() + 
  draw_label(
    "Corridor establishment\n(E=4WL)",
    fontface = 'bold',
    x = 0.5,
    hjust = 0.5
  ) +
  theme(
    # add margin on the left of the drawing canvas,
    # so title is aligned with left edge of first plot
    plot.margin = margin(0, 0, 0, 7)
  )

u2 <- ggdraw() + 
  draw_label(
    "Assisted migration\n(F=8)",
    fontface = 'bold',
    x = 0.5,
    hjust = 0.5
  ) +
  theme(
    # add margin on the left of the drawing canvas,
    # so title is aligned with left edge of first plot
    plot.margin = margin(0, 0, 0, 7)
  )



varSquare<-plot_grid(u1,u2,cdm3,amm3,cdm5,amm5,labels=c("","","a","b","c","d"),ncol=2,rel_heights = c(0.2,1,1))
varWords<-plot_grid(t1,legs2,t2,ncol=1,rel_heights = c(0.3,1,0.3))
varFull<-plot_grid(varWords,varSquare,rel_widths = c(0.3,2))
varFull
