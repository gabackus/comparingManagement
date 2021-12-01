#### Initial ####
rm(list=ls())
require(plyr)
require(ggplot2)
require(cowplot)

#### Load data ####
load("specDFA.RData")
load("specDFB.RData")
specDF <- rbind(specDFA,specDFB)
rm(specDFA,specDFB)

#### Figure 4 ####
specSummDF<-ddply(specDF,c("landscape","treatment","effort","targ"),summarize,N=length(community),
            gS=mean(gammaSimps2,na.rm=T),aS=mean(alphaSimps2,na.rm=T),
            gR=mean(gammaRich2,na.rm=T),aR=mean(alphaRich2,na.rm=T),
            tN=mean(totN2,na.rm=T),rM=mean(rangeMean2,na.rm=T),
            dN=mean(dispMean2,na.rm=T),tolM=mean(tolerMean2,na.rm=T),
            ex=mean(extinctions2,na.rm=T),g0=mean(gammaRich1,na.rm=T),
            d0=mean(dispMean1,na.rm=T),a0=mean(alphaSimps1,na.rm=T),
            ex2=mean(extinction,na.rm=T),tte=mean(timeToExtinction,na.rm=T),
            N2=mean(NTarg2,na.rm=T),NMax=mean(TargMax,na.rm=T),
            NMin=mean(NTargMin,na.rm=T),NMean=mean(NTargMean,na.rm=T),
            pers=1-mean(extinction,na.rm=T))
specSummDF2<-specSummDF
#specSummDF2$exD<-specSummDF2[,6:23]-rbind(specSummDF[specSummDF$treatment==0,6:23],specSummDF[specSummDF$treatment==0,6:23],specSummDF[specSummDF$treatment==0,6:23])
specSummDF2$exD<-specSummDF2$pers-c(rep(specSummDF$pers[specSummDF$treatment==0 & specSummDF$landscape==3],3),rep(specSummDF$pers[specSummDF$treatment==0 & specSummDF$landscape==5],3))
specSummDF2$exPrint<-  paste(sprintf("%.2f",specSummDF2$exD),'decrease')



g1 <- ggplot(specSummDF[specSummDF$landscape==3 & specSummDF$treatment!=8,],aes(y=factor(-targ),x=1-pers,color=factor(effort),shape=factor(effort)))+
  geom_segment(data=specSummDF2[specSummDF2$landscape==3 & specSummDF2$treatment==1,],aes(xend=1-(pers-exD),yend=factor(-targ)),color="black",linetype=2)+
  geom_point(size=4)+
  geom_text(data=specSummDF2[specSummDF2$landscape==3 & specSummDF2$treatment==1,],aes(label=exPrint,x=1-(pers-exD/2)),nudge_y=0.25,color="black")+
  scale_x_continuous(name="Extinction likelihood")+
  scale_y_discrete(name="",labels=c("Random","Strongest equatorward\ncompetition","Strongest poleward\ncompetition","Narrowest\nthermal toerlance","Shortest\ndisperser"))+
  coord_cartesian(xlim=c(0,1))+
  scale_color_manual(name="Management\nstrategy",values=c("black","red","cyan"),labels=c("No action","Corridor","Assisted migration"),guide=F)+
  scale_shape_manual(name="Management\nstrategy",values=c("circle","triangle","square"),labels=c("No action","Corridor","Assisted migration"),guide=F)+
  theme_classic()

g2 <- ggplot(specSummDF[specSummDF$landscape==5 & specSummDF$treatment!=8,],aes(y=factor(-targ),x=1-pers,color=factor(effort),shape=factor(effort)))+
  geom_segment(data=specSummDF2[specSummDF2$landscape==5 & specSummDF2$treatment==1,],aes(xend=1-(pers-exD),yend=factor(-targ)),color="black",linetype=2)+
  geom_point(size=4)+
  geom_text(data=specSummDF2[specSummDF2$landscape==5 & specSummDF2$treatment==1,],aes(label=exPrint,x=1-(pers-exD/2)),nudge_y=0.25,color="black")+
  scale_x_continuous(name="Extinction likelihood")+
  scale_y_discrete(name="",labels=c("Random","Strongest equatorward\ncompetition","Strongest poleward\ncompetition","Narrowest\nthermal toerlance","Shortest\ndisperser"))+
  coord_cartesian(xlim=c(0,1))+
  scale_color_manual(name="Management\nstrategy",values=c("black","red","cyan"),labels=c("No action","Corridor","Assisted migration"),guide=F)+
  scale_shape_manual(name="Management\nstrategy",values=c("circle","triangle","square"),labels=c("No action","Corridor","Assisted migration"),guide=F)+
  theme_classic()

g3 <- ggplot(specSummDF[specSummDF$landscape==3 & specSummDF$treatment!=1,],aes(y=factor(-targ),x=1-pers,color=factor(effort),shape=factor(effort)))+
  geom_segment(data=specSummDF2[specSummDF2$landscape==3 & specSummDF2$treatment==8,],aes(xend=1-(pers-exD),yend=factor(-targ)),color="black",linetype=2)+
  geom_point(size=4)+
  geom_text(data=specSummDF2[specSummDF2$landscape==3 & specSummDF2$treatment==8,],aes(label=exPrint,x=1-(pers-exD/2)),nudge_y=0.25,color="black")+
  scale_x_continuous(name="Extinction likelihood")+
  scale_y_discrete(name="",labels=c("Random","Strongest equatorward\ncompetition","Strongest poleward\ncompetition","Narrowest\nthermal toerlance","Shortest\ndisperser"))+
  coord_cartesian(xlim=c(0,1))+
  scale_color_manual(name="Management\nstrategy",values=c("black","cyan"),labels=c("No action","Assisted migration"),guide=F)+
  scale_shape_manual(name="Management\nstrategy",values=c("circle","square"),labels=c("No action","Assisted migration"),guide=F)+
  theme_classic()

g4 <- ggplot(specSummDF[specSummDF$landscape==5 & specSummDF$treatment!=1,],aes(y=factor(-targ),x=1-pers,color=factor(effort),shape=factor(effort)))+
  geom_segment(data=specSummDF2[specSummDF2$landscape==5 & specSummDF2$treatment==8,],aes(xend=1-(pers-exD),yend=factor(-targ)),color="black",linetype=2)+
  geom_point(size=4)+
  geom_text(data=specSummDF2[specSummDF2$landscape==5 & specSummDF2$treatment==8,],aes(label=exPrint,x=1-(pers-exD/2)),nudge_y=0.25,color="black")+
  scale_x_continuous(name="Extinction likelihood")+
  scale_y_discrete(name="",labels=c("Random","Strongest equatorward\ncompetition","Strongest poleward\ncompetition","Narrowest\nthermal toerlance","Shortest\ndisperser"))+
  coord_cartesian(xlim=c(0,1))+
  scale_color_manual(name="Management\nstrategy",values=c("black","cyan"),labels=c("No action","Assisted migration"),guide=F)+
  scale_shape_manual(name="Management\nstrategy",values=c("circle","square"),labels=c("No action","Assisted migration"),guide=F)+
  theme_classic()

legs <- ggplot(specSummDF[specSummDF$landscape==3,],aes(y=factor(-targ),x=pers,color=factor(effort),shape=factor(effort)))+
  geom_point(size=4)+
  scale_x_continuous(name="Persistence likelihood")+
  scale_y_discrete(name="",labels=c("Random","Strongest equatorward\ncompetition","Strongest poleward\ncompetition","Narrowest\nthermal toerlance","Shortest\ndisperser"))+
  coord_cartesian(xlim=c(0,1))+
  scale_color_manual(name="Management strategy",values=c("black","red","cyan"),labels=c("No action\n","Corridor establishment\n(E=4WL)\n","Assisted migration\n(F=8)\n"))+
  scale_shape_manual(name="Management strategy",values=c("circle","triangle","square"),labels=c("No action\n","Corridor establishment\n(E=4WL)\n","Assisted migration\n(F=8)\n"))+
  theme_classic()+
  theme(legend.title = element_text(size = 14),legend.text = element_text(size = 10),legend.box.margin = ggplot2::margin(0, 30, 30, 30))
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

gSquare<-plot_grid(g1,g3,g2,g4,labels=c("a","b","c","d"), rel_heights = c( 1,1),ncol=2)
gWords<-plot_grid(t1,legs2,t2,ncol=1,rel_heights = c(0.3,1,0.3))
gFull<-plot_grid(gWords,gSquare,rel_widths = c(0.4,2))

gFull