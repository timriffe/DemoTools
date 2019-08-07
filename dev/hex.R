
# Author: tim
###############################################################################
# messing around w hex sticker idea.

# install.packages("hexSticker")
library(hexSticker)
library(DemoTools)

Age <- 0:99

V5 <- groupAges(pop1m_pasex, Age=Age)
Age5 <- as.integer(names(V5))
cf2 <- agesmth(Value = pop1m_pasex, 
	  Age = Age, 
	  method = "Carrier-Farrag", 
	  OAG = TRUE)
  st2 <- agesmth(Value = pop1m_pasex, 
		  Age = Age, 
		  method = "Strong", 
		  OAG = TRUE)

  plot(Age,pop1m_pasex,pch=16)
  lines(Age,splitUniform(V5,Age=Age5,OAG=FALSE), lty=2, lwd = 2)
  lines(Age,splitUniform(cf2,Age=Age5,OAG=FALSE),col="blue")
  lines(Age,splitUniform(st2,Age=Age5,OAG=FALSE),col="red")
  legend("topright",
		  pch=c(16,NA,NA,NA),
		  lty=c(NA,2,1,1),
		  col=c("black","black","blue","red"),
		  lwd=c(NA,2,1,1),
		  legend=c("orig single","orig 5","Carrier-Farrag","Strong"))
spr1 <- sprague(pop1m_pasex, Age=Age,OAG=FALSE)
spr2 <- sprague(cf2, Age=Age5,OAG=FALSE)
spr3 <- sprague(st2, Age=Age5,OAG=FALSE)

plot(Age,pop1m_pasex,pch=16, main = "",
		cex=.6,col = gray(.5),axes=FALSE,
		xlab="",ylab="")
lines(Age,spr1,lty=2,lwd=2)
lines(Age,spr2,col="blue",lwd=2)
lines(Age,spr3,col="red",lwd=2)

sticker(expression(plot(Age,pop1m_pasex,pch=16, main = "",
						cex=.5,col = gray(.5),axes=FALSE,
						xlab="",ylab=""),
lines(Age,spr1,lty=2,lwd=1.5),
lines(Age,spr2,col="blue",lwd=1.5),
lines(Age,spr3,col="red",lwd=1.5)),
		package="DemoTools", p_size=40, 
		p_family = "mono",
		s_x=.7, s_y=.7, s_width=2.5, s_height=2.5,
		h_fill = "#FFFFFF",
		h_color = "royalblue",
		h_size = 1.2,
		filename="DemoTools.png",
		p_color = "black",
		white_around_sticker = TRUE,
		dpi=600)
args(sticker)


# save svg blank:
x_vertices <- 1+c(rep(-sqrt(3)/2, 2), 0, rep(sqrt(3)/2, 2), 0)
y_vertices <- 1+c(0.5, -0.5, -1, -0.5, 0.5, 1)
range(x_vertices);range(y_vertices)

pdf("DemoToolsBorder.pdf")
par(mai=c(0,0,0,0),xaxs="i",yaxs="i")
plot(NULL, type = 'n', asp=1,xlim=c(-.1,2.1),ylim=c(-.1,2.1),
		axes=FALSE, xlab="",ylab="")
polygon(x_vertices, y_vertices,lwd=3)
dev.off()

getwd()
pdf("DemoToolsInner.pdf")
par(mai=c(0,0,0,0),xaxs="i",yaxs="i")
plot(Age,pop1m_pasex,pch=16, main = "",
		cex=1,col = gray(.5),axes=FALSE,
		xlab="",ylab="")
lines(Age,spr1,lty=2,lwd=1.5)
lines(Age,spr2,col="blue",lwd=1.5)
lines(Age,spr3,col="red",lwd=1.5)
dev.off()


barplot(pop1m_pasex,space=0,border=NA)

library(RcolorBrewer)
display.brewer(11,"BrBG")

library(HMDHFDplus)
fem<-readHMDweb("USA","fltper_1x1",password=Sys.getenv("pw"),username=Sys.getenv("us"))
mal <-readHMDweb("USA","mltper_1x1",password=Sys.getenv("pw"),username=Sys.getenv("us"))

qxf <- subset(fem, Year == 2017 & Age <= 85)$qx
qxm <- subset(mal, Year == 2017 & Age <= 85)$qx

pdf("qx.pdf")
par(mai = c(0,0,0,0),xaxs='i',yaxs='i')
plot(NULL,type='l',log='y',
	 axes=FALSE, xlab="",ylab="",xlim=c(0,85),ylim=c(0.000001,.8), lwd = 2)
polygon(c(0:85,85:0),c(qxf,rev(qxm)))
dev.off()
getwd()


# get bars
pop <-readHMDweb("USA","Population",password=Sys.getenv("pw"),username=Sys.getenv("us"))
m2030 <- subset(pop, Year == 2000 & Age >= 20 & Age <= 29)$Male1
barplot(m2030,space=0,width=rep(1,length(m2030)),col = "#eeb820ff",border="#e08819ff")
library(DemoTools)
barplot(pop1m_pasex[21:30],space=0,width=rep(1,length(m2030)),col = "#eeb820ff",border="#e08819ff")
pop1m_pasex
