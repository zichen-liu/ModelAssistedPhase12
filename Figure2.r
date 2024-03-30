############# BOINET

setEPS()
postscript(paste0("results/Figure2_1.eps"))

par(xpd=FALSE)

plot(c(0, 1), c(0, 1),type= "n", xlab = "Toxicity Probability", ylab = "Efficacy Probability",xaxt = "n",yaxt="n",cex.lab=1.4,cex.main=1.4)

axis(1, at=c(0,0.3,0.7,1), labels=c(0,expression(lambda[e]),expression(paste(" ",lambda[d] )),1),cex.axis = 1.4)
axis(2, at=c(0,0.45,1), labels=c(0,expression(eta),1),cex.axis = 1.4)
rect(-1, -1,0.3, 0.45, col = "green", border = "transparent") 
rect(-1, 0.45,0.7, 2, col = "turquoise", border = "transparent") 
rect(0.7, -1,2, 2, col = "pink", border = "transparent") 
rect(0.3, -1,0.7, 0.45, col = "burlywood", border = "transparent") 


abline(h = 0.45)
abline(v  = 0.3)
abline(v=0.7)
text(0.13,0.2,"Escalate",cex=1.4)
text(0.13,0.7,"Stay",cex=1.4)
text(0.5,0.7,"Stay",cex=1.4)
text(0.5,0.2,"Choose from",cex=1.4)
text(0.5,0.1,"{ j-1, j, j+1 }",cex=1.4)
text(0.85,0.2,"De-escalate",cex=1.4)
text(0.85,0.7,"De-escalate",cex=1.4)
par(xpd=TRUE)
text(0.02,1.1,"(a) BOIN-ET",cex=1.6)

dev.off()

############# STEIN


setEPS()
postscript(paste0("results/Figure2_2.eps"))

par(xpd=FALSE)
plot(c(0, 1), c(0, 1),type= "n", xlab = "Toxicity Probability", ylab = "Efficacy Probability",xaxt = "n",yaxt="n",cex.lab=1.4,cex.main=1.4)

axis(1, at=c(0,0.3,0.7,1), labels=c(0,expression(lambda[e]),expression(paste(" ",lambda[d] )),1),cex.axis=1.4)
axis(2, at=c(0,0.45,1), labels=c(0,expression(eta),1),cex.axis=1.4)

rect(-1, -1,0.3, 0.45, col = "burlywood", border = "transparent") 
rect(-1, 0.45,0.7, 2, col = "turquoise", border = "transparent") 
rect(0.7, -1,2, 2, col = "pink", border = "transparent") 
rect(0.3, -1,0.7, 0.45, col = "orange", border = "transparent") 


abline(h = 0.45)
abline(v  = 0.3)
abline(v=0.7)
text(0.13,0.2,"Choose from",cex=1.4)
text(0.13,0.1,"{ j-1, j, j+1 }",cex=1.4)
text(0.13,0.7,"Stay",cex=1.4)
text(0.5,0.7,"Stay",cex=1.4)
text(0.5,0.2,"Choose from",cex=1.4)
text(0.5,0.1,"{ j-1, j }",cex=1.4)
text(0.85,0.2,"De-escalate",cex=1.4)
text(0.85,0.7,"De-escalate",cex=1.4)

par(xpd=TRUE)
text(0.02,1.1,"(b) STEIN",cex=1.6)

dev.off()
