##########TEPI

par(xpd=FALSE)
plot(c(0, 1), c(0, 1),type= "n", xlab = "Toxicity Probability", ylab = "Efficacy Probability",xaxt = "n",yaxt="n",cex.lab = 1.4)

axis(1, at=c(0,1), labels=c(0,1))
axis(2, at=c(0,1), labels=c(0,1))
axis(2, at=c(0.1,0.3,0.5,0.7), labels=c("Low","Moderate","High","Superb"),tck=0)
axis(1, at=c(0.1,0.3,0.7), labels=c("Low","Moderate High","Unacceptable"),tck=0)

rect(0, 0,0.2, 1, col = "green", border = "transparent") 
rect(0.2, 0,0.3, 0.6, col = "green", border = "transparent") 
rect(0.2, 0.6,0.3, 1, col = "turquoise", border = "transparent") 
rect(0.3, 0,0.4, 1, col = "turquoise", border = "transparent") 
rect(0.4, 0,1, 1, col = "pink", border = "transparent") 
rect(0.3, 0,0.4, 0.2, col = "pink", border = "transparent") 


abline(v  = 0.2,lty=2)
abline(v  = 0.3,lty=2)
abline(v  = 0.4,lty=2)
abline(v  = 0,lty=2)
abline(v  = 1,lty=2)

abline(h  = 0,lty=2)
abline(h  = 1,lty=2)
abline(h  = 0.2,lty=2)
abline(h  = 0.6,lty=2)
abline(h  = 0.4,lty=2)
par(xpd=TRUE)
legend(0.3,1.2, legend=c("E","S","D"), pch=rep(15,3), col=c("green","turquoise","pink"),ncol=3)

text(0.02,1.1,"(a) TEPI",cex=1.6)



########## PRINTE
par(xpd=FALSE)

plot(c(0, 1), c(0, 1),type= "n", xlab = "Toxicity Probability", ylab = "Efficacy Probability",xaxt = "n",yaxt="n",cex.lab=1.4)

axis(1, at=c(0,0.25,0.35,1), labels=c(0,expression(phi[T] - epsilon[1]),expression(paste("       ",phi[T] + epsilon[2])),1),cex.axis = 1.3)

axis(2, at=c(0,0.4,1), labels=c(0,expression(phi[E]),1),cex.axis=1.3)

rect(0, 0,0.35, 0.4, col = "green", border = "transparent") 

rect(0, 0.4,0.35, 1, col = "turquoise", border = "transparent") 
rect(0.35, 0,1, 1, col = "pink", border = "transparent") 

 
abline(h  = 0.0,lty=2)
abline(h  = 0.2,lty=2)
abline(h  = 0.4,lty=2)
abline(h  = 0.6,lty=2)
abline(h  = 0.8,lty=2)
abline(h  = 1.0,lty=2)

abline(v  = 0.0,lty=2)
abline(v  = 0.05,lty=2)
abline(v  = 0.15,lty=2)
abline(v  = 0.25,lty=2)
abline(v  = 0.35,lty=2)
abline(v  = 0.45,lty=2)
abline(v  = 0.55,lty=2)
abline(v  = 0.65,lty=2)
abline(v  = 0.75,lty=2)
abline(v  = 0.85,lty=2)
abline(v  = 0.95,lty=2)
abline(v  = 1,lty=2)

par(xpd=TRUE)

legend(0.3,1.2, legend=c("E","S","D"), pch=rep(15,3), col=c("green","turquoise","pink"),ncol=3)

text(0.05,1.1,"(b) PRINTE",cex=1.6)


