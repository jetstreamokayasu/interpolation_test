plot.circle <- function(x, y, r){
  theta <- seq(-pi, pi, length=100)
  polygon(x + r*cos(theta), y + r*sin(theta), col="gray")
}

xp<-c(-1, 2, -3, 3, -2, 2)
yp<-c(4, 4, 2, 2, 0, 0)
plot(xp, yp, xlim=c(-5, 5), ylim=c(-2, 6), type="n")
sapply(1:6, function(k)plot.circle(xp[k], yp[k], sqrt(17)/2))
points(xp, yp, pch=16, col=2, cex=2)

resp<-deldir(xp, yp)
p.tiles <- tile.list(resp)
plot(p.tiles)

for(i in 1:6){polygon(p.tiles[[i]], lwd=2)}

plot(xp, yp, xlim=c(0,2), ylim=c(0,4), type="n", xlab="", ylab="")
lines(c(0, 2.5-sqrt(2)), c(0, 0), lwd=2, col=4)
lines(c(0, 2.5-2), c(1, 1), lwd=2, col=4)

#スイスロールを作る
#可視化及び後に用いる次元削除のためRDRToolboxのインストール
source("http://bioconductor.org/biocLite.R")
biocLite("RDRToolbox")

library(RDRToolbox)

N=1000
p = sqrt(2 + 2 * seq(-1, 1 - 2/N, 2/N))
y = 2 * runif(N, -1, 1)
d_sr = cbind(p * cos(2*pi*p), y, p * sin(2*pi*p))

color <- vector()
for(i in 1:N){
  color[i] = rainbow(ceiling(N/50))[ceiling(i/50)]
}
labels=seq(N)

plotDR(data=d_sr, labels=labels, col=color)
rgl.snapshot("./data/swisroll.png")
