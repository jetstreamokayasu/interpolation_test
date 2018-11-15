library(ade4)
data(rpjdl)
coa1 <- dudi.coa(rpjdl$fau, scannf = FALSE, nf = 4)
s.hist(coa1$li)
s.hist(coa1$li, cgrid = 2, cbr = 3, adj = 0.5, clab = 0)

coa<- dudi.coa(vics.pca[["x"]][,1:2], scannf = FALSE, nf = 4)


require(gregmisc)
h2d <- hist2d(vics.pca[["x"]][,1:2],show=FALSE, same.scale=TRUE, nbins=c(20,30))
 persp( h2d$x, h2d$y, h2d$counts,
               ticktype="detailed", theta=60, phi=30,
                 expand=0.5, shade=0.5, col="cyan", ltheta=-30)

 
#データ点1の近傍で実験
plot(vics.pca[["x"]][,1], vics.pca[["x"]][,2], col=3, pch=16)
gridLine(vics.pca[["x"]], 4)

vic1s.potential<-calcPotential(vics.pca[["x"]][,1:2], x, y)
image(x, y, vic1s.potential, col = terrain.colors(100))
persp(x, y, vic1s.potential, theta = 30, phi = 30, expand = 0.5, col = heat.colors(100), border=NA)
contour(x, y, vic1s.potential, method = "edge", vfont = c("sans serif", "plain"))
