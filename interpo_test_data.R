library(ade4)
require(lle)
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


#LLEを使ってみる
kneibor<-calc_k(torus.300, 2)
torus300.vic1<-get.vicinity(torus.300.dist, 1, 10)
vic1s.line<-line.vics(1, torus300.vic1)
torus300.vic1s<-lle(torus.300[vic1s.line, ], 2, 10)
plot(torus300.vic1s[["Y"]])
figurePlot.coloredVic(torus.300, torus300.vic1, 1)

#save2File関数試し
#現在失敗中
save2File(torus.300, dir="./")

#PCAで貼られた面を表示
torus1.coefs<-confirmPlane(vics.pca)
figurePlot.coloredVic(torus.300, torus.vic1, 1)
plot3d(torus.300[-torus.vic1.line, ], size=4)
aspect3d("iso")
planes3d(torus1.coefs[1], torus1.coefs[2], torus1.coefs[3], torus1.coefs[4], col="blue", alpha=0.5)
points3d(torus.300[torus.vic1.line, ], col=3, size=4)
points3d(torus.vics1.oricord, col=2)
rgl.snapshot("./torus_300_1_plane.png")

figurePlot.coloredVic(torus.300, torus.vic17, centr = 17)
torus17.coefs<-confirmPlane(vic17s.pca)
planes3d(torus17.coefs[1], torus17.coefs[2], torus17.coefs[3], torus17.coefs[4], col="blue", alpha=0.5)
points3d(torus.vic17s.oricord, col=2)
