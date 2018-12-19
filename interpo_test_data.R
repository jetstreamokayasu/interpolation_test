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
rgl.snapshot("./torus_300_1_plane.eps")

figurePlot.coloredVic(torus.300, torus.vic17, centr = 17)
torus17.coefs<-confirmPlane(vic17s.pca)
planes3d(torus17.coefs[1], torus17.coefs[2], torus17.coefs[3], torus17.coefs[4], col="blue", alpha=0.5)
points3d(torus.vic17s.oricord, col=2)


#PCAで写された点で作る凸包内に補間点が有るかどうかを判定
plot(vics.pca[["x"]][,1], vics.pca[["x"]][,2], col=3, pch=16, cex=2)
res1<-deldir(vics.pca$x[,1], vics.pca$x[,2])
tiles <- tile.list(res1)
plot(range(vics.pca[["x"]][,1]),range(vics.pca[["x"]][,2]),type="n")
for(i in 1:res1$n.data){	polygon(tiles[[i]], lwd=2) }
points(vics.pca[["x"]][1,1], vics.pca[["x"]][1,2], col=2, pch=16)
chul<-chull(vics.pca[["x"]][,1:2])
polygon(vics.pca[["x"]][chul,1:2])
text(vics.pca[["x"]][chul,1], vics.pca[["x"]][chul,2], chul)

tst.x<--0.5
tst.y<-0.5
points(tst.x, tst.y, pch=13, col=4)
cross.mem<-chul[which(vics.pca[["x"]][chul,1]>=tst.x)]
sides1<-convex_hull_vertx(chul, 16)
sides<-sapply(cross.mem, function(k)convex_hull_vertx(chul, k))
cross1.side<-sidesSet(sides)
hline<-matrix(c(tst.x, tst.y, max(vics.pca[["x"]][chul,1][which(vics.pca[["x"]][chul,1]>=tst.x)]), tst.y), 2, 2, byrow=T)
c.ncross1<-convex_hull_check(vics.pca, hline, t(cross1.side))

p <- ggplot(vics.pca[["x"]][,1:2]) +   geom_point()

exist1<-exist_convexhull_check(vics.pca, t(as.matrix(c(tst.x, tst.y))))
tstxy<-matrix(c(-0.5, 0.5, 0.5, -0.5), 2,2, byrow = T)
exist1a<-exist_convexhull_check(vics.pca, tstxy)

voron.oricord_B<-voronoiBorder(torus.vic1.line, torus.300)
points(voron.oricord_B[[2]], pch=16, col=2)
points3d(voron.oricord_B, col=2)

#17点の近傍で試し
voron.oricord_17<-voronoiBorder(torus.vic17.line, torus.300)
points(voron.oricord_17[[2]], pch=16, col=2)
points(tile17s[[1]][["x"]], tile17s[[1]][["y"]], pch=16, col=2, cex=2)
points3d(voron.oricord_17[[1]], col=2)
