require(TDA)
require(myfs)
require(rgl)

torus.300<-torusUnif(300, 1, 2.5)
figurePlot(torus.300)

torus.dist<-distance(torus.300)

torus.vic1<-get.vicinity(torus.dist, 1, 15)

figurePlot.coloredVic(torus.300, torus.vic1, 1)

torus.vic1.line<-line.vics(1, torus.vic1)

vics.pca<-prcomp(torus.300[torus.vic1.line,])
plot(vics.pca[["x"]][,1], vics.pca[["x"]][,2], col=3, pch=16)
gridLine(vics.pca[["x"]], 4)

existCheck(c(-1.0026292, -0.7591696), c(-1.0026292+(0.4591493*1), -0.7591696+(1.450482*(1/4))), vics.pca[["x"]])
existCheck(c(-1.0026292, -0.7591696+(1.450482*(1/4))), c(-1.0026292+(0.4591493*1), -0.7591696+(1.450482*(2/4))), vics.pca[["x"]])

points(vics.pca[["x"]][8,1], vics.pca[["x"]][8,2], pch=16)
points(vics.pca[["x"]][1,1], vics.pca[["x"]][1,2], pch=16)

torus.vics1.pic<-pixelConvert(vics.pca[["x"]], 4)
test.pic<-matrix(0, 5, 5)

torus.vics1.cppic<-insertElement(torus.vics1.pic)

torus.vics1.incord<-pcaCoord.set(vics.pca[["x"]], torus.vics1.cppic, 4)
