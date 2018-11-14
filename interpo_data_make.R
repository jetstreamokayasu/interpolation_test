require(TDA)
require(myfs)
require(rgl)

torus.300<-torusUnif(300, 1, 2.5)
figurePlot(torus.300)

rgl.snapshot("./data/torus_300_1.png")

torus.300.dist<-distance(torus.300)


#データ点1の近傍で実験
torus.vic1<-get.vicinity(torus.dist, 1, 15)

figurePlot.coloredVic(torus.300, torus.vic1, 1)
rgl.snapshot("./data/torus_300.png")
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
points(torus.vics1.incord, pch=16, col=2)

torus.vics1.oricord<-originCoodinate(vics.pca, torus.vics1.incord)
points3d(torus.vics1.oricord, col=2)
rgl.snapshot("./data/torus_300_intered.png")
#データ点17の近傍で実験
torus.vic17<-get.vicinity(torus.dist, 17, 15)

figurePlot.coloredVic(torus.300, torus.vic17, centr = 17)

torus.vic17.line<-line.vics(centr =17, torus.vic17)

vic17s.pca<-prcomp(torus.300[torus.vic17.line,])
plot(vic17s.pca[["x"]][,1], vic17s.pca[["x"]][,2], col=3, pch=16)
gridLine(vic17s.pca[["x"]], 4)

torus.vic17s.pic<-pixelConvert(vic17s.pca[["x"]], 4)

torus.vic17s.cppic<-insertElement(torus.vic17s.pic)

torus.vic17s.incord<-pcaCoord.set(vic17s.pca[["x"]], torus.vic17s.cppic, 4)
points(torus.vic17s.incord, col=2, pch=16)

torus.vic17s.oricord<-originCoodinate(vic17s.pca, torus.vic17s.incord)
points3d(torus.vic17s.oricord, col=2)


torus.vic210<-get.vicinity(torus.dist, 210, 15)
torus.vic210.line<-line.vics(210, torus.vic210)

torus.covic210<-coveredVic(torus.vic210.line, torus.300, 5)

#データ点100の近傍で実験
torus.vic100s<-get.vicinity(torus.dist, 100, 15)

torus.vic100s.line<-line.vics(centr =100, torus.vic100s)

vic100s.pca<-prcomp(torus.300[torus.vic100s.line,])
plot(vic100s.pca[["x"]][,1], vic100s.pca[["x"]][,2], col=3, pch=16)
gridLine(vic100s.pca[["x"]], 4)

torus.vic100s.pic<-pixelConvert(vic100s.pca[["x"]], 4)

torus.vic100s.cppic<-insertElement(torus.vic100s.pic)

torus.vic100s.incord<-pcaCoord.set(vic100s.pca[["x"]], torus.vic100s.cppic, 4)
points(torus.vic100s.incord, col=2, pch=16)

torus.vic100s.oricord<-expandProcess(torus.vic100s, torus.vic100s.line, torus.300, torus.dist, 4)

figurePlot.coloredVic(torus.300, torus.vic100s, centr =100)
points3d(torus.vic100s.oricord, col=2)
##########################

#データ補間テスト
inter.oricord<-interPolation_test(torus.300, 15, 4)
figurePlot(torus.300)

points3d(inter.oricord, col="orange")
rgl.snapshot("./data/torus_300_intered_fin.png")

inter.oricord20<-interPolation_test(torus.300, 20, 4)
figurePlot(torus.300)
points3d(inter.oricord20, col="orange")

sphere<-sphereUnif(200, 2, 1)
plot3d(sphere)
sphere.inoricord<-interPolation_test(sphere, 15, 4)
points3d(sphere.inoricord, col="orange")
sphere.inoricord2<-interPolation_test(sphere, 10, 3)

torus.320<-torusUnif(320, 1, 2.5)
figurePlot(torus.320)
torus320.dist<-distance(torus.320)
torus320.vic20s<-get.vicinity(torus320.dist, center = 20, nvic = 20)
torus320.vic20s.line<-line.vics(20, torus320.vic20s)
figurePlot.coloredVic(torus.320, torus320.vic20s, 20)
figurePlot(torus.320[-torus320.vic20s.line, ])
torus320.no20<-torus.320[-torus320.vic20s.line, ]

torus320.no20.incord<-interPolation_test(torus320.no20, 15, 4)
points3d(torus320.no20.incord, col="orange")

#サイクル/ノイズ判別閾値テスト
sphere.400<-sphereUnif(400, 2, 1)
plot3d(sphere.400)
pre_thresh<-meanVicsDestance(sphere.400, 15)

#近傍3点の平均による補間
torus300.inter<-meanInterPolation(torus.300, 2)
figurePlot(torus.300)
points3d(torus300.inter, col="orange")
