require(TDA)
require(myfs)
require(rgl)
require(deldir)
require(ggplot2)
require(plyr)
require(reshape2)
require(ggmap)

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

#ボロノイ図試し
plot(vics.pca[["x"]][,1], vics.pca[["x"]][,2], col=3, pch=16)
res1<-deldir(vics.pca$x[,1], vics.pca$x[,2])
tiles <- tile.list(res1)
plot(range(vics.pca[["x"]][,1]),range(vics.pca[["x"]][,2]),type="n")
for(i in 1:res1$n.data){	polygon(tiles[[i]]) }
points(vics.pca[["x"]][,1:2], col=3, pch=16)
points(vics.pca[["x"]][1,1], vics.pca[["x"]][1,2], col=2, pch=16)
#text(d$longitude,d$latitude+0.0005,d$id,col=as.numeric(d$type)+1)
points(tiles[[1]][["x"]], tiles[[1]][["y"]], pch=15, col=2)
text(tiles[[1]][["x"]], tiles[[1]][["y"]], c(1:length(tiles[[1]][["x"]])))
points(tiles[[1]][["x"]], tiles[[1]][["y"]], pch=15, col=2)
points(vics.pca[["x"]][6,1], vics.pca[["x"]][6,2], col=5, pch=17)
points(tiles[[6]][["x"]], tiles[[6]][["y"]], pch=15, col=2)
text(tiles[[3]][["x"]], tiles[[3]][["y"]], c(1:length(tiles[[3]][["x"]])))

#ボロノイ領域内にランダムに点を打つ
set.seed(100)
ranx<-runif(1, min = min(tiles[[1]][["x"]]), max = max(tiles[[1]][["x"]]))
rany<-runif(1, min = min(tiles[[1]][["y"]]), max = max(tiles[[1]][["y"]]))
points(ranx, rany, pch=13, col=4)
cross.mem<-which(tiles[[1]][["x"]]>=ranx)
sides1<-vertex.side(tiles[[1]], 1)
sides<-sapply(cross.mem, function(k)vertex.side(tiles[[1]], k))
#ver1.set<-vertexSet(tiles[[1]], sides)
cross1.side<-sidesSet(sides)
hline<-matrix(c(ranx, rany, max(tiles[[1]][["x"]][which(tiles[[1]][["x"]]>=ranx)]), rany), 2, 2, byrow=T)
ncross1<-crossCheck(tiles[[1]], hline, t(cross1.side))
ncross1.1<-crossCheck(tiles[[1]], hline, sides)

ranpoint1<-randomPointVoronoi(tiles[[1]])
points(ranpoint1[1], ranpoint1[2], pch=13, col=4)

neibor1<-neighbourVoronoi(tiles, 1)
points(vics.pca[["x"]][2,1], vics.pca[["x"]][2,2], col=2, pch=17)
points(vics.pca[["x"]][3,1], vics.pca[["x"]][3,2], col=2, pch=17)
points(vics.pca[["x"]][4,1], vics.pca[["x"]][4,2], col=2, pch=17)
points(vics.pca[["x"]][7,1], vics.pca[["x"]][7,2], col=2, pch=17)
ranpoint2<-randomPointVoronoi(tiles[[2]])
points(ranpoint2[1], ranpoint2[2], pch=13, col=4)
ranpoint3<-randomPointVoronoi(tiles[[3]])
points(ranpoint3[1], ranpoint3[2], pch=13, col=4)


voron.oricord<-voronoiProcess(torus.vic1.line, torus.300)
points3d(voron.oricord, col=2)

figurePlot(torus.300)
in.oricord.vo<-voronoiInterpo(torus.300, 10)
points3d(in.oricord.vo, col=2)

#データ点17の近傍で実験
torus.vic17<-get.vicinity(torus.dist, 17, 15)
figurePlot.coloredVic(torus.300, torus.vic17, centr = 17)
torus.vic17.line<-line.vics(centr =17, torus.vic17)
vic17s.pca<-prcomp(torus.300[torus.vic17.line,])
plot(vic17s.pca[["x"]][,1], vic17s.pca[["x"]][,2], col=3, pch=16)

res17<-deldir(vic17s.pca$x[,1], vic17s.pca$x[,2])
tile17s <- tile.list(res17)
for(i in 1:res17$n.data){	polygon(tile17s[[i]]) }
points(vic17s.pca[["x"]][1,1], vic17s.pca[["x"]][1,2], col=2, pch=16)
neibor17<-neighbourVoronoi(tile17s, 1)
cenpoints17<-sapply(c(1,neibor17), function(k)centerVoronoi(tile17s[[k]]))
for(i in 1:(length(neibor17)+1)){
  points(cenpoints17[1,i], cenpoints17[2,i], pch=13, col=4)
}

voron17.oricord<-voronoiProcess(torus.vic17.line, torus.300)
points3d(voron17.oricord, col=2)

#ボロノイ補間がうまくいくか試す用のトーラス状データセット
torus.collect14 <- lapply(1:100, function(i){
  nsample <- 300
  #var <- runif(1, var.min, var.max)
  #noize.torus <- matrix(rnorm(nsample * 3, 0, var), nrow = nsample)
  torus <- torusUnif(nsample, 1, 2.5)
  return(list(nsample = nsample, noizyX = torus, diag = 0))
})
save(torus.collect14, file = "./data/torus.collect14")
torus.collect14.inted<-torus.collect14
for (i in 1:length(torus.collect14)) {
  inter.oricord<-voronoiInterpo(torus.collect14[[i]][["noizyX"]], 15)
  torus.collect14.inted[[i]][["noizyX"]]<-conbineInterOrigin(torus.collect14[[i]][["noizyX"]], inter.oricord)
  torus.collect14.inted[[i]][["nsample"]]<-nrow(torus.collect14[[i]][["noizyX"]])
}
save(torus.collect14.inted, file = "./data/torus.collect14.inted")
torus14inted.aggr<-proposedMethodOnly(torus.collect14.inted, 2, 3, 10)

save(torus14inted.aggr, file = "./data/torus14inted.aggr")

torus14inted.dim1<-cyclenumber(torus14inted.aggr[[1]], compare=F)
torus14inted.dim2<-cyclenumber(torus14inted.aggr[[2]], compare=F)

figurePlot(torus.collect14[[1]][["noizyX"]])
inter14.1.oricord<-voronoiInterpo(torus.collect14[[1]][["noizyX"]], 15)
points3d(inter14.1.oricord, col="orange")
figurePlot(torus.collect14.inted[[10]][["noizyX"]])

torus14inted1.diag<-ripsDiag(torus.collect14.inted[[1]][["noizyX"]], 2, 3, printProgress = T)
plot(torus14inted1.diag[[1]])

torus14.aggr<-proposedMethodOnly(torus.collect14, 2, 3, 10)
torus14.dim1<-cyclenumber(torus14.aggr[[1]], compare=F)
torus14.dim2<-cyclenumber(torus14.aggr[[2]], compare=F)
