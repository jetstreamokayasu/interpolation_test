require(TDA)
require(myfs)
require(rgl)

torus.300<-torusUnif(300, 1, 2.5)
figurePlot(torus.300)

torus.300.dist<-distance(torus.300)

#データ点1の近傍で実験
torus.vic1<-get.vicinity(torus.300.dist, 1, 15)

figurePlot.coloredVic(torus.300, torus.vic1, 1)

torus.vic1.line<-line.vics(1, torus.vic1)

vics.pca<-prcomp(torus.300[torus.vic1.line,])
plot(vics.pca[["x"]][,1], vics.pca[["x"]][,2], col=3, pch=16)
gridLine(vics.pca[["x"]], 4)

require(fields)
look<- image.count(vics.pca[["x"]][,1:2], nrow=4, ncol=4)
image.plot(look)

vic1s.plot.limit<-plotWithLimit(vics.pca[["x"]], 16, 3)

existCheck(c(-1.0026292, -0.7591696), c(-1.0026292+(0.4591493*1), -0.7591696+(1.450482*(1/4))), vics.pca[["x"]])
existCheck(c(-1.0026292, -0.7591696+(1.450482*(1/4))), c(-1.0026292+(0.4591493*1), -0.7591696+(1.450482*(2/4))), vics.pca[["x"]])

points(vics.pca[["x"]][8,1], vics.pca[["x"]][8,2], pch=16)
points(vics.pca[["x"]][1,1], vics.pca[["x"]][1,2], pch=16)

vic1s.pca.dist<-calcPcaDistance(vics.pca[["x"]])
vic1s.pca.dist2<-distance(vics.pca[["x"]])

plot(vics.pca[["x"]][,1], vics.pca[["x"]][,2], col=3, pch=16, xlim=c(-1.5, 1.5), ylim=c(-1.2, 1.2), asp = 1)
for (i in 1:length(vic1s.pca.dist)) {
  
  r<-vic1s.pca.dist[i]
  t = seq(0, 2*pi, length=100)
  plotWithLimit(cbind(r*cos(t), r*sin(t)), limit=vic1s.plot.limit)
  
}

r<-vic1s.pca.dist[1]
t = seq(0, 2*pi, length=100)
plotWithLimit(cbind(r*cos(t), r*sin(t)), limit=vic1s.plot.limit)

r2<-vic1s.pca.dist[2]
plotWithLimit(cbind(r2*cos(t), r2*sin(t)), limit=vic1s.plot.limit)


vic1.density<-circleDensity(vics.pca[["x"]][1,], vics.pca[["x"]][2,1:2], vic1s.pca.dist)
x<-seq(vic1s.plot.limit["x", "min"], vic1s.plot.limit["x", "max"], length=100)
y<-seq(vic1s.plot.limit["y", "min"], vic1s.plot.limit["y", "max"], length=100)

vic1.densityset<-calcDensitySet(vics.pca[["x"]][1,], pca.dist = vic1s.pca.dist, x, y)
image(x, y, vic1.densityset, col = terrain.colors(100))
persp(x, y, vic1.densityset, theta = 30, phi = 30, expand = 0.5, col = rainbow(50), border=NA)

vic1_2.dist<-get.vicinity(vic1s.pca.dist2, 2, length(vics.pca[["x"]][,1])-1)
vic1_2.densityset<-calcDensitySet(vics.pca[["x"]][2,], pca.dist = vic1_2.dist, x, y)
image(x, y, vic1_2.densityset+vic1.densityset, col = terrain.colors(100))

vic1.density.sum<-sumDensitySet(vics.pca[["x"]][,1:2], x, y)
image(x, y, vic1.density.sum, col = heat.colors(100))
persp(x, y, vic1.density.sum, theta = 30, phi = 30, expand = 0.5, col = heat.colors(100), border=NA)
contour(x, y, vic1.density.sum, method = "edge", vfont = c("sans serif", "plain"))

torus.vics1.pic<-pixelConvert(vics.pca[["x"]], 4)
test.pic<-matrix(0, 5, 5)

torus.vics1.cppic<-insertElement(torus.vics1.pic)

torus.vics1.incord<-pcaCoord.set(vics.pca[["x"]], torus.vics1.cppic, 4)
points(torus.vics1.incord, col=2, pch=16)

torus.vics1.oricord<-originCoodinate(vics.pca, torus.vics1.incord)
points3d(torus.vics1.oricord, col=2)

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
points3d(inter.oricord, col=2)
intered.torus.300<-conbineInterOrigin(torus.300, inter.oricord)
figurePlot(intered.torus.300)
ined.torus.300.diag<-ripsDiag(intered.torus.300, maxdimension = 2, maxscale = 3)
plot(ined.torus.300.diag[[1]])

#関数がうまくいくか試す用のトーラス状データセット
torus.collect11<- lapply(1:100, function(i){
  torus <- torusUnif(300, 1, 2.5)
  #cat("nsample=", nsample, "data", i, "\n")
  return(list(nsample = 300, noizyX = torus, diag = 0))
})
save(torus.collect11, file = "./data/torus.collect11")
figurePlot(torus.collect11[[1]][["noizyX"]])

for (i in 1:length(torus.collect11)) {
  
  inter.oricord<-interPolation_test(torus.collect11[[i]][["noizyX"]], 15, 4)
  torus.collect11[[i]][["noizyX"]]<-conbineInterOrigin(torus.collect11[[i]][["noizyX"]], inter.oricord)
  torus.collect11[[i]][["nsample"]]<-nrow(torus.collect11[[i]][["noizyX"]])
}

#トーラス状データのサイクル数推定
torus11.aggr<-homMethodsComp2compari3(torus.collect11, 2, 3, 10)
save(torus11.aggr, file = "./data/torus11.aggr")

torus11.cycle.dim1<-cyclenumber(torus11.aggr[[1]])
torus11.cycle.dim2<-cyclenumber(torus11.aggr[[2]])

torus.collect10.inted<-torus.collect10
save(torus.collect10.inted, file = "./data/torus.collect10.inted")

for (i in 1:length(torus.collect10.inted)) {
  
  inter.oricord<-interPolation_test(torus.collect10.inted[[i]][["noizyX"]], 15, 4)
  torus.collect10.inted[[i]][["noizyX"]]<-conbineInterOrigin(torus.collect10.inted[[i]][["noizyX"]], inter.oricord)
  torus.collect10.inted[[i]][["nsample"]]<-nrow(torus.collect10.inted[[i]][["noizyX"]])
}

torus10inted.aggr<-homMethodsComp2compari3(torus.collect10.inted, 2, 3, 10)
save(torus10inted.aggr, file = "./data/torus10inted.aggr")

torus10inted.dim1<-cyclenumber(torus10inted.aggr[[1]])
torus10inted.dim2<-cyclenumber(torus10inted.aggr[[2]])
