#パッケージを使って書き換え
##中心点が作るボロノイ領域に隣接するボロノイ領域を作る点の中で、中心点に最も遠い点が作るボロノイ領域
#頂点に点を補う

require(TDA)
require(myfs)
require(rgl)
require(phacm)
require(pracma)
require(deldir)
require(ggplot2)
require(plyr)
require(reshape2)
require(ggmap)
require(tidyverse)
require(interpo3d)
require(seephacm)

torus300<-torusUnif(300, 1, 2.5)

torus300_dist<-dist(torus300)
trs300_vics1<-interpo3d:::get_vicinity(torus300_dist, 1, 15)
figurePlot3d(torus300[-trs300_vics1, ])
points3d(torus300[trs300_vics1, ], col=2)


trs300_pca1<-prcomp(torus300[trs300_vics1, ])
plot(trs300_pca1[["x"]][,1], trs300_pca1[["x"]][,2], col=3, pch=16)
points(trs300_pca1[["x"]][1,1], trs300_pca1[["x"]][1,2], col=4, pch=16)

trs300_res1<-deldir(trs300_pca1$x[,1], trs300_pca1$x[,2])
trs300_tiles1<-tile.list(trs300_res1)
for(i in 1:trs300_res1$n.data){	polygon(trs300_tiles1[[i]]) }


#中心点が作るボロノイ領域の頂点の中で、中心点から最も遠い頂点に補うように変更する
#変更。中心点からの平均以上の点のみに補う
trs300_vertx1<-cbind(trs300_tiles1[[1]][["x"]], trs300_tiles1[[1]][["y"]])
trs300_vertx1_dist<-dist(rbind(trs300_pca1[["x"]][1,1:2], trs300_vertx1))

trs300_idx1<-cbind(trs300_tiles1[[1]][["x"]], trs300_tiles1[[1]][["y"]]) %>% 
             rbind(trs300_pca1[["x"]][1,1:2], .) %>% dist() %>% as.matrix(.) %>% 
              max.col(.)

trs300_vertx1_dist2<-cbind(trs300_tiles1[[1]][["x"]], trs300_tiles1[[1]][["y"]]) %>% 
  rbind(trs300_pca1[["x"]][1,1:2], .) %>% dist() %>% as.matrix(.)

trs300_idx1_2nd<-which(trs300_vertx1_dist2[,1] >= mean(trs300_vertx1_dist2[2:8, 1]))

trs300_inp1<-cbind(trs300_tiles1[[1]][["x"]][trs300_idx1_2nd-1], trs300_tiles1[[1]][["y"]][trs300_idx1_2nd-1])
points(trs300_inp1, col=2, pch=16)
trs300_inpo1<-interpo3d:::origin_coordinate(trs300_pca1, trs300_inp1, torus300[1,])
points3d(trs300_inpo1, col=3)


#最も遠い頂点に補うようにした関数を試す
trs300_inps<-voronoi_interpo3(torus300, 15)

#中心点から遠い順に補うようにした関数を試す
trs300_inps2<-voronoiInterpo2(torus300, 15, 3)

trs300_idx1_3rd<-order(-trs300_vertx1_dist2[-1, 1])
trs300_inp1_2nd<-cbind(trs300_tiles1[[1]][["x"]][trs300_idx1_3rd[1:3]], trs300_tiles1[[1]][["y"]][trs300_idx1_3rd[1:3]])
points(trs300_inp1_2nd, col=2, pch=16)
