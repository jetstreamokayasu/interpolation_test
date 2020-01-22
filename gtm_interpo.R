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
require(gtm)

#GTMを用いた補間

torus300<-torusUnif(300, 1, 2.5)

torus300_dist<-dist(torus300)
trs300_vics1<-interpo3d:::get_vicinity(torus300_dist, 1, 30)
figurePlot3d(torus300[-trs300_vics1, ])
points3d(torus300[trs300_vics1, ], col=2)
spheres3d(rbind(c(0, 0, 0), torus300[trs300_vics1[1], ]), radius = 0.1, col=3)

X<-torus300[trs300_vics1, ]
# (分散が0の変数を削除した後に) 1. オートスケーリング
Var0Variable <- which(apply(X,2,var) == 0)
if (length(Var0Variable) == 0) {
  print("分散が0の変数はありません")
} else {
  sprintf("分散が0の変数が %d つありました", length(Var0Variable))
  print( "変数:" )
  print( Var0Variable )
  print( "これらを削除します" )
  X <- X[,-Var0Variable]
}
X <- scale(X, center = TRUE, scale = TRUE)
# 2. マップサイズ
MapsizeColumn = 15 #横 10
MapsizeRow = 15 #縦 10 、
# 3. 動径基底関数 (Radial Basis Function, RBF) の数
RBFsizeColumn = 3 #横 3
RBFsizeRow = 3 #縦 3
# 4. ガウス関数の分散
RBFVariance = 1
# 5. EMアルゴリズムにおけるパラメータλの値
Lambda = 0
# 6. データ空間における分散の逆数βの初期値
# 7. 重みWの初期値
XGrid = gtm.rctg( MapsizeColumn, MapsizeRow)#写像先のグリッド
RBFGrid = gtm.rctg( RBFsizeColumn, RBFsizeRow)#基底関数のグリッド
RBFSetup = gtm.gbf( RBFGrid, RBFVariance^(1/2), XGrid)
InitnalWBeta = gtm.pci.beta( X, XGrid, RBFSetup)
#Beta = 0.01
Beta = InitnalWBeta$beta
# 8. 学習回数
NumOfTraining = 100 #Training cycleとlog-likelihoodとのプロットが収束していないようでしたら学習回数を大きくしてください
# 9. GTMマップ作成
GTMResults = gtm.trn( X, RBFSetup, InitnalWBeta$W, Lambda, NumOfTraining, Beta)
# 10. 二次元のマップ上でサンプルの位置関係を確認
GTMDist = gtm.dist(X, RBFSetup %*% GTMResults$W)
minDist = apply(GTMDist, 2, which.min)
GTMMode = matrix(XGrid[minDist, ], ncol = 2)
GTMR = gtm.resp3(GTMDist, GTMResults$beta, ncol(X))$R
GTMMean = t(GTMR) %*% XGrid
par(pty = "s")
plot( GTMMode[,1], GTMMode[,2], xlim=c(-1.2, 1.2), ylim=c(-1.2, 1.2), col = 'blue', xlab = "z1", ylab = "z2")
text( GTMMode[,1], GTMMode[,2], labels = rownames(X), pos=3, offset = 0.1)
par(pty = "m")
par(pty = "s")
plot( GTMMean[,1], GTMMean[,2], xlim=c(-1.2, 1.2), ylim=c(-1.2, 1.2), col = 'green', xlab = "z1", ylab = "z2", pch=16, cex=2)
text( GTMMean[,1], GTMMean[,2], labels = 1:72, pos=3, offset = 0.1)
par(pty = "m")

#次元削減後のデータ点にラベルを付ける
library(maptools)
pointLabel(GTMMean[,1], GTMMean[,2], as.character(trs300_vics1))
pointLabel(GTMMean[,1], GTMMean[,2], as.character(1:31))

#逆写像を求める
ori_GTMMean = t(GTMR) %*% (RBFSetup %*% GTMResults$W)
##スケーリングを戻す
inv_GTMMean<-sapply(1:nrow(ori_GTMMean), function(k){ori_GTMMean[k, ] * attr(X, "scaled:scale") + attr(X, "scaled:center")})
inv_GTMMean<-t(inv_GTMMean)

points3d(inv_GTMMean, col=3)


#次元削減後にボロノイ分割
points(GTMMean[1,1], GTMMean[1,2], pch=16, col=4, cex=2)
trs300_res1<-deldir(GTMMean[,1], GTMMean[,2])
trs300_tiles1<-tile.list(trs300_res1)
for(i in 1:trs300_res1$n.data){	polygon(trs300_tiles1[[i]], lwd=2) }

trs300_vertx1<-cbind(trs300_tiles1[[1]][["x"]], trs300_tiles1[[1]][["y"]])
points(trs300_tiles1[[1]][["x"]], trs300_tiles1[[1]][["y"]], pch=16, col=2, cex=2)
points(nei_interpo, pch=16, col=2, cex=2)
pointLabel(trs300_tiles1[[1]][["x"]], trs300_tiles1[[1]][["y"]], as.character(1:31))

RBF_inter = gtm.gbf( RBFGrid, RBFVariance^(1/2), cbind(trs300_tiles1[[1]][["x"]], trs300_tiles1[[1]][["y"]]))
inter_dist<-gtm.dist(X, RBF_inter %*% GTMResults$W)
inter_R = gtm.resp3(inter_dist, GTMResults$beta, ncol(X))$R
inter_mean = t(inter_R) %*% (RBF_inter %*% GTMResults$W)
inter_inv<-sapply(1:nrow(inter_mean), function(k){inter_mean[k, ] * attr(X, "scaled:scale") + attr(X, "scaled:center")})
inter_inv<-t(inter_inv)
points3d(inter_inv, col=4)


#中心点のボロノイ領域に隣接するボロノイ領域を見つける
equal<-trs300_tiles1[[1]][["x"]] %in% trs300_tiles1[[2]][["x"]]
x_eq<-t(sapply(trs300_tiles1[[1]][["x"]][equal], function(x){x==trs300_tiles1[[2]][["x"]]}))
y_eq<-sapply(1:length(trs300_tiles1[[1]][["x"]][equal]), function(k){trs300_tiles1[[2]][["y"]][x_eq[k,]]==trs300_tiles1[[1]][["y"]][equal][k]})

tst<-neibor_voronoi(trs300_tiles1)
nei_tile_list<-which(unlist(tst)==T)+1

nei_inter<-lapply(trs300_tiles1[c(1, nei_tile_list)], function(tile){cbind(tile[["x"]], tile[["y"]])})
nei_interpo<-unique(do.call(rbind, nei_inter))

RBF_inter2 = gtm.gbf( RBFGrid, RBFVariance^(1/2), nei_interpo)
inter_dist2<-gtm.dist(X, RBF_inter2 %*% GTMResults$W)
inter_R2 = gtm.resp3(inter_dist2, GTMResults$beta, ncol(X))$R
inter_mean2 = t(inter_R2) %*% (RBF_inter2 %*% GTMResults$W)
inter_mean2<- RBF_inter2 %*% GTMResults$W
inter_inv2<-sapply(1:nrow(inter_mean2), function(k){inter_mean2[k, ] * attr(X, "scaled:scale") + attr(X, "scaled:center")})
inter_inv2<-t(inter_inv2)
points3d(inter_inv2, col=4)

#sphere3dで描画
spheres3d(torus300[-trs300_vics1, ], radius = 0.08)
aspect3d("iso")
spheres3d(torus300[trs300_vics1, ], col=3, radius = 0.08)
spheres3d(inter_inv2, col=2, radius = 0.08)

rgl.snapshot("torus300_vics.png") 

#GTM全体補間試し
#近傍30点
trs300_incolle_set1<-gtm_interpolate(torus300_colle_set[[1]][1:5], 30)

figurePlot3d(trs300_incolle_set1[[1]][["noizyX"]][1:300, ])
points3d(trs300_incolle_set1[[1]][["noizyX"]][301:755, ], col=2)
spheres3d(trs300_incolle_set1[[1]][["noizyX"]][1:300, ], radius = 0.08)
spheres3d(trs300_incolle_set1[[1]][["noizyX"]][301:889, ], radius = 0.08, col=2)

rgl.snapshot("torus300.png") 
rgl.snapshot("torus300_ind_all.png") 

in_trs300_1_w1_5errs<-lapply(trs300_incolle_set1, function(X)torus_disterror(X[["noizyX"]], maxr = 2.5, minr = 1, nps = 300))
par(mgp=c(2.4,1,0))
boxplot(in_trs300_1_w1_5errs, xlab="Data Set", ylab="Error", cex.lab=1.6, cex.axis=1.6, lwd=2)

#誤差比較
trs300_incolle_set2<-all_interpolate(torus300_colle_set[[1]][1:5], 15)
in_trs300_2_w1_5errs<-lapply(trs300_incolle_set2, function(X)torus_disterror(X[["noizyX"]], maxr = 2.5, minr = 1, nps = 300))
par(mgp=c(2.4,1,0))
boxplot(in_trs300_2_w1_5errs, xlab="Data Set", ylab="Error", cex.lab=1.6, cex.axis=1.6, lwd=2)

##近傍50点
trs300_incolle_set1b<-gtm_interpolate(torus300_colle_set[[1]][1:5], 50)

figurePlot3d(trs300_incolle_set1b[[1]][["noizyX"]][1:300, ])
points3d(trs300_incolle_set1b[[1]][["noizyX"]][301:608, ], col=2)

in_trs300b_1_w1_5errs<-lapply(trs300_incolle_set1b, function(X)torus_disterror(X[["noizyX"]], maxr = 2.5, minr = 1, nps = 300))
par(mgp=c(2.4,1,0))
boxplot(in_trs300b_1_w1_5errs, xlab="Data Set", ylab="Error", cex.lab=1.6, cex.axis=1.6, lwd=2)

as.data.frame(cbind(1:31,inter_inv2)) %>%
  group_by(x) %>%
  filter(n()>1)


#全ボロノイ領域頂点へ補間
vic1_inter<-voronoi_vertex2(vics = trs300_vics1, figure = torus300)
points3d(vic1_inter, col=3)

in_torus300<-voronoi_gtm_interpo(figure = torus300, 30)
points3d(in_torus300, col=2)

#補間後に点を減らす
in_torus300B<-rbind(torus300, in_torus300)
dist_thre<-quantile_threshold(x = 0.8, X = in_torus300B)
in_trs300_cell<-cell_set2(x = in_torus300B, thresh = dist_thre)
in_trs300_cnct<-connect2(i = 1, cell_p = in_trs300_cell, all = 1:nrow(in_torus300B))
in_trs300_red<-reduce_inter_points(x = in_torus300B, conect = in_trs300_cnct, n_ori = 300)

in_trs300_red2<-reduce_intered(intered_X = in_torus300B, ratio = 0.8, n_ori = 300)
figurePlot3d(in_trs300_red2[["y"]][in_trs300_red2[["remain"]] <= 300, ])
points3d(in_trs300_red2[["y"]][in_trs300_red2[["remain"]] > 300, ], col=2)

#全ボロノイ領域の頂点へ補間、点数削減を行う一括補間
trs300_incolle_set1c<-gtm_inter_reduce(torus300_colle_set[[1]][1:3], 30, ratio = 0.7)

figurePlot3d(torus300_colle_set[[1]][[3]][["noizyX"]])
in_trs300_in1_3<-voronoi_gtm_interpo(figure = torus300_colle_set[[1]][[3]][["noizyX"]], 30)
points3d(in_trs300_in1_3, col=2)
