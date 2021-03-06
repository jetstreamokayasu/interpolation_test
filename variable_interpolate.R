#補間点数を変えてみる

#-------------------------
#ボロノイ分割図の再確認-----

trs300_cole1_43<-TDAdataset$new(torus300_colle_set[[1]][[43]][["noizyX"]])

trs300_cole1_43_201vics<-get_vicinity(dist = trs300_cole1_43$distmat, center = 201, nvic = 10)

trs300_cole1_43_201pca<-prcomp(trs300_cole1_43$data[trs300_cole1_43_201vics, ])
plot(trs300_cole1_43_201pca$x[, 1:2])
points(trs300_cole1_43_201pca$x[1,1], trs300_cole1_43_201pca$x[1,2], pch = 16, col = 2)

trs300_cole1_43_201res<-deldir(trs300_cole1_43_201pca$x[,1], trs300_cole1_43_201pca$x[,2])
trs300_cole1_43_201tiles<-tile.list(trs300_cole1_43_201res)

for(i in 1:trs300_cole1_43_201res$n.data){	polygon(trs300_cole1_43_201tiles[[i]]) }

#VicsPcaVoronoiクラス試し
trs300_cole1_43_201voro<-VicsPcaVoronoi$new(figure = trs300_cole1_43$data, center = 201, nvics = 10)

#ボロノイ領域プロット
plot(trs300_cole1_43_201voro$vics_pca$x[, 1:2])
points(trs300_cole1_43_201voro$vics_pca$x[1,1], trs300_cole1_43_201voro$vics_pca$x[1,2], pch = 16, col = 2)

trs300_cole1_43_201voro$plot_tiles()

points(trs300_cole1_43_201res$dirsgs$x1[14], trs300_cole1_43_201res$dirsgs$y1[14], pch = 16, col = "blue")
points(trs300_cole1_43_201pca$x[3,1], trs300_cole1_43_201pca$x[3,2], pch = 16, col = "green")

#PCA後の点2,6,8を含むボロノイ領域の頂点で、点1のボロノイ領域の頂点ではないものを選出
tst_idx<-(is.element(trs300_cole1_43_201res$dirsgs[["ind1"]], c(2, 6, 8)) | is.element(trs300_cole1_43_201res$dirsgs[["ind2"]], c(2, 6, 8))) %>%  
  `&`(., (trs300_cole1_43_201res$dirsgs[["ind1"]] != 1)) %>% `&`((trs300_cole1_43_201res$dirsgs[["ind2"]] != 1)) %>% which()

tst_idx2<-trs300_cole1_43_201res[["dirsgs"]] %$% ((ind1 %in% c(2, 6, 8)) | (ind2 %in% c(2, 6, 8))) %>% 
  and(trs300_cole1_43_201res[["dirsgs"]]$ind1 != 1) %>% and(trs300_cole1_43_201res[["dirsgs"]]$ind2 != 1)

points(trs300_cole1_43_201res$dirsgs$x1[!trs300_cole1_43_201res$dirsgs$bp1], trs300_cole1_43_201res$dirsgs$y1[!trs300_cole1_43_201res$dirsgs$bp1], pch = 16, col = "blueviolet")

#------------------------------
#重複するボロノイ領域の頂点を消す-------
voro_vertx<-cbind(trs300_cole1_43_201res$dirsgs$x1, trs300_cole1_43_201res$dirsgs$y1)

dupl_tst<-lapply(seq_len(nrow(voro_vertx)-1), function(i){sapply((i+1):nrow(voro_vertx), function(j){
  
  dupl<-c()
  
  if(identical(voro_vertx[i, ], voro_vertx[j, ])){dupl<-c(dupl, j)}
  
  if(length(dupl)==0){dupl<-0}
  
  return(dupl)
  
})})
no_dupl<-unlist(dupl_tst) %>% unique() %>% magrittr::extract(.!=0)


#--------------------------------------------
#隣接するボロノイ領域の頂点への補間を試す--------

#まずは全ボロノイ領域の頂点へ補う

trs300_cole1_43_201inted<-neighbor_voronoi_vertex(vicinity = trs300_cole1_43_201voro, figure = trs300_cole1_43$data, neighbor = "all")

points(trs300_cole1_43_201inted$voro_vertx[, 1], trs300_cole1_43_201inted$voro_vertx[, 2], pch = 16, col = "orange")

figurePlot3d(trs300_cole1_43$data[-trs300_cole1_43_201voro$vics_idx, ])
points3d(trs300_cole1_43$data[trs300_cole1_43_201voro$vics_idx, ], col = "green")
points3d(trs300_cole1_43_201inted$vertx_oricord, col = "red")

#データ全体の補間
trs300_cole1_43_inted<-neighbor_voronoi_interpol(figure = trs300_cole1_43$data, nvics = 10)
trs300_cole1_43$plot_data()
points3d(trs300_cole1_43_inted, col = "red")


#隣接するボロノイ領域のいくつ対象とするか指定

#指定数0(中央のボロノイ領域頂点にのみ補う)
trs300_cole1_43_201_toadd2<-neighbor_voronoi_vertex(vicinity = trs300_cole1_43_201voro, figure = trs300_cole1_43$data)
points(trs300_cole1_43_201_toadd2$voro_vertx[, 1], trs300_cole1_43_201_toadd2$voro_vertx[, 2], pch = 16, col ="blueviolet")

#指定数1(中央および隣接する1つのボロノイ領域頂点にのみ補う)
trs300_cole1_43_201_toadd_nei1<-neighbor_voronoi_vertex(vicinity = trs300_cole1_43_201voro, figure = trs300_cole1_43$data, neighbor = 1)
points(trs300_cole1_43_201_toadd_nei1$voro_vertx[, 1], trs300_cole1_43_201_toadd_nei1$voro_vertx[, 2], pch = 16, col ="blueviolet")

#指定数2(中央および隣接する2つのボロノイ領域頂点にのみ補う)
trs300_cole1_43_201_toadd_nei2<-neighbor_voronoi_vertex(vicinity = trs300_cole1_43_201voro, figure = trs300_cole1_43$data, neighbor = 2)
points(trs300_cole1_43_201_toadd_nei2$voro_vertx[, 1], trs300_cole1_43_201_toadd_nei2$voro_vertx[, 2], pch = 16, col ="orange")

points3d(trs300_cole1_43_201_toadd_nei2$vertx_oricord, col = "red")

#指定数3(中央および隣接する3つのボロノイ領域頂点にのみ補う)
trs300_cole1_43_201_toadd_nei3<-neighbor_voronoi_vertex(vicinity = trs300_cole1_43_201voro, figure = trs300_cole1_43$data, neighbor = 3)
points(trs300_cole1_43_201_toadd_nei3$voro_vertx[, 1], trs300_cole1_43_201_toadd_nei3$voro_vertx[, 2], pch = 16, col ="blueviolet")

#全ボロノイ領域の頂点へ補う
trs300_cole1_43_201_toadd_allnei<-neighbor_voronoi_vertex(vicinity = trs300_cole1_43_201voro, figure = trs300_cole1_43$data, neighbor = "all")
points(trs300_cole1_43_201_toadd_allnei$voro_vertx[, 1], trs300_cole1_43_201_toadd_allnei$voro_vertx[, 2], pch = 16, col ="blueviolet")

points3d(trs300_cole1_43_201_toadd_allnei$vertx_oricord, col = "red")

#指定数3(中央および隣接する3つのボロノイ領域頂点にのみ補う)でデータ全体を補間
trs300_cole1_43_toadd_nei3<-neighbor_voronoi_interpol(figure = trs300_cole1_43$data, nvics = 10, neighbor = 3)
trs300_cole1_43$plot_data()
points3d(trs300_cole1_43_toadd_nei3, col = "red")

#指定数1(中央および隣接する1つのボロノイ領域頂点にのみ補う)でデータ全体を補間
trs300_cole1_43_toadd_nei1<-neighbor_voronoi_interpol(figure = trs300_cole1_43$data, nvics = 10, neighbor = 1)
trs300_cole1_43$plot_data()
points3d(trs300_cole1_43_toadd_nei1, col = "red")

#指定数6(中央および隣接する6つのボロノイ領域頂点にのみ補う)でデータ全体を補間
#隣接数を過大に設定するとどうなるか確かめる
trs300_cole1_43_toadd_nei6<-neighbor_voronoi_interpol(figure = trs300_cole1_43$data, nvics = 10, neighbor = 6)
trs300_cole1_43$plot_data()
points3d(trs300_cole1_43_toadd_nei6, col = "red")

#全ボロノイ領域の頂点への補間をデータ全体に
trs300_cole1_43_toadd_allnei<-neighbor_voronoi_interpol(figure = trs300_cole1_43$data, nvics = 10, neighbor = "all")
trs300_cole1_43$plot_data()
points3d(trs300_cole1_43_toadd_allnei, col = "red")
