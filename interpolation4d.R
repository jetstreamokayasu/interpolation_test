#4次元データのに対してボロノイ補間
#3次元トーラスt3orus4_list3使用

t3orus4_list3_65<-t3orus4_list3[[65]][["noizyX"]]
t3orus4_list3_65_inst<-TDAdataset$new(t3orus4_list3_65)

#t3orus4_list3_65の1点目の近傍20点のインデックスを取得
t3orus4_list3_65_p1vics<-get_vicinity(dist = t3orus4_list3_65_inst$distmat, center = 1, nvic = 20)

figurePlot3d(t3orus4_list3_65_inst$data[-t3orus4_list3_65_p1vics, ])
points3d(t3orus4_list3_65_inst$data[t3orus4_list3_65_p1vics, ], col = 3)

vics_pca<-prcomp(t3orus4_list3_65_inst$data[t3orus4_list3_65_p1vics, ])
plot3d(vics_pca$x[, 1:3])
spheres3d(vics_pca$x[1, 1:3], col=2, radius = 0.03)

dela_tetra_list<-delaunay_triangulation3d(data = vics_pca$x[, 1:3])

for (t in 1:length(dela_tetra_list)) {
  dela_tetra_list[[t]]$plot_tetra_edges(edge_col = rainbow(length(dela_tetra_list))[t])
}

p1_dela_idx<-sapply(dela_tetra_list, function(tetra){tetra$is_vertex(vics_pca$x[1, 1:3])}) %>% which()

for (t in p1_dela_idx) {
  dela_tetra_list[[t]]$plot_tetra_edges(edge_col = rainbow(length(dela_tetra_list))[t])
}

p1_voro_verts<-sapply(dela_tetra_list, function(tetra){tetra$is_vertex(vics_pca$x[1, 1:3])}) %>% which() %>% 
  sapply(., function(idx)dela_tetra_list[[idx]]$cntr) %>% t()

spheres3d(p1_voro_verts, col = 4, radius = 0.03)

eigen01<-as.matrix(vics_pca$rotation[, 1])
eigen02<-as.matrix(vics_pca$rotation[, 2])
eigen03<-as.matrix(vics_pca$rotation[, 3])

p1_inted<-(p1_voro_verts %*% t(vics_pca$rotation[, 1:3])) %>% 
  apply(., 1, function(p)p + t3orus4_list3_65_inst$data[t3orus4_list3_65_p1vics[1], ]) %>% t()

points3d(p1_inted, col = 2)

vics_pca_max<-apply(vics_pca$x, 2, max)[1:3]
vics_pca_min<-apply(vics_pca$x, 2, min)[1:3]

p1_voro_idxs<-apply(p1_voro_verts, 1, function(p){
  
  all(map2_lgl(p, vics_pca_max, ~{.x <= .y}), map2_lgl(p, vics_pca_min, ~{.x >= .y}))
  
}) %>% which()

#------------------------
#関数テスト----------
p1_inted2<-voronoi_border4d(figure = t3orus4_list3_65_inst$data, vics = t3orus4_list3_65_p1vics)
points3d(p1_inted2[[2]], col = 4)

t3orus4_list3_65_inted<-voronoi_interpo4d(figure = t3orus4_list3_65_inst$data, n_vics = 20) %>% rbind(t3orus4_list3_65, .)

t3orus4_list3_65_inst$create_subsample(sub_size = t3orus4_list3_65_inst$n_points*0.8, n_subs = 1)


{inter_start<-Sys.time()
t3rs4_lst3_65_sub_inted<-t3orus4_list3_65_inst$subsamples[[1]]$data %>% 
  voronoi_interpo4d(figure = ., n_vics = 20) %>% rbind(t3orus4_list3_65_inst$subsamples[[1]]$data, .)
inter_end<-Sys.time()}

t3rs4_lst3_65_sub_inted_inst<-TDAdataset$new(t3rs4_lst3_65_sub_inted)
t3rs4_lst3_65_sub_inted_inst$calc_pd(maxdim = 3, maxscale = 9)
calc.landscape.peak(X = t3rs4_lst3_65_sub_inted_inst$get_pl()[[2]], dimension = 2, 
                    thresh = t3rs4_lst3_65_sub_inted_inst$get_pl()[["thresh"]]*(2*pi)/usephacm:::surface_nshpere(2))

t3rs4_lst3_65_sub_inted_inst$peaks<-
  map_dbl(seq_len(3), ~{calc.landscape.peak(X = t3rs4_lst3_65_sub_inted_inst$get_pl()[[.]], dimension = ., 
                                            thresh = t3rs4_lst3_65_sub_inted_inst$get_pl()[["thresh"]]*(2*pi)/usephacm:::surface_nshpere(.))
  })




calc.landscape.peak(X = t3rs4_lst3_65_sub_inted_inst$get_pl()[["2-land"]], dimension = 2, 
                    thresh = t3rs4_lst3_65_sub_inted_inst$get_pl()[["thresh"]]/2, tseq = t3rs4_lst3_65_sub_inted_inst$get_pl()[["tseq"]], show = T)

calc.landscape.peak(X = t3rs4_lst3_65_sub_inted_inst$get_pl()[["1-land"]], dimension = 1, 
                    thresh = t3rs4_lst3_65_sub_inted_inst$get_pl()[["thresh"]], tseq = t3rs4_lst3_65_sub_inted_inst$get_pl()[["tseq"]], show = T)

#-------------------------
#t3orus450_list2を用いて関数テスト----

t3rs450_lst2_77_inst<-TDAdataset$new(t3orus450_list2[[77]])
t3rs450_lst2_77_inst$calc_pd(maxdim = 3, maxscale = 9)

t3rs450_lst2_77_inst$plot_data()
points3d(t3rs450_lst2_77_inted_inst$data[451:t3rs450_lst2_77_inted_inst$n_points, ], col = 4)
movie3d(spin3d(axis = c(0, 1, 0), rpm = 5), duration = 30, 
        movie = "3d_torus_inted2", dir = "./pics", dev = 1, fps = 10)

{
  t3rs450_lst2_77_start_time<-Sys.time()
  
  t3rs450_lst2_77_inted<-t3rs450_lst2_77_inst$data %>% 
    voronoi_interpo4d(figure = ., n_vics = 30) %>% rbind(t3rs450_lst2_77_inst$data, .)
  
  t3rs450_lst2_77_end_time<-Sys.time()
}

t3rs450_lst2_77_inted_inst<-TDAdataset$new(t3rs450_lst2_77_inted)
t3rs450_lst2_77_inted_inst$calc_pd(maxdim = 3, maxscale = 9)

t3rs450_lst2_77_inted_inst$peaks<-
  map_dbl(seq_len(3), ~{calc.landscape.peak(X = t3rs450_lst2_77_inted_inst$get_pl()[[.]], dimension = ., 
                                            thresh = t3rs450_lst2_77_inted_inst$get_pl()[["thresh"]]*(2*pi)/usephacm:::surface_nshpere(.))
  })

#------------------------------
#補間点を減らしてみる-----

t3rs4_lst3_65_sub<-t3orus4_list3[[1]][["noizyX"]][sample(1:500, 400), ]

t3rs4_lst3_65_sub_inst<-TDAdataset$new(t3rs4_lst3_65_sub)

{inter_start<-Sys.time()
  t3rs4_lst3_65_sub_inted<-t3rs4_lst3_65_sub %>% 
    voronoi_interpo4d(figure = ., n_vics = 30) %>% rbind(t3rs4_lst3_65_sub, .)
  inter_end<-Sys.time()}

t3rs4_lst3_65_sub_inted_inst2<-TDAdataset$new(t3rs4_lst3_65_sub_inted)

t3rs4_lst3_65_sub_inted_inst<-TDAdataset$new(t3rs4_lst3_65_sub_inted[c(1:400, sample(401:nrow(t3rs4_lst3_65_sub_inted), (nrow(t3rs4_lst3_65_sub_inted)-400)*0.45) ), ])
t3rs4_lst3_65_sub_inted_inst$calc_pd(maxdim = 3, maxscale = 9)

#データ点数削減関数を使ってみる
{red_start<-Sys.time()
t3rs4_lst3_65_sub_red<-reduce_intered(intered_X = t3rs4_lst3_65_sub_inted, ratio = 0.75, n_ori = 400)
red_end<-Sys.time()
}

t3rs4_lst3_65_sub_red_inst<-TDAdataset$new(t3rs4_lst3_65_sub_red[["y"]])
t3rs4_lst3_65_sub_red_inst$calc_pd(maxdim = 3, maxscale = 9)

{red_start<-Sys.time()
  t3rs4_lst3_65_sub_red2<-reduce_intered(intered_X = t3rs4_lst3_65_sub_inted, ratio = 0.5, n_ori = 400)
  red_end<-Sys.time()
}

t3rs4_lst3_65_sub_red2_inst<-TDAdataset$new(t3rs4_lst3_65_sub_red2[["y"]])
