#補間を一括して行う関数
#補間点数を変更できるように
#中心点から遠い順に補う
variable_interpo<-function(collect, nvics, a=0){
  incollect<-collect
  for (l in 1:length(collect)) {
    inter_oricord<-voronoiInterpo2(collect[[l]][["noizyX"]], nvics, a)
    incollect[[l]][["noizyX"]]<-rbind(incollect[[l]][["noizyX"]], inter_oricord)
    incollect[[l]][["nsample"]]<-nrow(incollect[[l]][["noizyX"]])
    debugText(l, incollect[[l]][["nsample"]])
  }
  
  return(incollect)
}


#ボロノイ領域の頂点に補間を一括処理
#補間点数を変更できるように
#中心点から遠い順に補う
voronoiInterpo2<-function(figure, nvics, a){
  
  element<-rep(0, length = nrow(figure))
  
  dist<-dist(figure)
  
  for (i in 1:nrow(figure)) {
    
    if(element[i]==0){
      
      vics<-interpo3d:::get_vicinity(dist, i, nvics)
      
      element[vics]<-element[vics]+1
      
      #vics.oricord<-voronoiProcess(vics.line, figure)
      vics.oricord<-voronoiBorder2(vics, figure, a)[[1]]
      
      if(i==1){oricord<-vics.oricord}
      else{oricord<-rbind(oricord, vics.oricord)}
      
    }
    
  }
  
  #debugText(element)
  
  return(oricord)
  
}

#ボロノイ領域の頂点に点を打つ
#PCAで写された点による凸包内に入っていない点は除く
#補間点数を変更できるように
#中心点から遠い順に補う
voronoiBorder2<-function(vics, figure, a){
  
  require(deldir)
  
  vics.pca<-prcomp(figure[vics,])
  
  res<-deldir(vics.pca$x[,1], vics.pca$x[,2])
  
  tiles<-tile.list(res)
  
  insecs<-cbind(tiles[[1]][["x"]], tiles[[1]][["y"]])
  
  exist<-interpo3d:::exist_convexhull_check(vics.pca, insecs)
  
  if(length(insecs[which(exist==T), ]) > 0){
    
    insecs2<-insecs[which(exist==T), ]
    
    if(!is.matrix(insecs2)){insecs2<-t(as.matrix(insecs2))}
    
    vics_dist<-rbind(vics.pca[["x"]][1,1:2], insecs2) %>% dist() %>% 
      as.matrix(.) 
    
    idx<-order(-vics_dist[-1, 1])
    
    #debugText(vics_dist[-1, 1], idx, nrow(insecs2))
    
    if((nrow(insecs2) >= a) && (a != 0)){
  
    vics.oricord<-interpo3d:::origin_coordinate(vics.pca, insecs2[idx[1:a], ], figure[vics[1],])

    #debugText(insecs[which(exist==T)[1:a], ])
    
    return(list(oricord=vics.oricord, pca.inter=insecs[which(exist==T), ]))
    
    }else{
      
      vics.oricord<-interpo3d:::origin_coordinate(vics.pca, insecs[which(exist==T), ], figure[vics[1],])
      
      return(list(oricord=vics.oricord, pca.inter=insecs[which(exist==T), ]))
      
    }
    
  }
  
}

#補間手法を改良
#ボロノイ領域の頂点にの中で、最も中心点より遠い頂点に補間
#変更。中心点からの平均以上の点のみに補う
#データのリスト全体を補間
all_interpolate2<-function(collect, nvic){
  incollect<-collect
  for (l in 1:length(collect)) {
    inter_oricord<-voronoi_interpo3(collect[[l]][[2]], nvic)
    incollect[[l]][[2]]<-rbind(incollect[[l]][[2]], inter_oricord)
    incollect[[l]][[1]]<-nrow(incollect[[l]][[2]])
    cat("dataset", l, "has", incollect[[l]][[1]], "points\n")
  }
  
  return(incollect)
}

#figureに近傍nvic個で補間
voronoi_interpo3<-function(figure, nvics){
  
  element<-rep(0, length = nrow(figure))
  
  dist<-dist(figure)
  
  for (i in 1:nrow(figure)) {
    
    if(element[i]==0){
      
      vics<-get_vicinity(dist, i, nvics)
      
      element[vics]<-element[vics]+1
      
      vics_oricord<-voronoi_border3(vics, figure)[[1]]
      
      if(i==1){oricord<-vics_oricord}
      else{oricord<-rbind(oricord, vics_oricord)}
      
    }
    
  }
  
  return(oricord)
  
}


#ボロノイ領域の頂点にの中で、最も中心点より遠い頂点に補間
#変更。中心点からの平均以上の点のみに補う
#補間点の元の座標系の座標を返す
voronoi_border3<-function(vics, figure){
  
  vics_pca<-stats::prcomp(figure[vics,])
  
  tiles<-deldir::deldir(vics_pca$x[,1], vics_pca$x[,2]) %>% deldir::tile.list(.)
  
  insecs<-cbind(tiles[[1]][["x"]], tiles[[1]][["y"]])
  
  exist<-interpo3d:::exist_convexhull_check(vics_pca, insecs)
  
  if(length(insecs[which(exist==T), ]) > 0){
    
    vics_dist<-rbind(vics_pca[["x"]][1,1:2], insecs[which(exist==T), ]) %>% dist() %>% 
           as.matrix(.) 
    
    idx<-which(vics_dist[2:nrow(vics_dist), 1] >= mean(vics_dist[2:nrow(vics_dist), 1]))
    
    vics_oricord<-interpo3d:::origin_coordinate(vics_pca, insecs[which(exist==T)[idx], ], figure[vics[1],])
    
    return(list(oricord=vics_oricord, pca_inter=insecs[which(exist==T), ]))
    
  }
  
}

#中心点のボロノイ領域に隣接するボロノイ領域を見つける
neibor_voronoi<-function(tiles){

  neib_tile<-lapply(tiles[2:length(tiles)], function(tile){
    
    equal<-tiles[[1]][["x"]] %in% tile[["x"]]

    if(any(equal)){
      
      x_eq<-t(sapply(tiles[[1]][["x"]][equal], function(x){x==tile[["x"]]}))
      y_eq<-sapply(1:length(tiles[[1]][["x"]][equal]), function(k){tile[["y"]][x_eq[k,]]==tiles[[1]][["y"]][equal][k]})

      if(class(y_eq)=="list"){y_eq<-unlist(y_eq)}
      
      if(any(y_eq)){return(T)}else{return(F)}
      
    }
    
    else{return(F)}
    
  })
  
  return(neib_tile)
  
}


#GTMを用いた補間
#ボロノイ領域の頂点へ補間。データセットを与えて一括補間
#voronoi_vertexをvoronoi_vertex2へ。ボロノイ領域の全頂点へ補間
voronoi_gtm_interpo<-function(figure, nvics){
  
  element<-rep(0, length = nrow(figure))
  
  dist<-dist(figure)
  
  for (i in 1:nrow(figure)) {
    
    if(element[i]==0){
      
      vics<-get_vicinity(dist, i, nvics)
      
      element[vics]<-element[vics]+1
      
      vics_oricord<-voronoi_vertex2(vics, figure)
      
      #ボロノイ図がうまく作成できなかった場合
      #一度選ばれた点のインデックスを収納するelementを戻す
      #そして次のステップへ
      if(is.null(vics_oricord)){
        
        element[vics]<-element[vics]-1
        
      }else{
      
      if(i==1){oricord<-vics_oricord}
      else{oricord<-rbind(oricord, vics_oricord)}
        
      }
      
    }
    
  }
  
  #debugText(element)
  
  return(oricord)
  
}

#GTMを用いた補間
#ボロノイ領域の頂点へ補間。
voronoi_vertex<-function(vics, figure, MapsizeRow=15, MapsizeColumn=15, RBFsizeRow=3, RBFsizeColumn=3, RBFVariance=1, Lambda = 0, NumOfTraining = 100){
  
  require(gtm)
  
  vics_pca<-stats::prcomp(figure[vics,])
  
  #GTMを用いた次元削減
  X<-figure[vics,]
  # (分散が0の変数を削除した後に) 1. オートスケーリング
  Var0Variable <- which(apply(X,2,var) == 0)
  if (length(Var0Variable) == 0) {
    #print("分散が0の変数はありません")
  } else {
    sprintf("分散が0の変数が %d つありました", length(Var0Variable))
    print( "変数:" )
    print( Var0Variable )
    print( "これらを削除します" )
    X <- X[,-Var0Variable]
  }
  X <- scale(X, center = TRUE, scale = TRUE)
  # 6. データ空間における分散の逆数βの初期値
  # 7. 重みWの初期値
  XGrid = gtm.rctg( MapsizeColumn, MapsizeRow)#写像先のグリッド
  RBFGrid = gtm.rctg( RBFsizeColumn, RBFsizeRow)#基底関数のグリッド
  RBFSetup = gtm.gbf( RBFGrid, RBFVariance^(1/2), XGrid)
  InitnalWBeta = gtm.pci.beta( X, XGrid, RBFSetup)
  #Beta = 0.01
  Beta = InitnalWBeta$beta
  # 9. GTMマップ作成
  GTMResults = gtm.trn( X, RBFSetup, InitnalWBeta$W, Lambda, NumOfTraining, Beta, quiet = T)
  # 10. 二次元のマップ上でサンプルの位置関係を確認
  GTMDist = gtm.dist(X, RBFSetup %*% GTMResults$W)
  GTMR = gtm.resp3(GTMDist, GTMResults$beta, ncol(X))$R
  GTMMean = t(GTMR) %*% XGrid
  
  #次元削減後のボロノイ分割
  res<-deldir(GTMMean[,1], GTMMean[,2])
  tiles<-tile.list(res)
  
  #中心点のボロノイ領域に隣接するボロノイ領域を見つける
  nei_vo<-neibor_voronoi(tiles)
  nei_tile_list<-which(unlist(nei_vo)==T)+1
  
  nei_inter<-lapply(tiles[c(1, nei_tile_list)], function(tile){cbind(tile[["x"]], tile[["y"]])})
  nei_interpo<-unique(do.call(rbind, nei_inter))
  
  RBF_inter = gtm.gbf(RBFGrid, RBFVariance^(1/2), nei_interpo)
  #inter_dist<-gtm.dist(X, RBF_inter %*% GTMResults$W)
  # inter_R = gtm.resp3(inter_dist, GTMResults$beta, ncol(X))$R
  # inter_mean = t(inter_R) %*% (RBF_inter %*% GTMResults$W)
  inter_mean<-RBF_inter %*% GTMResults$W
  inter_inv<-sapply(1:nrow(inter_mean), function(k){inter_mean[k, ] * attr(X, "scaled:scale") + attr(X, "scaled:center")})

  return(t(inter_inv))
  
}

#GTMを用いた補間
#データセットのリストを一括補間
gtm_interpolate<-function(collect, nvic){
  incollect<-collect
  for (l in 1:length(collect)) {
    inter_oricord<-voronoi_gtm_interpo(collect[[l]][[2]], nvic)
    inter_oricord<-inter_oricord[!is.na(inter_oricord[,1]), ]
    incollect[[l]][[2]]<-rbind(incollect[[l]][[2]], inter_oricord)
    incollect[[l]][[1]]<-nrow(incollect[[l]][[2]])
    cat("dataset", l, "has", incollect[[l]][[1]], "points\n")
  }
  
  return(incollect)
}


voronoi_vertex2<-function(vics, figure, MapsizeRow=15, MapsizeColumn=15, RBFsizeRow=3, RBFsizeColumn=3, RBFVariance=1, Lambda = 0, NumOfTraining = 100){
  
  require(gtm)
  
  #vics_pca<-stats::prcomp(figure[vics,])
  
  #GTMを用いた次元削減
  X<-figure[vics,]
  # (分散が0の変数を削除した後に) 1. オートスケーリング
  Var0Variable <- which(apply(X,2,var) == 0)
  if (length(Var0Variable) == 0) {
    #print("分散が0の変数はありません")
  } else {
    sprintf("分散が0の変数が %d つありました", length(Var0Variable))
    print( "変数:" )
    print( Var0Variable )
    print( "これらを削除します" )
    X <- X[,-Var0Variable]
  }
  X <- scale(X, center = TRUE, scale = TRUE)
  # 6. データ空間における分散の逆数βの初期値
  # 7. 重みWの初期値
  XGrid = gtm.rctg( MapsizeColumn, MapsizeRow)#写像先のグリッド
  RBFGrid = gtm.rctg( RBFsizeColumn, RBFsizeRow)#基底関数のグリッド
  RBFSetup = gtm.gbf( RBFGrid, RBFVariance^(1/2), XGrid)
  InitnalWBeta = gtm.pci.beta( X, XGrid, RBFSetup)
  #Beta = 0.01
  Beta = InitnalWBeta$beta
  # 9. GTMマップ作成
  GTMResults = gtm.trn( X, RBFSetup, InitnalWBeta$W, Lambda, NumOfTraining, Beta, quiet = T)
  # 10. 二次元のマップ上でサンプルの位置関係を確認
  GTMDist = gtm.dist(X, RBFSetup %*% GTMResults$W)
  GTMR = gtm.resp3(GTMDist, GTMResults$beta, ncol(X))$R
  GTMMean = t(GTMR) %*% XGrid
  
  debugText(vics)
  
  #次元削減後のボロノイ分割
  res<-try(deldir(GTMMean[,1], GTMMean[,2]))
  
  #ボロノイ図がうまく作成できない場合
  #終了してNULLを返す
  if(is.null(res) || class(res)=="try-error"){
    
    return(NULL)
    
  }
  
  
  #debugText(res[["dirsgs"]][["x1"]], res[["dirsgs"]][["bp1"]])
  
  #境界領域に接しないボロノイ領域の頂点を選ぶ
  vertx<-cbind(res[["dirsgs"]][["x1"]][!res[["dirsgs"]][["bp1"]]], res[["dirsgs"]][["y1"]][!res[["dirsgs"]][["bp1"]]])
  RBF_inter = gtm.gbf(RBFGrid, RBFVariance^(1/2), vertx)
  inter_mean<-RBF_inter %*% GTMResults$W
  inter_inv<-sapply(1:nrow(inter_mean), function(k){inter_mean[k, ] * attr(X, "scaled:scale") + attr(X, "scaled:center")})
  
  return(t(inter_inv))
  
}

#GTMを用いた補間
#データセットのリストを一括補間
#全ボロノイ領域の頂点へ補間
#補間後に点を減らす
gtm_inter_reduce<-function(collect, nvic, ratio){
  incollect<-collect
  for (l in 1:length(collect)) {
    inter_oricord<-voronoi_gtm_interpo(collect[[l]][[2]], nvic)
    inter_oricord<-inter_oricord[!is.na(inter_oricord[,1]), ]
    red_oricord<-reduce_intered(intered_X = rbind(collect[[l]][[2]], inter_oricord), ratio = ratio, n_ori = nrow(collect[[l]][[2]]))
    incollect[[l]][[2]]<-red_oricord[["y"]]
    incollect[[l]][[1]]<-nrow(incollect[[l]][[2]])
    cat("dataset", l, "has", incollect[[l]][[1]], "points\n")
  }
  
  return(incollect)
}
