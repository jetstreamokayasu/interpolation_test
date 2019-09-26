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
    
    debugText(vics_dist[-1, 1], idx, nrow(insecs2))
    
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


