#補間を一括して行う関数
#補間点数を変更できるように
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
voronoiBorder2<-function(vics, figure, a){
  
  require(deldir)
  
  vics.pca<-prcomp(figure[vics,])
  
  res<-deldir(vics.pca$x[,1], vics.pca$x[,2])
  
  tiles<-tile.list(res)
  
  insecs<-cbind(tiles[[1]][["x"]], tiles[[1]][["y"]])
  
  exist<-exist_convexhull_check(vics.pca, insecs)
  
  if(length(insecs[which(exist==T), ]) > 0){
    
    insecs2<-insecs[which(exist==T), ]
    
    if(!is.matrix(insecs2)){insecs2<-t(as.matrix(insecs2))}
    
    #debugText(nrow(insecs2), a)
    
    if((nrow(insecs2) >= a) && (a != 0)){
  
    #vics.oricord<-interpo3d:::origin_coordinate(vics.pca, insecs[which(exist==T)[1:a], ], figure[vics[1],])
    vics.oricord<-originCoodinate(vics.pca, insecs[which(exist==T)[1:a], ])

    #debugText(insecs[which(exist==T)[1:a], ])
    
    return(list(oricord=vics.oricord, pca.inter=insecs[which(exist==T), ]))
    
    }else{
      
      vics.oricord<-interpo3d:::origin_coordinate(vics.pca, insecs[which(exist==T), ], figure[vics[1],])
      
      return(list(oricord=vics.oricord, pca.inter=insecs[which(exist==T), ]))
      
    }
    
  }
  
}
