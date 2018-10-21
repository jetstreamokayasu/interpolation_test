#start, endで区切られた範囲にデータ点が存在しするかしないかを判定する関数
#存在する場合はT、存在しない場合はFを返す
existCheck<-function(start, end, mapped){
  
  exie1.mem<-which(mapped[,1]>=start[1] & mapped[,1]<end[1])
  
  #debugText(exie1.mem, start, end)
  
  if(length(exie1.mem)>=1){
    
    exie2.mem<-which(mapped[exie1.mem, 2]>=start[2] & mapped[exie1.mem ,2]<end[2])
    
    #debugText(exie2.mem)
    
    if(length(exie2.mem)>=1) return(T)
    
    else return(F)
    
  }
  
  else return(F)
  
}

#格子状の線を引く
gridLine<-function(x, div){
  
  xlim<-range(x[,1])
  ylim<-range(x[,2])
  
  wid<-abs(xlim[2]-xlim[1])/div
  hei<-abs(ylim[2]-ylim[1])/div
  
  for (i in 0:div) {
    
    abline(v=xlim[1]+i*wid)
    abline(h=ylim[1]+i*hei)
    
  }
  
}


pixelConvert<-function(x, div){
  
  xlim<-range(x[,1])
  ylim<-range(x[,2])
  
  wid<-abs(xlim[2]-xlim[1])/div
  hei<-abs(ylim[2]-ylim[1])/div
  
  pixel<-matrix(0, div, div)
  
  for (i in 1:div) {
    
    for (j in 1:div) {
      
      if(existCheck(c(xlim[1]+(i-1)*wid, ylim[1]+(j-1)*hei), c(xlim[1]+i*wid, ylim[1]+j*hei), x)){
    
            pixel[j, i]<-1
        
      }
    }
    
  }
  
  return(pixel)
    
  }

#行列の注目した要素の8近傍に要素が入っているか
#入っていればTを返す
neighEleCheck<-function(pic, row, col){
  
  if(row!=1 && (col!=1) && (pic[row-1, col-1]>0)){return(T)}
  else if((row!=1) && pic[row-1, col]>0){return(T)}
  else if((row!=1) && col!=ncol(pic) && pic[row-1, col+1]>0){return(T)}
  else if(col!=1 && pic[row, col-1]>0){return(T)}
  else if(col!=ncol(pic) && pic[row, col+1]>0){return(T)}
  else if(row!=nrow(pic) && col!=1 && pic[row+1, col-1]>0){return(T)}
  else if(row!=nrow(pic) && pic[row+1, col]>0){return(T)}
  else if(row!=nrow(pic) && col!=ncol(pic) && pic[row+1, col+1]>0){return(T)}
  else{return(F)}
  
}  

#注目要素に要素が無く、注目要素の8近傍に要素があるとき
#注目要素に2を代入
insertElement<-function(pic){
  
  cp.pic<-pic
  
  for (i in 1:nrow(pic)) {
    
    for (j in 1:ncol(pic)) {
      
      if(pic[i,j]==0 &&neighEleCheck(pic, i, j)){
        
        cp.pic[i, j]<-2
        
      }
      
    }
    
  }
  
  return(cp.pic)
  
}

#指定された範囲にPCA後の座標を返す
pcaCoordinate<-function(xmin, ymin, wid, hei, row, col){
  
  return(c(xmin+(wid*(col-1))+wid/2, ymin+(hei*(row-1))+hei/2))
  
}

#データ点が存在しないピクセルの中央にデータ点を打つ
#座標はPCA後の座標
pcaCoord.set<-function(x, cppic, div){
  
  xlim<-range(x[,1])
  ylim<-range(x[,2])
  
  wid<-abs(xlim[2]-xlim[1])/div
  hei<-abs(ylim[2]-ylim[1])/div
  
  ele<-which(cppic==2, arr.ind=TRUE)
  
  coord<-sapply(1:nrow(ele), function(k)pcaCoordinate(xlim[1], ylim[1], wid, hei, ele[k,1], ele[k,2]))
  
  return(t(coord))
  
}

#PCA適用後の座標から元座標を算出
originCoodinate<-function(rpca, incord){
  
  eigen01<-as.matrix(rpca$rotation[,1])
  eigen02<-as.matrix(rpca$rotation[,2])
  
  oricord<-sapply(1:nrow(incord), function(l){
    
  return((incord[l, 1]*(eigen01)+incord[l, 2]*(eigen02))+rpca$center)
    
  })
  
  return(t(oricord))
  
}

#任意の点に近傍点の中から最も遠い点を抽出
#そのもっとも遠い点から近傍の中で最も近いn点を抽出
coveredVic<-function(vicsline, figure, n){
  
  distfarvic<-sapply(1:(length(vicsline)-2), function(k){
    
    dist.set<-c(vicsline[length(vicsline)], vicsline[k+1], sum((figure[vicsline[k+1],]-figure[vicsline[length(vicsline)],])^2))
    names(dist.set)<-c("start", "goal", "distance")
    
    return(dist.set)
    
  })
  
  distfarvic<-t(distfarvic)
  
  distfarvic<-distfarvic[order(distfarvic[,3]),]
  
  return(distfarvic[1:n,])
  
}

#指定されたデータ点のnvics点近傍をPCAで変換し
#膨張処理を行い、補間されたデータ点の下\元の座標系での座標を返す
expandProcess<-function(vics, vics.line, figure, dist, div){
  
  vics.pca<-prcomp(figure[vics.line,])
  
  vics.pic<-pixelConvert(vics.pca[["x"]], div)
  
  vics.cppic<-insertElement(vics.pic)
  
  vics.incord<-pcaCoord.set(vics.pca[["x"]], vics.cppic, div)
  
  vics.oricord<-originCoodinate(vics.pca, vics.incord)
  
  return(vics.oricord)
  
}

interPolation_test<-function(figure, nvics, div){
  
  element<-rep(0, length = nrow(figure))
  
  dist<-distance(figure)
  
  for (i in 1:nrow(figure)) {
    
    if(element[i]==0){
      
      vics<-get.vicinity(dist, i, nvics)
      
      vics.line<-line.vics(i, vics)
      
      element[vics.line]<-element[vics.line]+1
      
      vics.oricord<-expandProcess(vics, vics.line, figure, dist, div)
      
      if(i==1){oricord<-vics.oricord}
      else{oricord<-rbind(oricord, vics.oricord)}
      
    }
    
  }
  
  debugText(element)
  
  return(oricord)
  
}
