#start, endで区切られた範囲にデータ点が存在しするかしないかを判定する関数
#存在する場合はT、存在しない場合はFを返す
existCheck<-function(start, end, mapped){
  
  exie1.mem<-which(mapped[,1]>=start[1] & mapped[,1]<end[1])
  
  debugText(exie1.mem, start, end)
  
  if(length(exie1.mem)>=1){
    
    exie2.mem<-which(mapped[exie1.mem, 2]>=start[2] & mapped[exie1.mem ,2]<end[2])
    
    debugText(exie2.mem)
    
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
