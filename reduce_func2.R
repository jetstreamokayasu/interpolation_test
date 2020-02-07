#補間後のデータ点を減らす関数群
#元のデータ点は残すようにする

reduce_inter_points<-function(x, conect, n_ori){
  reduc<-0
  for (i in 1:length(conect)) {
    if(length(conect[[i]])!=1){
      debugText(conect[[i]])
      mean_cnct<-apply(x[conect[[i]], ], 2, mean)
      cnct_dist<-rbind(mean_cnct, x[conect[[i]], ]) %>% 
        dist()  %>% 
        as.matrix() %>% 
        .[1, -1] %>% 
        as.vector()
      
      if(sum(conect[[i]] <= n_ori) >= 1){
        
        remain<-conect[[i]][conect[[i]] <= n_ori][which.min(cnct_dist[conect[[i]] <= n_ori])]
      
      }else{
        remain<-conect[[i]][which.min(cnct_dist)]
        }
      
      debugText(setdiff(conect[[i]], remain), remain)
      if(length(reduc)==1){reduc<-setdiff(conect[[i]], remain)}
      else{reduc<-c(reduc, setdiff(conect[[i]], remain))}
      
    }
  }
  y<-x[-reduc, ]
  
  return(list(y=y, remain=setdiff(1:nrow(x), reduc)))
  
}

#点数削減関数まとめ
## X=減らす前のデータ, ratio=距離の閾値の下位から数える割合, n_ori=補間前のデータセットの点数

reduce_intered<-function(intered_X, ratio, n_ori){
  
  dist_th<-quantile_threshold(x = ratio, X = intered_X)
  cell<-cell_set2(x = intered_X, thresh = dist_th)
  cnct<-connect2(i = 1, cell_p = cell, all = 1:nrow(intered_X))
  red<-reduce_inter_points(x = intered_X, conect = cnct, n_ori = n_ori)
  
  return(red)
  
}


#点数削減関数まとめ
#補間されたデータのみを減らす
## X=減らす前のデータ, ratio=距離の閾値の下位から数える割合

reduce_intered2<-function(intered_X, ratio, n_ori){
  
  dist_th<-quantile_threshold(x = ratio, X = intered_X)
  cell<-cell_set2(x = intered_X, thresh = dist_th)
  cnct<-connect2(i = 1, cell_p = cell, all = 1:nrow(intered_X))
  red<-reduce_points(x = intered_X, conect = cnct)
  
  return(red)
  
}
