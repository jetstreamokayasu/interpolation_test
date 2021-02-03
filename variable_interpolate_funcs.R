#補間点数可変プログラム群

#------------------------------
#近傍点、近傍点のPCA、PCA適用後のボロノイ分割をまとめたクラス-----
#とりあえず3次元データ用。ボロノイ分割は2次元平面上

VicsPcaVoronoi<-R6Class(
  
  classname = "VicsPcaVoronoi", 
  
  public = list(
    
    vics_idx = NA, #近傍点のfigure上のインデックス
    vics = NA, #近傍点の座標
    vics_pca = NA, #近傍点へのPCA適用結果
    res = NA, #ドロネー/ボロノイ分割結果
    tiles = NA, #ボロノイ領域の集合
    
    initialize = function(figure, center, nvics){
      
      self$vics_idx<-dist(figure) %>% as.matrix() %>%  
        interpo3d::get_vicinity(dist = ., center = center, nvic = nvics)
      
      self$vics<-figure[self$vics_idx, ]
      
      self$vics_pca<-prcomp(self$vics)
      
      self$res<-deldir::deldir(self$vics_pca$x[, 1], self$vics_pca$x[, 2])
      
      self$tiles<-deldir::tile.list(self$res)
      
    },
    
    #ボロノイ領域をpolygonで描画する関数
    plot_tiles = function(...){
      
      for(i in seq_len(self$res$n.data)){	polygon(self$tiles[[i]], ...) }
      
    }
    
  )#publicの閉じ括弧
  
)#VicsPcaVoronoiクラスの閉じ括弧

#----------------------------
#隣接するボロノイ領域の頂点にも補間する関数----
#現時点では全ボロノイ領域の頂点に補間する
#vicinityはVicsPcaVoronoiクラス
#neighbor=いくつの隣接するボロノイ領域の頂点を対象とするか。"all"とすると全ての頂点を対象とする

neighbor_voronoi_vertex<-function(vicinity, figure, neighbor = 0){
  
  if(!is.numeric(neighbor) && neighbor != "all"){stop("neighbor is incorrect.")}
  
  if(neighbor == "all"){
    
    vertxs<-vicinity$res[["dirsgs"]] %$%  cbind(x1[!bp1], y1[!bp1])
    
  }else{
    
    vertxs<-cbind(vicinity$tiles[[1]][["x"]], vicinity$tiles[[1]][["y"]])
    
    if(neighbor > 0){
      
      voro_idx<-vicinity$res[["dirsgs"]] %$% c(ind1[ind2 == 1], ind2[ind1 == 1]) %>% sort() %>% magrittr::extract(1:min(neighbor, length(.)))
      
      candi_vertxs_idx1<-vicinity$res[["dirsgs"]] %$% (is.element(ind1, voro_idx) & (ind2 != 1) & not(bp1))
      candi_vertxs_idx2<-vicinity$res[["dirsgs"]] %$% (is.element(ind2, voro_idx) & (ind1 != 1) & not(bp2))

      vertxs<-vicinity$res[["dirsgs"]] %$% rbind(cbind(x1[candi_vertxs_idx1], y1[candi_vertxs_idx1]), cbind(x2[candi_vertxs_idx2], y2[candi_vertxs_idx2]) ) %>% 
        rbind(vertxs, .)
      
      dupl_idx<-lapply(seq_len(nrow(vertxs)-1), function(i){sapply((i+1):nrow(vertxs), function(j){
        
        dupl<-c()
        
        if(identical(vertxs[i, ], vertxs[j, ])){dupl<-c(dupl, j)}
        
        if(length(dupl)==0){dupl<-0}
        
        return(dupl)
        
      })}) %>%  unlist() %>% unique() %>% magrittr::extract(.!=0)
      
      vertxs<-vertxs[-dupl_idx, ]
      
    }
    
  }
  
  vertx_oricord<-(vertxs %*% t(vicinity$vics_pca$rotation[, 1:2])) %>%
    apply(., 1, function(p){p + figure[vicinity$vics_idx[1], ]}) %>% t()

  return(lst(vertx_oricord = vertx_oricord, voro_vertx = vertxs))
  
}

#---------------------------------------
#隣接するボロノイ領域の頂点へも補間する処理を、データ全体に行う関数----
#現時点では全ボロノイ領域の頂点に補間する
#figure=補間対象データ、nvics=近傍点数

neighbor_voronoi_interpol<-function(figure, nvics, neighbor = 0){
  
  element<-rep(0, length = nrow(figure))
  
  inted_point<-c()
  
  for (i in seq_len(nrow(figure))) {
    
    if(element[i]==0){
      
      vicinity<-VicsPcaVoronoi$new(figure, i, nvics)
      
      element[vicinity$vics_idx]<-element[vicinity$vics_idx]+1
      
      inted_point<-neighbor_voronoi_vertex(vicinity, figure, neighbor)[[1]] %>% rbind(inted_point, .)
      
    }
    
  }
  
  return(inted_point)
  
}
