#ボロノイ分割に関する関数群

#--------------------
#四面体(Tetrahedron)クラス

Tetrahedron<-R6Class(
  classname = "Tetrahedron",
  
  public = list(
    
    #四面体の頂点
    p1 = NA,
    p2 = NA,
    p3 = NA,
    p4 = NA,
    
    #頂点のリスト
    vertexes = lst(),
    
    #外接円の中心
    cntr = NA,
    
    #外接円の半径
    radius = NA,
    
    initialize = function(p1, p2, p3, p4){
      
      self$p1<-p1
      self$p2<-p2
      self$p3<-p3
      self$p4<-p4
      
      self$vertexes<-lst(self$p1, self$p2, self$p3, self$p4)
      
      self$calc_center()
      self$calc_radius()
      
    },
    
    calc_center = function(){
      
      self$cntr<-(1/2) * solve( rbind((self$p1-self$p2), (self$p1-self$p3), (self$p3-self$p4)) ) %*% 
        c((self$p1%*%self$p1 - self$p2%*%self$p2), (self$p1%*%self$p1 - self$p3%*%self$p3), (self$p3%*%self$p3 - self$p4%*%self$p4)) %>% as.vector()
      
    },
    
    calc_radius = function(){
      
      self$radius<-(self$p1 - self$cntr) %*% (self$p1 - self$cntr) %>% sqrt()
      
    },
    
    sphere_include_point = function(check_p){
      
      distance<-(check_p - self$cntr) %*% (check_p - self$cntr) %>% sqrt()
      
      if(distance < self$radius){return(TRUE)}else{return(FALSE)}
      
    },
    
    is_equal_tetra = function(comp_tetra){
      
      if(comp_tetra$radius == self$radius && all(comp_tetra$cntr == self$cntr)){return(TRUE)}
      else{return(FALSE)}
      
    },
    
    plot_tetra_edges = function(edge_col = 1, ...){
      
      if(rgl.cur() != 0){
        
        points3d(rbind(self$p1, self$p2, self$p3, self$p4), ...)
        lines3d(rbind(self$p1, self$p2), col = edge_col)
        lines3d(rbind(self$p1, self$p3), col = edge_col)
        lines3d(rbind(self$p1, self$p4), col = edge_col)
        lines3d(rbind(self$p2, self$p3), col = edge_col)
        lines3d(rbind(self$p2, self$p4), col = edge_col)
        lines3d(rbind(self$p3, self$p4), col = edge_col)
        
      }
      
    },
    
    #check_pがこの四面体の頂点か否か
    is_vertex = function(check_p){
      
      return(some(self$vertexes, ~{is.element(check_p, .) %>% all()}))
      
    }
    
  )#publicの閉じ括弧
  
)

#--------------------------
#Tetrahedronを与えて、四面体の辺をrglで描画する関数-----
#lines3dに任意変数を与えるための関数
plot_tetrahedron<-function(tetra, ...){
  
  lines3d(rbind(tetra$p1, tetra$p2), ...)
  lines3d(rbind(tetra$p1, tetra$p3), ...)
  lines3d(rbind(tetra$p1, tetra$p4), ...)
  lines3d(rbind(tetra$p2, tetra$p3), ...)
  lines3d(rbind(tetra$p2, tetra$p4), ...)
  lines3d(rbind(tetra$p3, tetra$p4), ...)
  
}

#------------------------
#３次元データに対してドロネー四面体分割を行う関数
#data=3次元データ、leave_first=最初の四面体を残すか否か
delaunay_triangulation3d<-function(data, leave_first = F){
  
  #初期四面体を作る
  
  max_x<-max(data[, 1])*1.1
  max_y<-max(data[, 2])*1.1
  max_z<-max(data[, 3])*1.1
  
  min_x<-min(data[, 1])*1.1
  min_y<-min(data[, 2])*1.1
  min_z<-min(data[, 3])*1.1
  
  #データ点をすべて含む立方体の頂点
  Q<-c(min_x, min_y, min_z)
  F_m<-c(max_x, max_y, max_z)
  C<-c(max_x, min_y, min_z)
  
  #データ点をすべて含む立方体の中心
  G<-(Q+F_m)/2
  
  #データ点をすべて含む立方体の外接球の半径
  #ドロネー分割の初期四面体の内接球とする
  r_i<-sqrt((F_m-G) %*% (F_m-G)) %>% as.vector()
  
  #初期四面体の底面正三角形の重心
  G2<-G
  G2[3]<-G[3]-r_i
  
  #初期四面体の底面正三角形における一辺の中心点
  H<-G2
  H[2]<-G2[2]-sqrt(2)*r_i
  
  #データ点をすべて含む立方体の辺QC方向の単位ベクトル
  #要するにx軸単位ベクトル
  e_QC<-c(1, 0, 0)
  
  #初期四面体の底面正三角形の頂点
  S<-H+as.vector(sqrt(6)*r_i*e_QC)
  T_v<-H-sqrt(6)*r_i*e_QC
  U<-G2
  U[2]<-G2[2]+2*sqrt(2)*r_i
  
  #初期四面体の底面正三角形に相対する頂点
  V<-G2
  V[3]<-G2[3]+4*r_i
  
  #初期四面体S-(T_v)-U-V
  first_tetra<-Tetrahedron$new(S, T_v, U, V)
  
  #ドロネー分割を実現する四面体のリスト
  div_tetra_list<-lst(first_tetra)
  
  #ドロネー分割
  for (i in 1:nrow(data)) {
    
    #追加するデータ点の選択
    select_p<-data[i, ]
    
    # 四面体の外接球内に追加したポイントが内包されていないか判定。
    temp_div_tetra<-lst()
    
    cnt<-0
    
    for(j in 1:length(div_tetra_list)){
      
      tri<-div_tetra_list[[j+cnt]]
      
      check<-tri$sphere_include_point(select_p)
      
      if(check){
        
        # 内包していた四面体の頂点を使用して再分割する。
        new_tetra1<-Tetrahedron$new(tri$p1, tri$p2, tri$p3, select_p)
        new_tetra2<-Tetrahedron$new(tri$p1, tri$p2, tri$p4, select_p)
        new_tetra3<-Tetrahedron$new(tri$p1, tri$p3, tri$p4, select_p)
        new_tetra4<-Tetrahedron$new(tri$p2, tri$p3, tri$p4, select_p)
        
        temp_div_tetra<-append(temp_div_tetra, lst(new_tetra1, new_tetra2, new_tetra3, new_tetra4))
        
        div_tetra_list<-div_tetra_list[-(j+cnt)]
        cnt<-cnt-1
        
      }
      
    }
    
    #重複している四面体を削除
    del_tetras<-lst()
    
    for (k in 1:length(temp_div_tetra)) {
      for (l in 1:length(temp_div_tetra)) {
        
        if(k == l){next}
        
        if(temp_div_tetra[[k]]$radius == temp_div_tetra[[l]]$radius && 
           all(temp_div_tetra[[k]]$cntr == temp_div_tetra[[l]]$cntr)){
          
          #削除候補四面体をリストに追加
          del_tetras<-append(del_tetras, lst(temp_div_tetra[[k]], temp_div_tetra[[l]]))
          
        }
        
      }
    }
    
    for (d_tet in del_tetras) {
      
      eq_tetra<-sapply(temp_div_tetra, function(tetra)d_tet$is_equal_tetra(tetra))
      
      if(any(eq_tetra)){
        
        eq_idx<-which(eq_tetra == TRUE)
        temp_div_tetra<-temp_div_tetra[-eq_idx]
        
      }
      
    }
    
    # 重複四面体を削除した後にappend
    div_tetra_list<-append(div_tetra_list, temp_div_tetra)
    
  }
  
  if(!leave_first){
    
    # はじめに作成した四面体の母点をもつ四面体分割を削除------
    del_tetras<-lst()
    
    for (tetra_check in div_tetra_list) {
      
      first_points<-lst(first_tetra$p1, first_tetra$p2, first_tetra$p3, first_tetra$p4)
      
      if(some(first_points, ~{all(. == tetra_check$p1)}) || some(first_points, ~{all(. == tetra_check$p2)}) ||
         some(first_points, ~{all(. == tetra_check$p3)}) || some(first_points, ~{all(. == tetra_check$p4)})){
        
        del_tetras<-append(del_tetras, lst(tetra_check))
        
      }
      
    }
    
    for (d_tet in del_tetras) {
      
      eq_tetra<-sapply(div_tetra_list, function(tetra)d_tet$is_equal_tetra(tetra))
      
      if(any(eq_tetra)){
        
        eq_idx<-which(eq_tetra == TRUE)
        div_tetra_list<-div_tetra_list[-eq_idx]
        
      }
      
    }
    
  }
  
  return(div_tetra_list)
  
}

#---------------------------------------
#3次元データのボロノイ分割を行う関数------
#ドロネー四面体分割のリストを受け取り、ボロノイ領域の頂点のリストを返す
voronoi_diagram3d<-function(div_tetra_list){
  
  #ボロノイの頂点と頂点のインデックスのリスト
  voronoi_verts<-lst()
  
  for(s in 1:length(div_tetra_list)){
    
    tri_1<-div_tetra_list[[s]]
    tri_1_verts<-lst(tri_1$p1, tri_1$p2, tri_1$p3, tri_1$p4)
    
    for (t in 1:length(div_tetra_list)) {
      
      if(s==t){next}
      
      tri_2<-div_tetra_list[[t]]
      
      cnt2<-0
      
      if(some(tri_1_verts, ~{all(. == tri_2$p1)})){cnt2<-cnt2+1}
      
      if(some(tri_1_verts, ~{all(. == tri_2$p2)})){cnt2<-cnt2+1}
      
      if(some(tri_1_verts, ~{all(. == tri_2$p3)})){cnt2<-cnt2+1}
      
      if(some(tri_1_verts, ~{all(. == tri_2$p4)})){cnt2<-cnt2+1}
      
      if(cnt2 == 3){
        
        if(none(voronoi_verts, ~{is.element(c(s,t), .$tetra_idx) %>% all()})){
          
          voronoi_verts[[length(voronoi_verts)+1]]<-rbind(tri_1$cntr, tri_2$cntr) %>% as_tibble() %>% 
            lst(voro_v = ., tetra_idx = c(s,t)) 
          
        }
      }
      
    }
    
  }
  
  return(voronoi_verts)
  
}

#----------------------------
#4次元データへの補間関数。データ全体に対して行う-----
#figure=補間対象データ、n_vics=近傍点の数
#補間された点を返す
voronoi_interpo4d<-function(figure, n_vics){
  
  element<-rep(0, length = nrow(figure))
  
  fig_dist<-dist(figure) %>% as.matrix()
  
  inted_points<-c()
  
  for (i in 1:nrow(figure)) {
    
    if(element[i]==0){
      
      vics<-get_vicinity(dist = fig_dist, center = i, nvic = n_vics)
      
      element[vics]<-element[vics]+1
      
      inted_points<-voronoi_border4d(figure, vics)[[1]] %>% rbind(inted_points, .)
      
    }
    
  }
  
  return(inted_points)
  
}

#------------------------------------------
#4次元データfigureにおいて、vicsで指定した近傍点に関して補間を行う関数-----
voronoi_border4d<-function(figure, vics){
  
  vics_pca<-prcomp(figure[vics, ])
  
  dela_tetra_list<-delaunay_triangulation3d(data = vics_pca$x[, 1:3])
  
  voro_verts<-sapply(dela_tetra_list, function(tetra){tetra$is_vertex(vics_pca$x[1, 1:3])}) %>% which() %>% 
    sapply(., function(idx)dela_tetra_list[[idx]]$cntr) %>% t()
  
  #vics_pcaの最大値以下、または最小値以上の要素を持つ点のみ残す
  vics_pca_max<-apply(vics_pca$x, 2, max)[1:3]
  vics_pca_min<-apply(vics_pca$x, 2, min)[1:3]
  
  voro_vert_idxs<-apply(voro_verts, 1, function(p){
    
    all(map2_lgl(p, vics_pca_max, ~{.x <= .y}), map2_lgl(p, vics_pca_min, ~{.x >= .y}))
    
  }) %>% which()
  
  vics_oricord<-(voro_verts[voro_vert_idxs, ] %*% t(vics_pca$rotation[, 1:3])) %>% 
    apply(., 1, function(p){p + figure[vics[1], ]}) %>% t()
  
  return(lst(vics_oricord=vics_oricord, voro_verts=voro_verts[voro_vert_idxs, ]))
  
}
