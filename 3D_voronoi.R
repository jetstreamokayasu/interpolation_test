#3次元ドロネー三角形分割

#----------------
#テストデータの準備------
t_x<-runif(10, min = -2, max = 2)
t_y<-runif(10, min = -2, max = 2)
t_z<-runif(10, min = -2, max = 2)

t_data<-cbind(t_x, t_y, t_z)

max_x<-max(t_x)*1.1
max_y<-max(t_y)*1.1
max_z<-max(t_z)*1.1

min_x<-min(t_x)*1.1
min_y<-min(t_y)*1.1
min_z<-min(t_z)*1.1

#データ点をすべて含む立方体の頂点
Q<-c(min_x, min_y, min_z)
F_m<-c(max_x, max_y, max_z)
C<-c(max_x, min_y, min_z)

#データ点をすべて含む立方体の中心
G<-(Q+F_m)/2

#データ点をすべて含む立方体の外接球の半径
#ドロネー分割の初期立方体の内接球とする
r_i<-sqrt((F_m-G) %*% (F_m-G))

#初期四面体の底面正三角形の重心
G2<-G
G2[3]<-G[3]-r_i

#初期四面体の底面正三角形の一辺の中心点
H<-G2
H[2]<-G2[2]-sqrt(2)*r_i

#データ点をすべて含む立方体の辺QC方向の単位ベクトル
#要するにx軸単位ベクトル
e_QC<-c(1, 0, 0)

#初期四面体の底面正三角形の頂点
S<-H+sqrt(6)*r_i*e_QC
T_v<-H-sqrt(6)*r_i*e_QC
U<-G2
U[2]<-G2[2]+2*sqrt(2)*r_i

#初期四面体の底面正三角形に相対する頂点
V<-G2
V[3]<-G2[3]+4*r_i

# V2<-G2
# V2[3]<-sqrt((S-T_v) %*% (S-T_v))*(sqrt(6)/3)+G2[3]

#初期四面体S-(T_v)-U-V

#点群をプロット
plot3d(t_data)
spheres3d(t_data, radius = 0.05)
#初期四面体確認用プロット
points3d(rbind(Q, F_m, G), col=2)
#points3d(rbind(G2, H), col=3)
points3d(rbind(G2, H, S, T_v, U, V), col=3)
spheres3d(x = G, radius = r_i, alpha = 0.1)
lines3d(rbind(S, T_v))
lines3d(rbind(U, T_v))
lines3d(rbind(S, U))
lines3d(rbind(S, V))
lines3d(rbind(T_v, V))
lines3d(rbind(U, V))

tetra_test<-Tetrahedron$new(S, T_v, U, V)
points3d(rbind(c(10, 10, 10), tetra_test$cntr), col=4)
spheres3d(x = tetra_test$cntr, radius = tetra_test$radius, col=4, alpha = 0.1)

tetra_test2<-Tetrahedron$new(t_data[1, ], t_data[2, ], t_data[3, ], t_data[4, ])

# M1<-rbind((S-T_v), (S-U), (U-V)) %>% solve()
# vec1<-c((S%*%S - T_v%*%T_v), (S%*%S - U%*%U), (U%*%U - V%*%V))
# 
# r1<-(S-tetra_test$cntr)%*%(S-tetra_test$cntr) %>% sqrt()
# spheres3d(x = G, radius = r1, col=4, alpha = 0.1)

#ドロネー分割を実現する四面体のリスト
div_tetra_list<-lst(tetra_test)

#tetra_listX<-lst(tetra_test, tetra_test2)

#ドロネー分割
for (i in 1:nrow(t_data)) {
  
  debugText(i)
  #追加するデータ点の選択
  select_p<-t_data[i, ]
  
  # 四面体の外接球内に追加したポイントが内包されていないか判定。
  temp_div_tetra<-lst()
  
  cnt<-0
  
  for(j in 1:length(div_tetra_list)){
    
    debugText(j)
    debugText(cnt)
    
    tri<-div_tetra_list[[j+cnt]]
    
    check<-tri$sphere_include_point(select_p)
    debugText(check)
    
    if(check){
      
      # 内包していた三角形の頂点を使用して再分割する。
      new_tetra1<-Tetrahedron$new(tri$p1, tri$p2, tri$p3, select_p)
      new_tetra2<-Tetrahedron$new(tri$p1, tri$p2, tri$p4, select_p)
      new_tetra3<-Tetrahedron$new(tri$p1, tri$p3, tri$p4, select_p)
      new_tetra4<-Tetrahedron$new(tri$p2, tri$p3, tri$p4, select_p)
      
      temp_div_tetra<-append(temp_div_tetra, lst(new_tetra1, new_tetra2, new_tetra3, new_tetra4))
      
      #debugText(temp_div_tetra)
      
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
  
  #debugText(del_tetras)
  
  for (d_tet in del_tetras) {
    
    eq_tetra<-sapply(temp_div_tetra, function(tetra)d_tet$is_equal_tetra(tetra))
    
    debugText(eq_tetra)
    
    if(any(eq_tetra)){
      
      eq_idx<-which(eq_tetra == TRUE)
      temp_div_tetra<-temp_div_tetra[-eq_idx]
      
    }
    
  }
  
  # 重複四面体を削除した後にappend
  div_tetra_list<-append(div_tetra_list, temp_div_tetra)
  
}

#初期四面体を含めたドロネー分割リスト
div_tetra_list_with_first_tetra<-div_tetra_list

#-------------------------
# はじめに作成した四面体の母点をもつ四面体分割を削除------
del_tetras<-lst()

for (tetra_check in div_tetra_list) {
  
  first_points<-lst(tetra_test$p1, tetra_test$p2, tetra_test$p3, tetra_test$p4)
  
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

#---------------------
#四面体描画----
for (m in 1:length(div_tetra_list)) {
  
  div_tetra_list[[m]]$plot_tetra_edges(edge_col = rainbow(length(div_tetra_list))[m])
  
}

t_num<-1

plot3d(t_data)
spheres3d(t_data, radius = 0.05)
writeWebGL(filename = "./webGL/3d_delaunay.html", width = 800, height = 800)

#------------------------------
#ボロノイ分割を求める----
voronoi_verts<-lst()

for(s in 1:length(div_tetra_list)){
  
  tri_1<-div_tetra_list[[s]]
  tri_1_verts<-lst(tri_1$p1, tri_1$p2, tri_1$p3, tri_1$p4)
  
  debugText(s)
  
  for (t in 1:length(div_tetra_list)) {
    
    debugText(t)
    
    if(s==t){next}
    
    tri_2<-div_tetra_list[[t]]
    
    cnt2<-0
    
    if(some(tri_1_verts, ~{all(. == tri_2$p1)})){cnt2<-cnt2+1}
    
    if(some(tri_1_verts, ~{all(. == tri_2$p2)})){cnt2<-cnt2+1}
    
    if(some(tri_1_verts, ~{all(. == tri_2$p3)})){cnt2<-cnt2+1}
    
    if(some(tri_1_verts, ~{all(. == tri_2$p4)})){cnt2<-cnt2+1}
    
    debugText(cnt2)
    
    if(cnt2 == 3){
      
      if(none(voronoi_verts, ~{is.element(c(s,t), .$tetra_idx) %>% all()})){
      
        voronoi_verts[[length(voronoi_verts)+1]]<-rbind(tri_1$cntr, tri_2$cntr) %>% as_tibble() %>% 
          lst(voro_v = ., tetra_idx = c(s,t)) 
      
      }
    }
    
  }
  
}

#初期四面体を含めた場合

voronoi_verts_with_first_tetra<-lst()

for(s in 1:length(div_tetra_list_with_first_tetra)){
  
  tri_1<-div_tetra_list_with_first_tetra[[s]]
  tri_1_verts<-lst(tri_1$p1, tri_1$p2, tri_1$p3, tri_1$p4)
  
  debugText(s)
  
  for (t in 1:length(div_tetra_list_with_first_tetra)) {
    
    debugText(t)
    
    if(s==t){next}
    
    tri_2<-div_tetra_list_with_first_tetra[[t]]
    
    cnt2<-0
    
    if(some(tri_1_verts, ~{all(. == tri_2$p1)})){cnt2<-cnt2+1}
    
    if(some(tri_1_verts, ~{all(. == tri_2$p2)})){cnt2<-cnt2+1}
    
    if(some(tri_1_verts, ~{all(. == tri_2$p3)})){cnt2<-cnt2+1}
    
    if(some(tri_1_verts, ~{all(. == tri_2$p4)})){cnt2<-cnt2+1}
    
    debugText(cnt2)
    
    if(cnt2 == 3){
      
      if(none(voronoi_verts_with_first_tetra, ~{is.element(c(s,t), .$tetra_idx) %>% all()})){
        
        voronoi_verts_with_first_tetra[[length(voronoi_verts_with_first_tetra)+1]]<-rbind(tri_1$cntr, tri_2$cntr) %>% as_tibble() %>% 
          lst(voro_v = ., tetra_idx = c(s,t)) 
        
      }
    }
    
  }
  
}

#----------------------------
#ボロノイ描画
lines3d(voronoi_verts[[1]])

#初期四面体あり
for (n in 1:length(voronoi_verts_with_first_tetra)) {
  lines3d(voronoi_verts_with_first_tetra[[n]]$voro_v, col = "red")
  #text3d(voronoi_verts[[n]]$voro_v, texts = n)
}

#初期四面体なし
for (n in 1:length(voronoi_verts)) {
  lines3d(voronoi_verts[[n]]$voro_v)
  #text3d(voronoi_verts[[n]]$voro_v, texts = n)
}

plot3d(t_data)
spheres3d(t_data, col = 4, radius = 0.05)
writeWebGL(width = 800, height = 800)

v_num<-24
t_num1<-voronoi_verts[[v_num]]$tetra_idx[1]
t_num2<-voronoi_verts[[v_num]]$tetra_idx[2]

div_tetra_list[[t_num1]]$plot_tetra_edges(edge_col = rainbow(16)[t_num1])
div_tetra_list[[t_num2]]$plot_tetra_edges(edge_col = rainbow(16)[t_num2])
lines3d(voronoi_verts[[v_num]]$voro_v)

#------------------------------
#関数確認------
div_tetra_list2<-delaunay_triangulation3d(data = t_data)
for (u in 1:length(div_tetra_list2)) {
  div_tetra_list2[[u]]$plot_tetra_edges(edge_col = rainbow(16)[u])
}

voronoi_verts2<-voronoi_diagram3d(div_tetra_list2)
map_lgl(voronoi_verts2, function(voro){some(voronoi_verts, ~{all( voro$voro_v == .$voro_v )})})

for (voro in voronoi_verts2) {
  lines3d(voro$voro_v)
}

#--------------------------
#指定した点を含む四面体のインデックスを抽出
tetra_idx<-sapply(div_tetra_list2, function(tetra){tetra$is_vertex(t_data[2, ])}) %>% which()

for (u in tetra_idx) {
  div_tetra_list2[[u]]$plot_tetra_edges(edge_col = rainbow(16)[u])
}

voro_vert<-sapply(tetzra_idx, function(idx)div_tetra_list2[[idx]]$cntr) %>% t() %>% as_tibble()


#------------------------
#rglアニメーション---

play3d(spin3d(axis = c(0, 1, 0), rpm = 5))
rgl_anime<-spin3d(axis = c(0, 0, 1), rpm = 3)
movie3d(rgl_anime, duration = 30, movie = "delaunay", dir = "./pics", dev = 1, fps = 10)

