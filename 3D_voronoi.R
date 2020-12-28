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
G2[3]<-G2[3]-r_i

#初期四面体の底面正三角形の一辺の中心点
H<-G2
H[2]<-H[2]-sqrt(2)*r_i

#データ点をすべて含む立方体の辺QC方向の単位ベクトル
#要するにx軸単位ベクトル
e_QC<-c(1, 0, 0)

#初期四面体の底面正三角形の頂点
S<-H+sqrt(6)*r_i*e_QC
T_v<-H-sqrt(6)*r_i*e_QC
U<-G2
U[2]<-U[2]+2*sqrt(2)*r_i

#初期四面体の底面正三角形に相対する頂点
V<-G2
V[3]<-V[3]+4*r_i

#初期四面体S-(T_v)-U-V

#初期四面体確認用プロット
plot3d(t_data)
points3d(rbind(Q, F_m, G), col=2)
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

M1<-rbind((S-T_v), (S-U), (U-V)) %>% solve()
vec1<-c((S%*%S - T_v%*%T_v), (S%*%S - U%*%U), (U%*%U - V%*%V))

r1<-(S-tetra_test$cntr)%*%(S-tetra_test$cntr) %>% sqrt()
spheres3d(x = G, radius = r1, col=4, alpha = 0.1)

#ドロネー分割
for (i in 1:nrow(t_data)) {
  
}


