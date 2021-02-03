#スライドや論文用の図表を作る

#---------------------
#PCAの再構築先を平均点を通る面ではなく、中心点を通る面にしたことを表す図を作る-------

thetes<-seq(0, 8*pi/10, by = pi/10)

#中心点＆近傍点
#まずすべて緑で描画
plot(cos(thetes+pi/2), sin(thetes+pi/2), xlim = c(-1.2, 1.2), ylim = c(-1.2, 1.2), axes = F, xlab = "", ylab = "", 
     pch = 16, cex = 2.5, col = 3)
#近傍点以外の点を黒でプロット
points(cos(c(-pi/10, 9*pi/10)+pi/2), sin(c(-pi/10, 9*pi/10)+pi/2), pch = 16, cex = 2.5)

circ_pca<-prcomp(cbind(cos(thetes+pi/2), sin(thetes+pi/2)))

#中心点
central<-c(cos(4*pi/10+pi/2), sin(4*pi/10+pi/2))
points(cos(4*pi/10+pi/2), sin(4*pi/10+pi/2), cex = 2.5, pch = 16, col = "blue")
#points(circ_pca$center[1], circ_pca$center[2], col = 2)

#PCAの再構築先の面
draw_line(x = circ_pca$center+5*circ_pca$rotation[, "PC1"], y = circ_pca$center-5*circ_pca$rotation[, "PC1"], lwd = 4, col = "blueviolet")

#PCAの再構築先の面上に分布する元の点を描画
circ_pc1<-circ_pca$x[,1] %*% t(circ_pca[["rotation"]][, "PC1"]) %>% 
  apply(., 1, function(x){x+circ_pca$center}) %>% t() 
points(circ_pc1[,1], circ_pc1[,2], pch = 16, col = "blueviolet")

#擬似的な補間点を赤で描画
circ_added<-rbind((circ_pc1[5,]+circ_pc1[4,])/2, (circ_pc1[5,]+circ_pc1[6,])/2)
points(x = circ_added[, 1], y = circ_added[, 2], pch = 16, cex = 2.5, col = "red")

#接平面を描画
draw_line(x = c(cos(4*pi/10+pi/2), sin(4*pi/10+pi/2))+3*c(-sin(4*pi/10+pi/2), cos(4*pi/10+pi/2)), 
          y = c(cos(4*pi/10+pi/2), sin(4*pi/10+pi/2))-3*c(-sin(4*pi/10+pi/2), cos(4*pi/10+pi/2)), lwd = 3)

#--------------------------
#PCAの再構築先の面を、中心点を通るように平行移動-----
draw_line(x = c(cos(4*pi/10+pi/2), sin(4*pi/10+pi/2))+5*circ_pca$rotation[, "PC1"], y = c(cos(4*pi/10+pi/2), sin(4*pi/10+pi/2))-5*circ_pca$rotation[, "PC1"], lwd = 4, col = "blueviolet")

#平行移動したPCA再構築面上に分布する元の点を描画
circ_pc1_moved<-circ_pca$x[,1] %*% t(circ_pca[["rotation"]][, "PC1"]) %>% 
  apply(., 1, function(x){x+central}) %>% t() 
points(circ_pc1_moved[,1], circ_pc1_moved[,2], pch = 16, col = "blueviolet")

#擬似的な補間点を赤で描画
circ_added2<-rbind((circ_pc1_moved[5,]+circ_pc1_moved[4,])/2, (circ_pc1_moved[5,]+circ_pc1_moved[6,])/2)
points(x = circ_added2[, 1], y = circ_added2[, 2], pch = 16, cex = 2.5, col = "red")
