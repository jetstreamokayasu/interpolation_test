require(phacm)

anu<-anulusUnif(100, 1, 1.5)
plot3d(anu)
anu2<-anulusUnif(100)
plot(anu2)
plot(anu)

anu.diag<-compute_pd(anu, 1, 2)
plot(anu.diag)
autoplot(anu.diag)

anu.pl<-compute_pl(anu.diag)
autoplot(anu.pl)
anu.smpl<-compute_smooth_pl(anu.pl)
autoplot(anu.smpl, dimension = 1)

sphe<-sphereUnif(100, 2, 1)
figurePlot(sphe)

sphe.diag<-compute_pd(sphe, 2, 2)
autoplot(sphe.diag)

sphe.pl<-compute_pl(sphe.diag)
autoplot(sphe.pl)

sphe.smpl<-compute_smooth_pl(sphe.pl)
autoplot(sphe.smpl, dimension = 1)
local.max<-count_local_maximal(sphe.pl, 0.1)
losm.ma<-count_smooth_maximal(sphe.pl)

torus200<-torusUnif(200, 1, 2.5)
figurePlot(torus200)

torus200.diaga<-compute_pd(torus200, maxdimension = 2, maxscale = 3)
autoplot(torus200.diaga)

torus200.pl<-compute_pl(torus200.diaga)
autoplot(torus200.pl)

trus200.smpl<-compute_smooth_pl(torus200.pl)
autoplot(trus200.smpl)

trus200.betti<-count_smooth_maximal(torus200.pl)

torus200.in<-voronoiInterpo(figure = torus200, nvics = 15)
points3d(torus200.in, col="orange")
