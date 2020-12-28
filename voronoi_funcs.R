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
    
    #外接円の中心
    cntr = NA,
    
    #外接円の半径
    radius = NA,
    
    initialize = function(p1, p2, p3, p4){
      
      self$p1<-p1
      self$p2<-p2
      self$p3<-p3
      self$p4<-p4
      
      self$calc_center()
      self$calc_radius()
      
    },
    
    calc_center = function(){
      
      self$cntr<-(1/2) * solve( rbind((self$p1-self$p2), (self$p1-self$p3), (self$p3-self$p4)) ) %*% 
        c((self$p1%*%self$p1 - self$p2%*%self$p2), (self$p1%*%self$p1 - self$p3%*%self$p3), (self$p3%*%self$p3 - self$p4%*%self$p4)) %>% as.vector()
      
    },
    
    calc_radius = function(){
      
      self$radius<-(self$p1 - self$cntr) %*% (self$p1 - self$cntr) %>% sqrt()
      
    }
    
  )
  
)
