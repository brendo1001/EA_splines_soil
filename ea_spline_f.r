#EA_spline function
ea_spline<- function(mxd,dat,lam,d){
ndata=length(levels(as.factor(dat[,1])))

# matrix of the continous splines for each data point
m_fyfit<- matrix(NA,ncol=length(c(0:mxd)),nrow=ndata)

# Matrix in which the averaged values of the spline are fitted. The depths are specified in the (d) object
yave<- matrix(NA,ncol=length(d),nrow=ndata)

# Matrix in which the sum of square errors of each lamda iteration for the working profile are stored
sse<- matrix(NA,ncol=length(lam),nrow=1)

# Matrix in which the sum of square errors for eac h lambda iteration for each profile are stored
sset<- matrix(NA,ncol=length(lam),nrow=ndata)

# Matrix profile ids
mat_id<- matrix(NA,ncol=1,nrow=ndata)

# number of data


nl<- length(lam)  # Length of the lam matrix

s<- 0.05*mean(dat[,4])  # 5% of the standard deviation of the target attribute 
s2= s*s   # overall variance of soil attribute
    ##########################################################################################################  
d.dat<- as.data.frame(dat)
sp_dat<-split(d.dat,d.dat[,1]) 
for (st in 1:ndata) {
  
  subs<-sp_dat[[st]]  # subset the profile required
  subs<-as.matrix(subs)
  mat_id[st,1]<- subs[1,1]


    # manipulate the profile data to the required form
    ir<- c(1:length(subs[,1]))
    ir<-as.matrix(t(ir))
    u<- subs[ir,2]
    u<-as.matrix(t(u))   # upper 
    v<- subs[ir,3]
    v<-as.matrix(t(v))   # lower
    y<- subs[ir,4]
    y<-as.matrix(t(y))   # concentration 
    n<- length(y);       # number of observations in the profile
      
      ### routine for handling profiles with one observation
      
      if (n == 1)
      {xfit<- as.matrix(t(c(0:mxd)))
      nj<- max(v)+1
      if (nj > mxd)
      {nj<- mxd+1}
      yfit<- xfit 
      yfit[,1:nj]<- y   # values extrapolated onto yfit
      if (nj <= mxd)
              {yfit[,(nj+1):(mxd+1)]=-9999}
              m_fyfit[st,]<- yfit
              nd<- length(d)-1  # number of depth intervals
          dl<-d+1     #  increase d by 1
          
          for (cj in 1:nd) {
                            xd1<- dl[cj]
                            xd2<- dl[cj+1]     
                            if (nj>xd1 & nj<xd2)
                            {xd2<- nj
                            yave[st,cj]<- mean(yfit[,xd1:xd2])}
                            else
                            {yave[st,cj]<- mean(yfit[,xd1:xd2])}   # average of the yfit at the specified depth intervals
                            yave[st,cj+1]<- max(v)
                                              }
               }
      # End of single observation profile routine
      
      # Start of routine for fitting spline to profiles with multiple observations         
      
      else                                      
      {np1<- n+1;        # number of interval boundaries
      nm1<- n-1;
      delta<- v-u;      # depths of each layer
      del<- c(u[2:n],u[n])-v      # del is (u1-v0,u2-v1, ...)

      ####################################################################################################
      #create the (n-1)x(n-1) matrix r
      #first create r with 1's on the diagonal and upper diagonal, and
      # 0's elsewhere
      r<-matrix(0,ncol=nm1,nrow=nm1)
      for (dig in 1:nm1){
      r[dig,dig]<-1}
      for (udig in 1:nm1-1){
      r[udig,udig+1]<-1}

      #then create a diagonal matrix d2 of differences to premultiply
      #the current r
      d2<- matrix(0,ncol=nm1,nrow=nm1)
      diag(d2)<-delta[2:n]  # delta = depth of each layer

      # then premultiply and add the transpose; this gives half of r
      r<- d2 %*% r
      r<- r + t(r)

      # then create a new diagonal matrix for differences to add to the diagonal
      d1<- matrix(0,ncol=nm1,nrow=nm1)
      diag(d1)<-delta[1:nm1]  # delta = depth of each layer

      d3<- matrix(0,ncol=nm1,nrow=nm1)
      diag(d3)<-del[1:nm1]  # del =  differences

      r= r+2*d1 + 6*d3;    ###

      # create the (n-1)xn matrix q
      q<- matrix(0,ncol=n,nrow=n)
      for (dig in 1:n){
      q[dig,dig]<- -1 }
      for (udig in 1:n-1){
      q[udig,udig+1]<-1 }
      q<- q[1:nm1,1:n]
      dim.mat<-matrix(q[],ncol=length(1:n),nrow=length(1:nm1))

      # inverse of r
      rinv<- solve(r)

      # identity matrix i
      ind<- diag(n)

      # create the matrix coefficent z

      pr.mat<- matrix(0,ncol=length(1:nm1),nrow=length(1:n))
      pr.mat[]<-6*n*lam
      fdub<-pr.mat*t(dim.mat)%*%rinv
      z<-fdub%*%dim.mat+ind

      # solve for the fitted layer means
      sbar<- solve(z,t(y))

      # calculate the fitted value at the knots
      b<- 6*rinv%*%dim.mat%*% sbar
      b0<- rbind(0,b) # add a row to top = 0
      b1<- rbind(b,0) # add a row to bottom = 0
      gamma<-(b1-b0) / t(2*delta)
      alfa<- sbar-b0 * t(delta) / 2-gamma * t(delta)^2/3


      # fit the spline 
      xfit<- as.matrix(t(c(0:mxd))) # spline will be interpolated onto these depths (1cm res)
      nj<- max(v)+1
      if (nj > mxd)
      {nj<- mxd+1}
      yfit<- xfit
        for (k in 1:nj){
        xd<-xfit[k]
        if (xd< u[1])
          {p<- alfa[1]} else
          {for (its in 1:n) {
          if(its < n)
          {tf2=as.numeric(xd>v[its] & xd<u[its+1])} else {tf2<-0}
          if (xd>=u[its] & xd<=v[its])
              {p=alfa[its]+b0[its]*(xd-u[its])+gamma[its]*(xd-u[its])^2} else if (tf2)
              {phi=alfa[its+1]-b1[its]*(u[its+1]-v[its])
              p=phi+b1[its]*(xd-v[its])}
                            }}
              yfit[k]=p }
              if (nj <= mxd)
              {yfit[,nj:mxd+1]=-9999}
              m_fyfit[st,]<- yfit
        ##########################################################################################################
          # CALCULATION OF THE ERROR BETWEEN OBSERVED AND FITTED VALUES
          # calculate Wahba's estimate of the residual variance sigma^2
          ssq<- sum((t(y)-sbar)^2)
          g<- solve(z)
          ei<- eigen(g)
          ei<- ei$values
          df<- n-sum(ei)
          sig2w<- ssq/df
          # calculate the Carter and Eagleson estimate of residual variance
          dfc<- n-2*sum(ei)+sum(ei^2)
          sig2c<- ssq/dfc
          # calculate the estimate of the true mean squared error
          tmse= ssq/n-2*s2*df/n+s2
          sset[st]<- tmse
        ##########################################################################################################  
          
          # Averages of the spline at specified depths
          nd<- length(d)-1  # number of depth intervals
          dl<-d+1     #  increase d by 1
          for (cj in 1:nd) { 
                            xd2<- dl[cj+1]  
                            if (dl[cj] == 1)
                            {xd1<- 1}
                            else
                            {xd1<- dl[cj]+1} 
                            xd2<- dl[cj+1]     
                            
                            if (nj>=xd1 & nj<=xd2)
                            {xd2<- nj
                            yave[st,cj]<- mean(yfit[,xd1:xd2])}
                            else
                            {yave[st,cj]<- mean(yfit[,xd1:xd2])}   # average of the spline at the specified depth intervals
                            yave[st,cj+1]<- max(v)
          ##########################################################################################################  
                                              }
          }}     
          yave_d<- as.data.frame(yave)
          jmat<- matrix(NA,ncol=1,nrow=length(d))
          for (i in 1:length(d)-1) {
          a1<-paste(d[i],d[i+1],sep="-")
          a1<-paste(a1,"cm",sep="")
          jmat[i]<- a1 }
          jmat[length(d)]<- "soil depth"
          for (jj in 1:length(jmat)){
          names(yave_d)[jj]<- jmat[jj] 
          }
          yave_d<-data.frame(profile_id=mat_id, yave_d)
          retval <- list(sset,yave_d, t(m_fyfit))
                    return(retval)}
     names(z)[3] <- "c2"               
# END


# Function to convert alphabet types profile lables to numeric lables
# Start
alp_num <- function (dat_m) {
c1<-as.numeric(dat_m[,1])
c2<-as.numeric(dat_m[,2])
c3<-as.numeric(dat_m[,3])
c4<-as.numeric(dat_m[,4])
cc<-as.matrix(cbind(c1, c2,c3,c4))
return (cc) }
# END          
          

