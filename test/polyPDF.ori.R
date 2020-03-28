polyPDF<-function(moments,xmin=NULL,xmax=NULL,scale=1,disp=FALSE){
  #This function constructs a PD function from the "moments" input data frame,
  #containing the moment order ("n") in one column and their values ("moments")
  #in another. The limits of the function ("xmin","xmax") are required inputs.
  #A scale factor ("scale") can be used to transform data avoiding singularity.
  #By default, "scale" is 1. If "disp" is TRUE, the PDF will be plotted and the
  #polynomial coefficients are shown. The output of this function is the PDF.
  #Moments can be generated using the "momentset" function given in Appendix A.2.
  if (is.null(xmin)&is.null(xmax)){
    print("Please input xmin and xmax estimated values")
    return(NULL)
  }
  
  n=length(moments$n) #Number of moments available
  degree=n-1 #Degree of polynomial
  A=matrix(0,n,n) #Initialize matrix of coefficients
  B=matrix(0,n,1) #Initialize vector of independent terms
  for (i in 1:n){ #For each moment
    ni=moments$n[i] #ni-th moment
    for (j in 1:n){ #For each power term
      A[i,j]=(((xmax*scale)^(ni+j))-((xmin*scale)^(ni+j)))/(ni+j)
    }
    B[i]=moments$moments[i]*scale^ni #Scale moments
  }
  a=as.vector(solve(A,B)) #Find coefficients
  #Definition of PDF function
  rhof<-function(x){
    rho=a[1] #Independent term
    for (i in 2:(degree+1)){
      rho=rho+a[i]*(x*scale)^(i-1) #Polynomial terms
    }
    #Density is zero when rho is negative or X is beyond boundaries
    rho=rho*scale*as.integer(x>=xmin&x<=xmax)*as.integer(rho>0)
    return(rho)
  }
  
  if (disp==TRUE){
    #Set coefficients names
    names(a)[1:2]=c("(Intercept)","x")
    if (length(a)>2){
      for (i in 3:length(a)){
        names(a)[i]=paste("x^",toString(i-1))
      }
    }
    print("Polynomial model of the PDF:")
    print(a)
    #PDF plot
    i=1:1001
    y=xmin+(xmax-xmin)*(i-1)/1000
    rho=rhof(y) #Calculate density
    plot(y,rho,type="l",col="blue",xlim=c(xmin,xmax),ylim=c(0,max(rho)),
         xlab="Measurement values",ylab="Probability Density")
  }
  return(rhof)
}