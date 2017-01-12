##
## convolution algebra for chisquare convolutions, from https://arxiv.org/pdf/1208.2691v3.pdf
##


setClass("gammaconv", slots=c(coef="numeric",exp="numeric",power="numeric"))

setMethod("show","gammaconv",function(object) cat("Gamma convolution with", length(object@exp),"terms\n"))


isWholeNumber<-function(n) all((n-round(n))<1e-8)



## beta is inverse shape parameter: kernel is $x^{\alpha-1}\exp(-\beta x)$

elgamma<-function(alpha,beta){
    if ((length(alpha)!=1) || (length(beta)!=1)) stop("only one term allowed")
    if (!isWholeNumber(alpha)) stop("alpha must be an integer")

    new("gammaconv", coef=beta^alpha/gamma(alpha), exp = -beta, power=alpha-1)
}


BesselExpand<-function(nu, coef, truncation){
    m<-(0:(truncation/2))
    m<-m[2*m+nu<=truncation]

    a<-(coef/2)^(2*m+nu)/gamma(m+nu+1)/factorial(m)
    list(powers=2*m+nu,coefs=a)
}


BesselError<-function(N,a,b) {
	if(a<b){
		temp<-a;a<-b;b<-temp
	}
	lambda<-(a-b)/(4*a*b)
	Iv<- BesselExpand(0,lambda,N)
	g<-function(x) outer(x,Iv$powers,"^")%*%Iv$coefs
	function(x) besselI(x*lambda,0,expon.scaled=TRUE)-g(x)*exp(-lambda*x)
	}


OptM<-function(a,b,xmax,epsilon,memo=NA){
	lambda<-(a-b)/(4*a*b)
	N1<-if(is.na(memo)) 10*abs(log10(epsilon)) else memo
	N0<-5
	lowError <-integrate(BesselError(N0,a,b),lower=0,upper=xmax)$value
	if (abs(lowError)<epsilon) 
	    return(N0)
	repeat({ cat("!")
	    hiError <-integrate(BesselError(N1,a,b),lower=0,upper=xmax)$value
	   if (abs(hiError)>epsilon)
	     N1<-N1+5
	   else
	     break
	})
	repeat({ 
	     N<-round((N1+N0)/2)
	     cat(N)
	     if((N1==N) || (N0==N)) 
	         break
	     midError <-integrate(BesselError(N,a,b),lower=0,upper=xmax)$value
	     if (abs(midError)<epsilon)
	       N1<-N
	     else{
	     	 N0<-N
	     }   
	})     
	return(N)	
	
}


chisqtwo<-function(a,b,k=1, N){
    if(!isWholeNumber(k)) stop("degree must be an integer")
    if(a<=0) stop("a must be positive")
    if(b<=0) stop("b must be positive")
     
    
    front<- ((4*a*b)^-(k/2)) * ((a-b)/(8*a*b))^(1/2-k/2)* (gamma(1/2+k/2)/gamma(k))
    
    Iv<-BesselExpand(k/2-1/2, (b-a)/(4*a*b),N )
    
    new("gammaconv",
        coef=front*Iv$coefs,
        exp= rep(-(a+b)/(4*a*b),length(Iv$coefs)),
        power=round(Iv$powers+(k/2-1/2))
        )
}


setMethod("*", c("numeric","gammaconv"), function(e1,e2) {
    if (any(e1<0)) stop("coefficients must be positive")
    new("gammaconv", coef=e2@coef*e1, exp=e2@exp,power=e2@power)
})

pochhammer1<-function(x,n) prod(x-1+seq_len(n))
pochhammer<-Vectorize(pochhammer1)

oneFone<-function(r,s, multiplier){
    if(!isWholeNumber(r)) stop("r must be an integer")
    if(!isWholeNumber(s)) stop("s must be an integer")
    if(!(r>0)) stop("must have r>0")
    if(!(r<s)) stop("must have r<s")
    
    front<-pochhammer(1-s,r)*exp(lgamma(s-2+1)-lgamma(r-1+1))*multiplier^(1-s)
    frontpower<-1-s
    first<-list(); second<-list()
    first$powers<-0:(s-r-1)
    first$coef<-pochhammer(-s+r+1,first$powers)/factorial(first$powers)/pochhammer(2-s,first$powers)
    second$powers<-0:(r-1)
    second$coef<-(-1)^(second$powers)*pochhammer(1-r,second$powers)/factorial(second$powers)/pochhammer(2-s,second$powers)

    new("gammaconv",
        coef=front*c(first$coef*multiplier^first$powers,-(second$coef*multiplier^second$powers)),
        power=c(first$powers+frontpower,second$powers+frontpower),
        exp=rep(c(0,multiplier),c(length(first$powers),length(second$powers)))
        )
    
}


convone<-function(a,n, b,m){ 
## $f(x)=\exp(ax)x^n$, $g(x)=\exp(bx)x^m$, $f\star g$
    front<- gamma(m+1)*gamma(n+1)/gamma(m+n+2)
    confluent<-oneFone(m+1,n+m+2, (b-a))
    new("gammaconv",coef=confluent@coef*front, power=confluent@power+m+n+1,exp=confluent@exp+a)
}

setMethod("%*%", c("gammaconv","gammaconv"), function(x,y) {

    allcoef<-numeric(0)
    allpower<-integer(0)
    allexp<-numeric(0)
    
    for(i in 1:length(x@power)){
        for(j in 1:length(y@power)){
            if (x@exp[i]==y@exp[j]) {cat("#");next}
            term<-convone(x@exp[i],x@power[i],y@exp[j],y@power[j])
            allcoef<-c(allcoef,term@coef*x@coef[i]*y@coef[j])
            allpower<-c(allpower,term@power)
            allexp<-c(allexp,term@exp)
        }
    }
    
    new("gammaconv", coef=allcoef,power=allpower,exp=allexp)
})


simplify<-function(object){
    i<-paste(object@exp,object@power)
    u<-!duplicated(i)
    c<-rowsum(object@coef,i,reorder=FALSE)
    new("gammaconv",coef=as.vector(c), exp=object@exp[u],power=object@power[u])

}

          

evaldensity<-function(object,x){
    if (length(x)==1)
        sum(object@coef*exp(object@exp*x)*x^(object@power))
    else{
        e<-exp(outer(object@exp,x))
        p<-outer(object@power,x, function(a,b) b^a)
        colSums(e*p*object@coef)
    }
}



stcoef<- function(obj) {
    a<-obj@power+1
    obj@coef*gamma(a)/(-obj@exp)^(obj@power+1)
}

setMethod("as.data.frame","gammaconv", function(x) data.frame(coef=x@coef,exp=x@exp,power=x@power,stcoef=stcoef(x)))


setMethod("[",c("gammaconv","index","ANY","ANY"), 
	function(x,i,j,drop) new("gammaconv", coef=x@coef[i],power=x@power[i],exp=x@exp[i]) 
	)
	
	
	
sttail<- function(obj,x) {
    a<-obj@power+1
    b<- -obj@exp
    pgamma(x,a,b,lower=FALSE)*stcoef(obj)
}


evaltail<-function(object,x){
    s<-stcoef(object)
    if (length(x)==1)
        sum(sttail(object,x))
    else{
    	  a<-object@power+1
   	  b<- -object@exp
   	  n<-length(a)
       p<-outer(1:n,x,function(j,xi) pgamma(xi,a[j],b[j],lower=FALSE))
       colSums(p*stcoef(object))
    }
}

