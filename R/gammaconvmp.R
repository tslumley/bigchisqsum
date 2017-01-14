
PRECISION<-100

setClass("gammaconvmp", slots=c(coef="mpfr",exp="mpfr",power="integer"))

setMethod("show","gammaconvmp",function(object) cat("Multiple-precision Gamma convolution with", length(object@exp),"terms\n"))


elgammamp<-function(alpha,beta){
    if ((length(alpha)!=1) || (length(beta)!=1)) stop("only one term allowed")
    if (!isWholeNumber(alpha)) stop("alpha must be an integer")
    beta<-mpfr(beta, PRECISION)

    new("gammaconvmp", coef=beta^alpha/gamma(alpha), exp = -beta, power=as.integer(alpha-1))
}


chisqtwomp<-function(a,b,k=1, N){
    if(!isWholeNumber(k)) stop("degree must be an integer")
    if(a<=0) stop("a must be positive")
    if(b<=0) stop("b must be positive")
     
    
    front<- ((4*a*b)^-(k/2)) * ((a-b)/(8*a*b))^(1/2-k/2)* (gamma(1/2+k/2)/gamma(k))
    
    Iv<-BesselExpand(k/2-1/2, (b-a)/(4*a*b),N )
    
    new("gammaconvmp",
        coef=mpfr(front*Iv$coefs, PRECISION),
        exp= mpfr(rep(-(a+b)/(4*a*b),length(Iv$coefs)), PRECISION),
        power=as.integer(round(Iv$powers+(k/2-1/2)))
        )
}


setMethod("*", c("numeric","gammaconvmp"), function(e1,e2) {
    if (any(e1<0)) stop("coefficients must be positive")
    new("gammaconv", coef=e2@coef*e1, exp=e2@exp,power=e2@power)
})


oneFonemp<-function(r,s, multiplier){
    if(!isWholeNumber(r)) stop("r must be an integer")
    if(!isWholeNumber(s)) stop("s must be an integer")
    if(!(r>0)) stop("must have r>0")
    if(!(r<s)) stop("must have r<s")

    multiplier<-mpfr(multiplier, PRECISION)
    
    front<-pochMpfr(1-s,r)*exp(lgamma(s-2+1)-lgamma(r-1+1))*multiplier^(1-s)
    frontpower<-1-s
    first<-list(); second<-list()
    first$powers<-0:(s-r-1)
    first$coef<-pochMpfr(-s+r+1,first$powers)/factorial(first$powers)/pochMpfr(2-s,first$powers)
    second$powers<-0:(r-1)
    second$coef<-(-1)^(second$powers)*pochMpfr(1-r,second$powers)/factorial(second$powers)/pochMpfr(2-s,second$powers)

    new("gammaconvmp",
        coef=front*c(first$coef*multiplier^first$powers,-(second$coef*multiplier^second$powers)),
        power=as.integer(c(first$powers+frontpower,second$powers+frontpower)),
        exp=rep(c(mpfr(0,PRECISION),multiplier),c(length(first$powers),length(second$powers)))
        )
    
}




convonemp<-function(a,n, b,m){ 
## $f(x)=\exp(ax)x^n$, $g(x)=\exp(bx)x^m$, $f\star g$
    front<- gamma(m+1)*gamma(n+1)/gamma(m+n+2)
    confluent<-oneFonemp(m+1,n+m+2, (b-a))
    new("gammaconvmp",coef=confluent@coef*front, power=confluent@power+m+n+1L,exp=confluent@exp+a)
}

setMethod("nrow","gammaconvmp",function(x) length(x@coef))


setMethod("%*%", c("gammaconvmp","gammaconvmp"), function(x,y) {
 

    maxterms<-sum(outer(x@power,y@power,function(n1,n2) 2*pmax(n1+1,n2+1)))
    allcoef<-mpfr(rep(0,maxterms),PRECISION)
    allpower<-integer(maxterms)
    allexp<-mpfr(rep(0,maxterms),PRECISION)
    here<-0
    
    for(i in 1:length(x@power)){
        for(j in 1:length(y@power)){
            if (x@exp[i]==y@exp[j]) {cat("#");next}
            
            term<-convonemp(x@exp[i],x@power[i],y@exp[j],y@power[j])
            size<-nrow(term)
            
            if(here+size>maxterms) stop("thomas can't count")
            allcoef[here+(1:size)]<-term@coef*x@coef[i]*y@coef[j]
            allpower[here+(1:size)]<-term@power
            allexp[here+(1:size)]<-term@exp
            here<-here+size
        }
    }
    
    new("gammaconvmp", coef=allcoef[1:here],power=allpower[1:here],exp=allexp[1:here])
})

                  
simplifymp<-function(object){
    i<-paste(asNumeric(object@exp),asNumeric(object@power))
    u<-!duplicated(i)
    ul<-unique(i)
    i<-factor(i,levels=unique(i))
    c<-mpfr(rep(0,sum(u)),PRECISION)
    for(j in 1:sum(u)){
      c[i]<-sum(object@coef[i==ul[j]])
    }
    new("gammaconvmp",coef=c, exp=object@exp[u],power=object@power[u])

}

                  

evaldensitymp<-function(object,x){
    if (length(x)==1)
        asNumeric(sum(object@coef*exp(object@exp*x)*x^(object@power)))
    else{
        e<-exp(outer(object@exp,x))
        p<-outer(object@power,x, function(a,b) b^a)
        asNumeric(colSums(e*p*object@coef))
    }
}


stcoefmp<- function(obj) {
    a<-obj@power+1
    asNumeric(obj@coef*gamma(a)/(-obj@exp)^(obj@power+1))
}

setMethod("as.data.frame","gammaconvmp",
          function(x) data.frame(coef=asNumeric(x@coef),exp=asNumeric(x@exp),power=asNumeric(x@power),stcoef=stcoefmp(x))
          )


setMethod("[",c("gammaconvmp","index","ANY","ANY"), 
	function(x,i,j,drop) new("gammaconvmp", coef=x@coef[i],power=x@power[i],exp=x@exp[i]) 
	)



	
sttailmp<- function(obj,x) {
    a<-obj@power+1
    b<- -obj@exp
    asNumeric(pgamma(x,a,b,lower=FALSE)*stcoef(obj))
}


evaltailmp<-function(object,x){
    s<-stcoef(object)
    if (length(x)==1)
        asNumeric(sum(sttail(object,x)))
    else{
    	  a<-object@power+1
   	  b<- -object@exp
   	  n<-length(a)
       p<-outer(1:n,x,function(j,xi) pgamma(xi,a[j],b[j],lower=FALSE))
       asNumeric(colSums(p*stcoef(object)))
    }
}



