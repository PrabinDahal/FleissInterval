##########################################################################################################################
## Title          : Function to work out confidence interval for correlated proportions using method by Fleiss et al, 2003
## Reference Book : Fleiss, Levin and Paik.Statistical methods for rates and proportions. 2003. pp. 441-444.
## Equations used : Equations 15.2, 15.3 and 15.4 on pages 441-444 of the book
## Version Date   : 08/05/15
## Example data   : pups surivival data from Table 15.1 (p443) from the book
###########################################################################################################################
require(binom) # For computing Wilson's confidence interval

## k      = Number of clusters
## i      = Cluster index
## yi     = Number of samples with outcome of interest in cluster i
## ni     = Total number of samples in cluster i
## vif    = variannce inflation factor
## data   = Data frame containing yi and ni. Each row represents a cluster

fleiss<-function(data,yi,ni){
    if(!is.null(data)){
    p<-sum(data$yi)/sum(data$ni)
    data$ni_2<-data$ni^2
    n.w=sum(data$ni_2)/sum(data$ni)
    K=nrow(data) 						          
    n.a=(sum(data$ni))/K
    data$var.com<-(data$ni-n.a)^2
    s2<-sum(data$var.com)/K
    } 
  j<-as.data.frame(cbind(p,K,n.a,s2))

  if(j$K>1){			
  # Use Fliess's if the aggregate data comes from >1 cluster else report standard interval using Wilson's method
      if(!is.null(j)){
      data$rho.com<-(data$yi*(data$yi-1))-(2*j$p*(data$ni-1)*data$yi)+(data$ni*(data$ni-1)*j$p*j$p) # Top of Equation (15.4)
      data$rho.com.dn<-data$ni*(data$ni-1)*j$p*(1-j$p)						    # Bottom of Equation (15.4)	
      rho.hat<-sum(data$rho.com)/sum(data$rho.com.dn)						    # Equation (15.4) evaluated
    } 
  
    l<-as.data.frame(cbind(j,rho.hat))
    if(!is.null(l)){
      vif<-1+((l$s2/l$n.a)+(l$n.a-1))*l$rho.hat  # Equation (15.2)
    } 
    
    m<-as.data.frame(cbind(l,vif))
    if(!is.null(m)){
      var.p<-(m$p*(1-m$p)*m$vif)/sum(data$ni)	# Equation (15.3)
    } 
    
    q<-as.data.frame(cbind(m,var.p))
    if(!is.null(q)){
      lower<-q$p-1.96*sqrt(q$var.p)
      upper<-q$p+1.96*sqrt(q$var.p)
    } 
    result<-as.data.frame(cbind(p,lower,upper,method="fleiss"))
    return(data.frame(result))
  }else{ 
    # Wilson method outputs are numeric, whereas the fleiss output are characters. Solution using binom package. 
    result<-data.frame(p=as.character(binom.confint(data$yi, data$ni, conf.level = 0.95, methods ="wilson")$mean),         
                       lower=as.character(binom.confint(data$yi, data$ni, conf.level = 0.95, methods = "wilson")$lower),
                       upper=as.character(binom.confint(data$yi, data$ni, conf.level = 0.95, methods = "wilson")$upper),
                       method=as.character(binom.confint(data$yi, data$ni, conf.level = 0.95, methods = "wilson")$method))
    return(merge(result,data.frame(method="wilson"),all=T))
  }
} # End;Not run
#####################################
## Application to pups data
#####################################
## Pups data from the book in Table 15.1 (p. 443)
pups<-data.frame(litter=c(seq(1,16)),
                 yi=c(12,11,10,9,10,9,9,8,8,4,7,4,5,3,3,0), 
                 ni=c(12,11,10,9,11,10,10,9,9,5,9,7,10,6,10,7)
		)
attach(pups)
## Apply the fuction
fleiss(pups,yi,ni)

# End (Not Run)
