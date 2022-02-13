getari<-function(r1,r2,maxdist=200)
{	l=max(c(r1[,2],r2[,2]));
	A=B=array(0,dim=rep(l,2));
	for(i in order(r1[,2]-r1[,1],decreasing=T))
	{	A[r1[i,1]:r1[i,2],r1[i,1]:r1[i,2]]=i;
	}
	for(i in order(r2[,2]-r2[,1],decreasing=T))
	{	B[r2[i,1]:r2[i,2],r2[i,1]:r2[i,2]]=i;
	}
	a=b=ari=NULL;
	for(i in 0:min(maxdist,l-1))
	{	ta=A[(1:(l-i)-1)*l+(i+1):l];
		tb=B[(1:(l-i)-1)*l+(i+1):l];
		ari=c(ari,myrand(ta,tb)[2]);
		a=c(a,ta);
		b=c(b,tb);
	}
	rt=NULL;
	rt$ari=myrand(a,b);
	rt$dari=ari;

	return(rt);
}

#FP and FN are in the sense that x is the truth
myrand<-function(x,y)
{
	t=as.matrix(table(x,y));
	a=sum(choose(t,2));#TP
	c=sum(choose(apply(t,1,sum),2)) - a; #FN
	d=sum(choose(apply(t,2,sum),2)) -  a; #FP
	b=sum(choose(sum(t),2))-a-c-d;#TN
a=a+1;
b=b+1;
c=c+1;
d=d+1;

	ri=(a+b)/(a+b+c+d);
	ari=(a-(a+c)*(a+d)/(a+b+c+d))/((a+c+a+d)/2-(a+c)*(a+d)/(a+b+c+d));
	jac=a/(a+c+d);
	fm=sqrt(a/(a+d)*a/(a+c));
	return(c(ri,ari,jac,fm));#,a-1,b-1,c-1,d-1));
}

getadjr2<-function(hic,tad,maxdist=200)
{	l=dim(hic)[1];
	print(dim(tad))
	A=array(0,dim=rep(l,2));
	for(i in order(tad[,2]-tad[,1],decreasing=T))
	{	A[tad[i,1]:tad[i,2],tad[i,1]:tad[i,2]]=i;
	}
	adjr2=pr=rc=oh=NULL;
	for(i in 0:min(maxdist,l-1))
	{	t=(1:(l-i)-1)*l+(i+1):l;
		tt=unique(A[t]);
		#print (max(tt))
		#print (dim(tad)[1])
		#print (length(unique(tt)))
		y=hic[t];
		mm=rep(0,max(tt+1));
		for(j in tt)
		{	mm[j+1]=mean(y[A[t]==j]); 
		}
		yy=y-mm[A[t]+1];
		#adjr2=rbind(adjr2,c(1.-var(yy)/var(y)*(length(t)-1)/(length(t)-length(tt)),var(y),length(tt)));
		adjr2=rbind(adjr2,c(1.-sum(yy^2)/sum(y^2)*(length(t)-1)/(length(t)-length(tt)),var(y),length(tt)));
		#print(adjr2[i+1,]);
		cut=quantile(y,prob=1:19/20)
		tpr=trc=toh=NULL;
		for(j in cut)
		{	n=length(which(y>=j & A[t]>0));
			trc=c(trc,n/length(which(y>=j)));
			toh=c(toh,(length(which(y>=j & A[t]==0))+1)/(length(which(A[t]==0))+1));
			tpr=c(tpr,n/length(which(A[t]>0)));
		}
		pr=rbind(pr,tpr);
		rc=rbind(rc,trc);
		oh=rbind(oh,toh);
	}
	rt=NULL;
	rt$adjr2=adjr2;
	rt$pr=pr;
	rt$rc=rc;
	rt$oh=oh;
	return(rt);
}

getadjr3<-function(hic,tad,maxdist=200)
{       l=dim(hic)[1];
A=array(0,dim=rep(l,2));
tad = tad[order(tad[,2]-tad[,1],decreasing=T),]
for(i in 1:dim(tad)[1])
{       A[tad[i,1]:tad[i,2],tad[i,1]:tad[i,2]]=i;
}
for(i in order(tad[,2]-tad[,1],decreasing=T))
{       A[tad[i,1]:tad[i,2],tad[i,1]:tad[i,2]]=A[tad[i,1]:tad[i,2],tad[i,1]:tad[i,2]]+dim(tad)[1];
}
adjr2=pr=rc=oh=NULL;
for(i in 0:min(maxdist,l-1))
{       t=(1:(l-i)-1)*l+(i+1):l;
tt=sort(unique(A[t]));
y=hic[t];
mm=rep(0,max(tt+1));
                for(j in tt)
                {       mm[j+1]=mean(y[A[t]==j]);
                }
                yy=y-mm[A[t]+1];
#adjr2=rbind(adjr2,c(1.-var(yy)/var(y)*(length(t)-1)/(length(t)-length(tt)),var(y),length(tt)));
adjr2=rbind(adjr2,c(1.-sum(yy^2)/sum(y^2)*(length(t)-1)/(length(t)-length(tt)),var(y),length(tt)));

#print(adjr2[i+1,]);
cut=quantile(y,prob=1:19/20)
tpr=trc=toh=NULL;
for(j in cut)
{       n=length(which(y>=j & A[t]>0));
trc=c(trc,n/length(which(y>=j)));
toh=c(toh,(length(which(y>=j & A[t]==0))+1)/(length(which(A[t]==0))+1));
tpr=c(tpr,n/length(which(A[t]>0)));
}
pr=rbind(pr,tpr);
rc=rbind(rc,trc);
oh=rbind(oh,toh);
}
rt=NULL;
rt$adjr2=adjr2;
rt$pr=pr;
rt$rc=rc;
rt$oh=oh;
return(rt);
}
