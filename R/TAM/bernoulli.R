p1=0.3
p2=0.4
p3=0.5
p=c(p1,p2,p3)
a<-expand.grid(0:1,0:1,0:1)
b=rep(0,2^3)
b[which(apply(a,1,sum)==0)]=1
b[which(apply(a,1,sum)==1)]=0
b[which(apply(a,1,sum)==2)]=0.1125
b[which(apply(a,1,sum)==3)]=0
c<-kronecker(matrix(c(1-p3,p3,-1,1),nrow=2),kronecker(matrix(c(1-p2,p2,-1,1),nrow=2),matrix(c(1-p1,p1,-1,1),nrow=2)))%*%b
c


p1=0.3
p2=0.4
p3=0.5
p=c(p1,p2,p3)
a<-expand.grid(0:1,0:1,0:1)
b=rep(0,2^3)
b[which(apply(a,1,sum)==0)]=1
b[which(apply(a,1,sum)==1)]=0
b[which(apply(a,1,sum)==2)]=c(0.10,0.10,0.19)
b[which(apply(a,1,sum)==3)]=0
c<-kronecker(matrix(c(1-p3,p3,-1,1),nrow=2),kronecker(matrix(c(1-p2,p2,-1,1),nrow=2),matrix(c(1-p1,p1,-1,1),nrow=2)))%*%b
c

d<-array(c,dim=c(2,2,2))
d[,,1]+d[,,2]

a<-expand.grid(0:1,0:1)
b=rep(0,2^2)
b[which(apply(a,1,sum)==0)]=1
b[which(apply(a,1,sum)==2)]=p2*(1-p3)
kronecker(matrix(c(1-p3,p3,-1,1),nrow=2),matrix(c(1-p2,p2,-1,1),nrow=2))%*%b



p1=0.5
p2=0.5
p3=0.5
p=c(p1,p2,p3)
a<-expand.grid(0:1,0:1,0:1)
b=rep(0,2^3)
b[which(apply(a,1,sum)==0)]=1
b[which(apply(a,1,sum)==1)]=0
b[which(apply(a,1,sum)==2)]=0.25
b[which(apply(a,1,sum)==3)]=0
c<-kronecker(matrix(c(1-p3,p3,-1,1),nrow=2),kronecker(matrix(c(1-p2,p2,-1,1),nrow=2),matrix(c(1-p1,p1,-1,1),nrow=2)))%*%b
c


p1=0.4
p2=0.4
p3=0.4
p=c(p1,p2,p3)
a<-expand.grid(0:1,0:1,0:1)
b=rep(0,2^3)
b[which(apply(a,1,sum)==0)]=1
b[which(apply(a,1,sum)==1)]=0
b[which(apply(a,1,sum)==2)]=0.4*(0.6)*0.75
b[which(apply(a,1,sum)==3)]=0
c<-kronecker(matrix(c(1-p3,p3,-1,1),nrow=2),kronecker(matrix(c(1-p2,p2,-1,1),nrow=2),matrix(c(1-p1,p1,-1,1),nrow=2)))%*%b
array(c,dim=c(2,2,2))


p1=0.2
p2=0.2
p3=0.2
p=c(p1,p2,p3)
a<-expand.grid(0:1,0:1,0:1)
b=rep(0,2^3)
b[which(apply(a,1,sum)==0)]=1
b[which(apply(a,1,sum)==1)]=0
b[which(apply(a,1,sum)==2)]=0.2*(0.8)*0.5714
b[which(apply(a,1,sum)==3)]=0
c<-kronecker(matrix(c(1-p3,p3,-1,1),nrow=2),kronecker(matrix(c(1-p2,p2,-1,1),nrow=2),matrix(c(1-p1,p1,-1,1),nrow=2)))%*%b
array(c,dim=c(2,2,2))



p1=1/3
p2=1/3
p3=1/3
p=c(p1,p2,p3)
a<-expand.grid(0:1,0:1,0:1)
b=rep(0,2^3)
b[which(apply(a,1,sum)==0)]=1
b[which(apply(a,1,sum)==1)]=0
b[which(apply(a,1,sum)==2)]=p1*(1-p1)^2/(2-3*p1)
b[which(apply(a,1,sum)==3)]=0
c<-kronecker(matrix(c(1-p3,p3,-1,1),nrow=2),kronecker(matrix(c(1-p2,p2,-1,1),nrow=2),matrix(c(1-p1,p1,-1,1),nrow=2)))%*%b
array(c,dim=c(2,2,2))


p1=0.45
p2=0.45
p3=0.45
p=c(p1,p2,p3)
a<-expand.grid(0:1,0:1,0:1)
b=rep(0,2^3)
b[which(apply(a,1,sum)==0)]=1
b[which(apply(a,1,sum)==1)]=0
b[which(apply(a,1,sum)==2)]=p1*(1-p1)^2/(2-3*p1)
b[which(apply(a,1,sum)==3)]=0
c<-kronecker(matrix(c(1-p3,p3,-1,1),nrow=2),kronecker(matrix(c(1-p2,p2,-1,1),nrow=2),matrix(c(1-p1,p1,-1,1),nrow=2)))%*%b
array(c,dim=c(2,2,2))

p1=0.3
p2=0.4
p3=0.5
p=c(p1,p2,p3)
a<-expand.grid(0:1,0:1,0:1)
b=rep(0,2^3)
b[which(apply(a,1,sum)==0)]=1
b[which(apply(a,1,sum)==1)]=0
b[which(apply(a,1,sum)==2)]=p1*(1-p2)*(1-p3)/(2-sum(p))
b[which(apply(a,1,sum)==3)]=0
c<-kronecker(matrix(c(1-p3,p3,-1,1),nrow=2),kronecker(matrix(c(1-p2,p2,-1,1),nrow=2),matrix(c(1-p1,p1,-1,1),nrow=2)))
c<-c%*%b
c<-array(c,dim=c(2,2,2))
apply(c,c(1,2),sum)[2,2]-p1*p2
apply(c,c(1,3),sum)[2,2]-p1*p3
apply(c,c(2,3),sum)[2,2]-p2*p3




p1=0.1
p2=0.2
p3=0.2
p=c(p1,p2,p3)
a<-expand.grid(0:1,0:1,0:1)
c<-kronecker(matrix(c(1-p3,p3,-1,1),nrow=2),kronecker(matrix(c(1-p2,p2,-1,1),nrow=2),matrix(c(1-p1,p1,-1,1),nrow=2)))
d<-apply(c,1,function(x) c(-x[1]/sum(x[which(apply(a,1,sum)==2)]),sum(x[which(apply(a,1,sum)==2)])))
c(max(c(0,d[1,d[2,]>0])),min(d[1,d[2,]<0]))

b=rep(0,2^3)
b[which(apply(a,1,sum)==0)]=1
b[which(apply(a,1,sum)==1)]=0
b[which(apply(a,1,sum)==2)]=min(d[1,d[2,]<0])
b[which(apply(a,1,sum)==3)]=0
c<-c%*%b
c<-array(c,dim=c(2,2,2))
round(c*1000)/1000
apply(c,c(1,2),sum)[2,2]-p1*p2
apply(c,c(1,3),sum)[2,2]-p1*p3
apply(c,c(2,3),sum)[2,2]-p2*p3








p1=0.2
p2=0.3
p3=0.4
p=c(0.6,0,0,0,0,0.1,0.2,0.1)
p=c(0.6,0,0,0,0.1,0,0.1,0.2)
a<-expand.grid(0:1,0:1,0:1)

c<-kronecker(matrix(c(1,-p3,1,1-p3),nrow=2),kronecker(matrix(c(1,-p2,1,1-p2),nrow=2),matrix(c(1,-p1,1,1-p1),nrow=2)))%*%p
c<-array(c,dim=c(2,2,2))
round(c*1000)/1000
b=as.vector(c)

b0=b
d<-kronecker(matrix(c(1-p3,p3,-1,1),nrow=2),kronecker(matrix(c(1-p2,p2,-1,1),nrow=2),matrix(c(1-p1,p1,-1,1),nrow=2)))%*%b0
round(array(d,dim=c(2,2,2))*1000)/1000








p1=0.1
p2=0.2
p3=0.2
p4=0.3
p=c(p1,p2,p3)
a<-expand.grid(0:1,0:1,0:1,0:1)
c<-kronecker(matrix(c(1-p4,p4,-1,1),nrow=2),
             kronecker(matrix(c(1-p3,p3,-1,1),nrow=2),
                       kronecker(matrix(c(1-p2,p2,-1,1),nrow=2),matrix(c(1-p1,p1,-1,1),nrow=2))))
d<-apply(c,1,function(x) c(-x[1]/sum(x[which(apply(a,1,sum)==2)]),sum(x[which(apply(a,1,sum)==2)])))
c(max(c(0,d[1,d[2,]>0])),min(d[1,d[2,]<0]))

e=c[2,]
sum(e*b)

b=rep(0,2^4)
b[which(apply(a,1,sum)==0)]=1
b[which(apply(a,1,sum)==1)]=0
b[which(apply(a,1,sum)==2)]=min(d[1,d[2,]<0])
b[which(apply(a,1,sum)==3)]=0
c<-c%*%b
c<-array(c,dim=c(2,2,2,2))
round(c*1000)/1000
apply(c,c(1,2),sum)[2,2]-p1*p2
apply(c,c(1,3),sum)[2,2]-p1*p3
apply(c,c(2,3),sum)[2,2]-p2*p3





p1=0.1
p2=0.2
p3=0.2
p4=0.3
p5=0.3
p=c(p1,p2,p3,p4,p5)
a<-expand.grid(0:1,0:1,0:1,0:1,0:1)
a2<-a[which(apply(a,1,sum)==2),]
c<-kronecker(matrix(c(1-p5,p5,-1,1),nrow=2),kronecker(matrix(c(1-p4,p4,-1,1),nrow=2),
             kronecker(matrix(c(1-p3,p3,-1,1),nrow=2),
                       kronecker(matrix(c(1-p2,p2,-1,1),nrow=2),matrix(c(1-p1,p1,-1,1),nrow=2)))))
d<-apply(c,1,function(x) c(-x[1]/sum(x[which(apply(a,1,sum)==2)]*apply(a2,1,function(x) min(p[which(x==1)])*(1-max(p[which(x==1)])))),sum(x[which(apply(a,1,sum)==2)])))
c(max(c(0,d[1,d[2,]>0])),min(d[1,d[2,]<0]))

b=rep(0,2^5)
b[which(apply(a,1,sum)==0)]=1
b[which(apply(a,1,sum)==1)]=0
b[which(apply(a,1,sum)==2)]=min(d[1,d[2,]<0])*apply(a2,1,function(x) min(p[which(x==1)])*(1-max(p[which(x==1)])))
b[which(apply(a,1,sum)==3)]=0
c<-c%*%b
c<-array(c,dim=c(2,2,2,2,2))
round(c*1000)/1000
apply(c,c(1,2),sum)[2,2]-p1*p2
apply(c,c(1,3),sum)[2,2]-p1*p3
apply(c,c(2,3),sum)[2,2]-p2*p3


p1=0.1
p2=0.2
p3=0.2
p4=0.3
p5=0.3
p=c(p1,p2,p3,p4,p5)
a<-expand.grid(0:1,0:1,0:1,0:1,0:1)
a2<-a[which(apply(a,1,sum)==2),]
c<-kronecker(matrix(c(1-p5,p5,-1,1),nrow=2),kronecker(matrix(c(1-p4,p4,-1,1),nrow=2),
                                                      kronecker(matrix(c(1-p3,p3,-1,1),nrow=2),
                                                                kronecker(matrix(c(1-p2,p2,-1,1),nrow=2),matrix(c(1-p1,p1,-1,1),nrow=2)))))
d<-apply(c,1,function(x) c(-x[1]/sum(x[which(apply(a,1,sum)==2)]),sum(x[which(apply(a,1,sum)==2)])))
c(max(c(0,d[1,d[2,]>0])),min(d[1,d[2,]<0]))

b=rep(0,2^5)
b[which(apply(a,1,sum)==0)]=1
b[which(apply(a,1,sum)==1)]=0
b[which(apply(a,1,sum)==2)]=min(d[1,d[2,]<0])
b[which(apply(a,1,sum)==3)]=0
c<-c%*%b
c<-array(c,dim=c(2,2,2,2,2))
round(c*1000)/1000
apply(c,c(1,2),sum)[2,2]-p1*p2
apply(c,c(1,3),sum)[2,2]-p1*p3
apply(c,c(2,3),sum)[2,2]-p2*p3








p1=0.2
p2=0.2
p3=0.2
p4=0.3
p5=0.3
p6=0.3
p=c(p1,p2,p3,p4,p5,p6)
a<-expand.grid(0:1,0:1,0:1,0:1,0:1,0:1)
a2<-a[which(apply(a,1,sum)==2),]
c<-kronecker(matrix(c(1-p6,p6,-1,1),nrow=2),kronecker(matrix(c(1-p5,p5,-1,1),nrow=2),kronecker(matrix(c(1-p4,p4,-1,1),nrow=2),
                                                      kronecker(matrix(c(1-p3,p3,-1,1),nrow=2),
                                                                kronecker(matrix(c(1-p2,p2,-1,1),nrow=2),matrix(c(1-p1,p1,-1,1),nrow=2))))))
d<-apply(c,1,function(x) c(-x[1]/sum(x[which(apply(a,1,sum)==2)]*apply(a2,1,function(x) min(p[which(x==1)])*(1-max(p[which(x==1)])))),sum(x[which(apply(a,1,sum)==2)])))
c(max(c(0,d[1,d[2,]>0])),min(d[1,d[2,]<0]))

b=rep(0,2^6)
b[which(apply(a,1,sum)==0)]=1
b[which(apply(a,1,sum)==1)]=0
b[which(apply(a,1,sum)==2)]=min(d[1,d[2,]<0])*apply(a2,1,function(x) min(p[which(x==1)])*(1-max(p[which(x==1)])))
b[which(apply(a,1,sum)==3)]=0
c<-c%*%b
c<-array(c,dim=c(2,2,2,2,2,2))
round(c*1000)/1000
apply(c,c(1,2),sum)[2,2]-p1*p2
apply(c,c(1,3),sum)[2,2]-p1*p3
apply(c,c(2,3),sum)[2,2]-p2*p3






p1=0.1
p2=0.1
p3=0.1
p4=0.1
p5=0.1
p6=0.1
p=c(p1,p2,p3,p4,p5,p6)
a<-expand.grid(0:1,0:1,0:1,0:1,0:1,0:1)
a2<-a[which(apply(a,1,sum)==2),]
c<-kronecker(matrix(c(1-p6,p6,-1,1),nrow=2),kronecker(matrix(c(1-p5,p5,-1,1),nrow=2),kronecker(matrix(c(1-p4,p4,-1,1),nrow=2),
                                                                                               kronecker(matrix(c(1-p3,p3,-1,1),nrow=2),
                                                                                                         kronecker(matrix(c(1-p2,p2,-1,1),nrow=2),matrix(c(1-p1,p1,-1,1),nrow=2))))))
d<-apply(c,1,function(x) c(-x[1]/sum(x[which(apply(a,1,sum)==2)]*apply(a2,1,function(x) min(p[which(x==1)])*(1-max(p[which(x==1)])))),sum(x[which(apply(a,1,sum)==2)])))
c(max(c(0,d[1,d[2,]>0])),min(d[1,d[2,]<0]))

b=rep(0,2^6)
b[which(apply(a,1,sum)==0)]=1
b[which(apply(a,1,sum)==1)]=0
b[which(apply(a,1,sum)==2)]=0
b[which(apply(a,1,sum)==3)]=0
c<-c%*%b
c<-array(c,dim=c(2,2,2,2,2,2))
round(c*1000)/1000
apply(c,c(1,2),sum)[2,2]-p1*p2
apply(c,c(1,3),sum)[2,2]-p1*p3
apply(c,c(2,3),sum)[2,2]-p2*p3


p1=0.1
p2=0.2
p3=0.2
p=c(p1,p2,p3)
p=c()
a<-expand.grid(0:1,0:1,0:1)
c<-kronecker(matrix(c(1-p3,p3,-1,1),nrow=2),kronecker(matrix(c(1-p2,p2,-1,1),nrow=2),matrix(c(1-p1,p1,-1,1),nrow=2)))
d<-apply(c,1,function(x) c(-x[1]/sum(x[which(apply(a,1,sum)==2)]),sum(x[which(apply(a,1,sum)==2)])))
c(max(c(0,d[1,d[2,]>0])),min(d[1,d[2,]<0]))




p1=0.2
p2=0.3
p3=0.4
p=c(0.6,0,0,0,0,0.1,0.2,0.1)
#p=c(0.6,0,0,0,0.1,0,0.1,0.2)
a<-expand.grid(0:1,0:1,0:1)

c<-kronecker(matrix(c(1,-p3,1,1-p3),nrow=2),kronecker(matrix(c(1,-p2,1,1-p2),nrow=2),matrix(c(1,-p1,1,1-p1),nrow=2)))%*%p
c<-array(c,dim=c(2,2,2))
round(c*1000)/1000
b=as.vector(c)

b0=b
d<-kronecker(matrix(c(1-p3,p3,-1,1),nrow=2),kronecker(matrix(c(1-p2,p2,-1,1),nrow=2),matrix(c(1-p1,p1,-1,1),nrow=2)))%*%b0
round(array(d,dim=c(2,2,2))*1000)/1000

b0=b
b0[-1]=0.5*b0[-1]
d<-kronecker(matrix(c(1-p3,p3,-1,1),nrow=2),kronecker(matrix(c(1-p2,p2,-1,1),nrow=2),matrix(c(1-p1,p1,-1,1),nrow=2)))%*%b0
round(array(d,dim=c(2,2,2))*1000)/1000

b0=b
b0[-1]=0
d<-kronecker(matrix(c(1-p3,p3,-1,1),nrow=2),kronecker(matrix(c(1-p2,p2,-1,1),nrow=2),matrix(c(1-p1,p1,-1,1),nrow=2)))%*%b0
round(array(d,dim=c(2,2,2))*1000)/1000


p1=0.2
p2=0.3
p3=0.4
p=c(0.6,0,0,0,0,0.1,0.2,0.1)
#p=c(0.6,0,0,0,0.1,0,0.1,0.2)
a<-expand.grid(0:1,0:1,0:1)
c<-kronecker(matrix(c(1,-p3,1,1-p3),nrow=2),kronecker(matrix(c(1,-p2,1,1-p2),nrow=2),matrix(c(1,-p1,1,1-p1),nrow=2)))%*%p
b=as.vector(c)
res=c()
for(i in 1:101){
  b0=b
  b0[-1]=b0[-1]*(i-1)/100
  d<-kronecker(matrix(c(1-p3,p3,-1,1),nrow=2),kronecker(matrix(c(1-p2,p2,-1,1),nrow=2),matrix(c(1-p1,p1,-1,1),nrow=2)))%*%b0
  res[i]=d[1]
}
plot(1:101,res)



c<-kronecker(matrix(c(1-p3,p3,-1,1),nrow=2),kronecker(matrix(c(1-p2,p2,-1,1),nrow=2),matrix(c(1-p1,p1,-1,1),nrow=2)))%*%b
round(c*1000)/1000








p1=0.2
p2=0.3
p3=0.4
p=c(p1,p2,p3)
a<-expand.grid(0:1,0:1,0:1)
b=rep(0,2^3)
b[which(apply(a,1,sum)==0)]=1
b[which(apply(a,1,sum)==1)]=0
#b[which(apply(a,1,sum)==2)]=1
b[which(apply(a,1,sum)==2)]=c(p1*(1-p2),p1*(1-p3),p2*(1-p3))
b[which(apply(a,1,sum)==3)]=0
c<-kronecker(matrix(c(1-p3,p3,-1,1),nrow=2),kronecker(matrix(c(1-p2,p2,-1,1),nrow=2),matrix(c(1-p1,p1,-1,1),nrow=2)))%*%b
c<-array(c,dim=c(2,2,2))












apply(c,c(1,2),sum)[2,2]-p1*p2
apply(c,c(1,3),sum)[2,2]-p1*p3
apply(c,c(2,3),sum)[2,2]-p2*p3

plot(seq(0,1,length.out=100),(1-seq(0,1,length.out=100))/(0.5-seq(0,1,length.out=100)))




p1=1/3
p2=1/3
p3=1/3
p=c(p1,p2,p3)
a<-expand.grid(0:1,0:1,0:1)
b=rep(0,2^3)
c<-kronecker(matrix(c(1,-p3,1,1-p3),nrow=2),kronecker(matrix(c(1,-p2,1,1-p2),nrow=2),matrix(c(1,-p1,1,1-p1),nrow=2)))%*%b
array(c,dim=c(2,2,2))


























