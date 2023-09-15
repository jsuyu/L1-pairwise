library("ADMM")      #metodologia per determinar optimitzar Lasso
library("coda.base") #paquet de coordenades composicionals
library("dplyr")

#---------------------------------------------------------------------------------------
set.seed(333)
#carregem la matriu de dades
GEMAS<-read.csv("GEMAS.csv",sep=',')
ilr_gemas<-read.csv("ilr_gemas.csv",sep=';')
sel<-c("SiO2","TiO2","Al2O3","Fe2O3","MnO","MgO","CaO","Na2O","K2O","P2O5","LOI")
ilr<-c("ilr1","ilr2","ilr3","ilr4","ilr5","ilr6","ilr7","ilr8","ilr9","ilr10")
x<-GEMAS
#zero replacement
#x[x==0]<-1
# nomes algunes de les columnes que usen en el paper
x<-x[,c("sand","silt","clay","AnnPrec","pH_CaCl2","PM","SiO2","TiO2","Al2O3", "Fe2O3",
        "MnO","MgO","CaO","Na2O","K2O","P2O5","LOI")]
# mirem quantes i quines files tenen NA
which(apply(is.na(x),1,sum)>0)#191 308 445 568 592 937
# mirem quantes i quines files tenen zeros
which(apply((x==0),1,sum)>0)#138 1100 1218 1548 1634 1988
#
### ELIMINEM files amb NA i ZEROS
#carreguem la variable resposta
y<-as.matrix(GEMAS[,"pH_CaCl2"])
elim<-!((apply(is.na(x),1,sum)>0)|(apply((x==0),1,sum)>0))
y<-y[elim]
x<-x[elim,]

x<-as.matrix(x[,sel])
#---------------------------------------------------------------------------------------
#passem a coordenades olr la matriu de dades
n<-nrow(x)                            #n=nombre de individus de la composiciÃ³
D<-ncol(x)                            #D=nombre de parts de la composiciÃ³
#tenim tres bases olr
psi<-sbp_basis(as.matrix(ilr_gemas[,ilr]))
psi<-ilr_basis(D)
psi<-sbp_basis(matrix(c(1,-1,0,0,0,0,0,0,0,0,0,
                        1,1,-1,0,0,0,0,0,0,0,0,
                        1,1,1,-1,0,0,0,0,0,0,0,
                        0,0,0,0,1,0,0,0,-1,0,0,
                        -1,-1,-1,-1,1,0,0,0,1,0,0,
                        -1,-1,-1,-1,-1,1,1,1,-1,1,1,
                        0,0,0,0,0,0,1,-1,0,0,0,
                        0,0,0,0,0,1,0,0,0,-1,0,
                        0,0,0,0,0,1,0,0,0,1,-1,
                        0,0,0,0,0,-1,1,1,0,-1,-1),D,D-1))
olr_x<-coordinates(x,psi)

#---------------------------------------------------------------------------------------
#constuion la matriu G que ens defineix la mÃ¨trica L1-pairwise
G<-matrix(0,(D-1)*(D)/2,D)            #inicialitzem la matiu G
colini<-as.vector(1:(D-1))            #columna on comenÃ§a cada submatriu
dimsub<-as.vector((D-1):1)            #dimensiÃ³ de cada submatriu
rowini<-as.vector(rep(1,D-1))         #fila on comenÃ§a cada submatriu
for (i in 2:(D-1)){
  rowini[i]=rowini[i-1]+dimsub[i-1]
}
#emplenem la matiu G
for(i in 1:(D-1)){
  G[rowini[i]:(rowini[i]+dimsub[i]-1),(colini[i]+1):D]<--diag(dimsub[i])
  G[rowini[i]:(rowini[i]+dimsub[i]-1),colini[i]]<-1
}
G<-1/(D-1)*G
Gpsi<-G%*%psi                       #matriu que cal aplicar en coordenades olr per tenir la mÃ¨trica L1-pairwise

#---------------------------------------------------------------------------------------
# solve LASSO via reducing from Generalized LASSO
# Diagonal = diag(D-1) # matriu Identitat si volem fer LASSO
A = olr_x
A<-cbind(rep(1,n),A) #afegim el terme independent del model lineal
psi<-cbind(rep(0,D),psi) #eliminem el terme independent del model lineal
Gpsi<-cbind(rep(0,D*(D-1)/2),Gpsi) #eliminem el terme independent del model lineal
theta<-0.5
enetGpsi<-rbind(theta*Gpsi,(1-theta)*psi)
b = y
llista=list() #llista on guardarem les lambdes i parÃ metres beta

t1<-proc.time()
n.lambda<-100
vector.lambda<-matrix(c(seq(-2,8,length.out=n.lambda)))
exp(vector.lambda) #visualitzem la particiÃ³ de les lambdes
#definim el nombre de folders
n.fold<-10
MSE<-matrix(0,n.lambda,n.fold) #inicialitzem la matiu on guardarem els MSE

for(j in 1:n.lambda){ #fem correr les lambda
  yourdata<-cbind(b,A)
  #Randomly shuffle the data
  yourdata<-yourdata[sample(nrow(yourdata)),]
  #Create 10 equally size folds
  folds <- cut(seq(1,nrow(yourdata)),breaks=n.fold,labels=FALSE)
  #Perform 10 fold cross validation
  for(i in 1:n.fold){ #fem correr els folders
    #Segement your data by fold using the which() function 
    testIndexes <- which(folds==i,arr.ind=TRUE)
    testData <- yourdata[testIndexes, ]
    trainData <- yourdata[-testIndexes, ]
    A.Test<-testData[,-1,drop=FALSE]
    b.Test<-testData[,1,drop=FALSE]
    A.Train<-trainData[,-1,drop=FALSE]
    b.Train<-trainData[,1,drop=FALSE]
    #Use the test and train data partitions however you desire...
    output  = admm.genlasso(A.Train,b.Train,Gpsi,lambda=exp(vector.lambda[j]), 
                            rho=1.0, 
                            alpha=1.0,
                            abstol=1e-4, 
                            reltol=1e-4, 
                            maxiter=5000000)
    param<-output$x
    MSE[j,i]<-mean((b.Test-A.Test%*%param)^2)
  }
}
MSE<-cbind(MSE,rowMeans(MSE),apply(MSE[,1:n.fold],1, sd)/sqrt(n.fold))
MSE<-cbind(MSE,MSE[,n.fold+1]-MSE[,n.fold+2],MSE[,n.fold+1]+MSE[,n.fold+2],vector.lambda,exp(vector.lambda))
t2<-proc.time()

i.lambda.min<-which.min(MSE[,n.fold+1])
log.lambda.min<-MSE[i.lambda.min,n.fold+5]
se.lambda.min<-MSE[i.lambda.min,n.fold+2]
i.lambda.se<-max(which(sapply(MSE[,n.fold+1], function(x)x<MSE[i.lambda.min,n.fold+1]+se.lambda.min&x>MSE[i.lambda.min,n.fold+1]-se.lambda.min)))
log.lambda.se<-MSE[i.lambda.se,n.fold+5]

plot(MSE[,n.fold+5],MSE[,n.fold+1],type="p",col=1,xlab=expression(ln(lambda)),ylab="MSE")
lines(MSE[,n.fold+5],MSE[,n.fold+3],type="l",col=2)
lines(MSE[,n.fold+5],MSE[,n.fold+4],type="l",col=2)
abline(v=log.lambda.min)
abline(v=log.lambda.se)
