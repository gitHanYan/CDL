############ Codes for generating simulated data 
N=c(100,200,300)
m=1
n=1
p=N[m]# the dimension of variables
q=N[n]# the number of observations
b=p/50# number of blocks
inv_A=matrix(0,p,p)#precision matrix of A
inv_B=matrix(0,p,p)#precision matrix of B

#####setting1(Auto-regression)
A=matrix(0,50,50)#blocks of invA
B=matrix(0,50,50)#blocks of invB
for (i in 1:50) {
  for (j in 1:50) {
    A[i,j]=0.6^(abs(i-j))
  }
}
B=A
for (i in 1:50) {
  for (j in 1:50) {
    if(abs(i-j)==25){
      B[i,j]=0.15
    }
  }
}
for (k in 1:b) {
  inv_A[(50*k-49):(50*k),(50*k-49):(50*k)]=A
  inv_B[(50*k-49):(50*k),(50*k-49):(50*k)]=B
}
#####setting1----End
#####setting2
p1=p^2
data=1:p1
p2=sample(data,size=p1,replace=F)
p3=runif(p2,min = -0.09,max = 0.09)
A_mat=matrix(p3,p,p)
B_mat=matrix(0,p,p)
diag_element=c(1,1.2,1.4)# make sure the matrix is positive definite
for (i in 1:p) {
  B_mat[i,]=rbinom(p,1,0.01)
}
inv_A=A_mat*B_mat
inv_A[!upper.tri(inv_A,diag=TRUE)]=0
inv_A=inv_A+t(inv_A)

for (k in 1:b) {
  for (i in 1:50) {
    A[i,1:50]=runif(50,-0.1,0.1)
  }
  for (i in 1:50) {
    B[i,]=rbinom(50,1,0.4)#40% elements are nonzero
  }
  A=B*A
  A[!upper.tri(A,diag=TRUE)]=0
  A=A+t(A)
  for (i in 1:50) {
    A[i,i]=diag_element[m]
  }
  for (i in 1:50) {
    C[i,]=runif(50,-0.5,0.5)
  }
  for (i in 1:50) {
    D[i,]=rbinom(50,1,0.02)#2% elements are different
  }
  C=C*D
  C[!upper.tri(C,diag=TRUE)]=0
  C=C+t(C)
  B=C+A
  inv_A[(50*k-49):(50*k),(50*k-49):(50*k)]=A
  inv_B[(50*k-49):(50*k),(50*k-49):(50*k)]=B
}
#####setting2----End
#####setting3(Nearest-neighboorhood)
umin=0.4
umax=0.8
p2=50
inv_A=diag(p)
near=5
b=2
for (k in 1:b) {
  subi = diag(1,p2,p2)
  point.matrix = matrix(runif(2*p2,0,1),p2,2)
  num.nearest = matrix(0,p2,near)
  for (j in 1:p2) {
    distancej = apply((t(point.matrix[-j,]) - point.matrix[j,])^2,2,sum)
    corj = setdiff(1:p,j)[order(distancej)[1:near]]
    num.nearest[j,] = corj
  }
  for (oi in 1:dim(num.nearest)[1]) {
    for (oj in 1:dim(num.nearest)[2]) {
      i=oi;j=num.nearest[oi,oj]
      ij = sample(c(runif(1,umin,umax),runif(1,-umax,-umin)))[1]
      subi[i,j]=ij;subi[j,i]=ij
    }
  }
  for (i in 1:p2) {
    subi[i,i] = sum(abs(subi[i,setdiff(1:p2,i)]))+0.3
  }
  inv_A[(50*k-49):(50*k),(50*k-49):(50*k)]=subi
}
C=matrix(0,50,50)
D=matrix(0,50,50)
Delta=matrix(0,p,p)
for (k in 1:b) {
  for (i in 1:50) {
    C[i,]=runif(50,-0.5,0.5)
  }
  for (i in 1:50) {
    D[i,]=rbinom(50,1,0.02)#2% elements are different
  }
  C=C*D
  C[!upper.tri(C,diag=TRUE)]=0
  C=C+t(C)
  Delta[(50*k-49):(50*k),(50*k-49):(50*k)]=C
}
inv_B=inv_A+Delta
#####setting3----End

#####generate 0~5 % percent elements outsides the blocks
p1=p^2
data=1:p1
p2=sample(data,size=p1,replace=F)
p3=runif(p2,min = -0.09,max = 0.09)
A_mat=matrix(p3,p,p)
B_mat=matrix(0,p,p)

for (i in 1:p) {
  B_mat[i,]=rbinom(p,1,0.05)####determine the percentage of nonzero elements
}
error=A_mat*B_mat
error[!upper.tri(error,diag=TRUE)]=0
error=error+t(error)
for (k in 1:b) {
  error[(50*k-49):(50*k),(50*k-49):(50*k)]=0
}
inv_A=inv_A+error
