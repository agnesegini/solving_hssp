# Computes the right kernel of M using LLL.
# We assume that m>=2*n. This is only to take K proportional to M.height()
# We follow the approach from https://hal.archives-ouvertes.fr/hal-01921335/document
def kernelLLL(M):
  n=M.nrows()
  m=M.ncols()
  if m<2*n: return M.right_kernel().matrix()
  K=2^(m//2)*M.height()
  
  MB=Matrix(ZZ,m+n,m)
  MB[:n]=K*M
  MB[n:]=identity_matrix(m)
  
  MB2=MB.T.LLL().T
  
  assert MB2[:n,:m-n]==0
  Ke=MB2[n:,:m-n].T

  return Ke
  

def matNbits(M):
  return max([M[i,j].nbits() for i in range(M.nrows()) for j in range(M.ncols())])

def NZeroVectors(M):
  return sum([vi==0 and 1 or 0 for vi in M])


# Matrix rounding to integers
def roundM(M):
  M2=Matrix(ZZ,M.nrows(),M.ncols())
  for i in range(M.nrows()):
    for j in range(M.ncols()):
      M2[i,j]=round(M[i,j])
  return M2
  


def Step1_plus_gen(n,v,m=0,BKZ=False):
  if m==0 and n % 2==1:
    m=n*(n+3)/2 # n is odd
  elif m==0 and n % 2 ==0:
    m=n*(n+4)/2 # n is even
  k=4

  print "n=",n,"m=",m,"k=",k,
  if v>-1: print "kappa=",v,

  iota=0.035
  nx0=int(2*iota*n^2+n*log(n,2))
  print "nx0=",nx0

  x0,a,X,b=genParams(n,m,nx0,v)

  M=orthoLatticeMod(b,n,x0)

  print "Step 1",
  t=cputime()

  M[:n//k,:n//k]=M[:n//k,:n//k].LLL()
  
  M2=M[:2*n,:2*n].LLL()
  tprecomp=cputime(t)
  print "  LLL:%.1f" % tprecomp,

  RF=RealField(matNbits(M))

  M4i=Matrix(RF,M[:n//k,:n//k]).inverse()
  M2i=Matrix(RDF,M2).inverse()  
    
  ts1=cputime()
  while True:
    flag=True
    for i in range((m/n-2)*k):
      indf=2*n+n//k*(i+1)
      if i==(m/n-2)*k-1:
        indf=m
        
      mv=roundM(M[2*n+n//k*i:indf,:n//k]*M4i)
      if mv==0: 
        continue
      flag=False
      M[2*n+n//k*i:indf,:]-=mv*M[:n//k,:]
    if flag: break
  print "  Sred1:%.1f" % cputime(ts1),

  M[:2*n,:2*n]=M2

  ts2=cputime()
  while True:
    #print "  matNBits(M)=",matNbits(M[2*n:])
    mv=roundM(M[2*n:,:2*n]*M2i)
    if mv==0: break
    M[2*n:,:]-=mv*M[:2*n,:]
  tt10=cputime(t)
  print "  Sred2:%.1f" % cputime(ts2),
  
  # The first n vectors of M should be orthogonal
  northo=NZeroVectors(M[:n,:2*n]*X[:2*n])

  for i in range(2,m/n):
    northo+=NZeroVectors(M[i*n:(i+1)*n,:2*n]*X[:2*n]+X[i*n:(i+1)*n])

  print "  #ortho vecs=",northo,"out of",m-n,
  
  # Orthogonal of the orthogonal vectors
  # We compute modulo 3 if multivariate
  if BKZ: KK=ZZ
  else: KK=GF(3)
  MO=Matrix(KK,n,m)
  
  tk=cputime()
  MO[:,:2*n]=kernelLLL(M[:n,:2*n])
  print "  Kernel LLL: %.1f" % cputime(tk),

  for i in range(2,m/n):
    MO[:,i*n:(i+1)*n]=-(M[i*n:(i+1)*n,:2*n]*MO[:,:2*n].T).T
  #print "Total kernel computation",cputime(tk)
  tt1=cputime(t)
  tt1O=cputime(tk)
  print "  Total Step 1: %.1f" % tt1
  #L=IntegerLattice(MO)
  #Y=X.T
  #R=Y.rows()
  #print all(w in L for w in R)
  #LX=IntegerLattice(Y)
  #print all(w in LX for w in MO.rows())
  
  del M,mv
  if BKZ:return MO,x0,a,X,b,int(m),tt1,tt10,tt1O
  else: return MO,x0,a,X,b,int(m),tt1
  

def Step1(n,v,x0,a,X,b,m,BKZ=False):
  k=4
  M=orthoLatticeMod(b,n,x0)
  print "Step 1",
  t=cputime()

  M[:n//k,:n//k]=M[:n//k,:n//k].LLL()
  
  M2=M[:2*n,:2*n].LLL()
  tprecomp=cputime(t)
  print "  LLL:%.1f" % tprecomp,

  RF=RealField(matNbits(M))

  M4i=Matrix(RF,M[:n//k,:n//k]).inverse()
  M2i=Matrix(RDF,M2).inverse()  
    
  ts1=cputime()
  while True:
    flag=True
    for i in range((m/n-2)*k):
      indf=2*n+n//k*(i+1)
      if i==(m/n-2)*k-1:
        indf=m
        
      mv=roundM(M[2*n+n//k*i:indf,:n//k]*M4i)
      if mv==0: 
        continue
      flag=False
      M[2*n+n//k*i:indf,:]-=mv*M[:n//k,:]
    if flag: break
  print "  Sred1:%.1f" % cputime(ts1),

  M[:2*n,:2*n]=M2

  ts2=cputime()
  while True:
    #print "  matNBits(M)=",matNbits(M[2*n:])
    mv=roundM(M[2*n:,:2*n]*M2i)
    if mv==0: break
    M[2*n:,:]-=mv*M[:2*n,:]
  tsr=cputime(ts2)
  tt10=cputime(t)
  print "  Sred2:%.1f" % cputime(ts2),
  
  # The first n vectors of M should be orthogonal
  northo=NZeroVectors(M[:n,:2*n]*X[:2*n])

  for i in range(2,m/n):
    northo+=NZeroVectors(M[i*n:(i+1)*n,:2*n]*X[:2*n]+X[i*n:(i+1)*n])

  print "  #ortho vecs=",northo,"out of",m-n,
  
  # Orthogonal of the orthogonal vectors
  # We compute modulo 3 if multivariate
  if BKZ: KK=ZZ
  else: KK=GF(3)
  MO=Matrix(KK,n,m)
  
  tk=cputime()
  MO[:,:2*n]=kernelLLL(M[:n,:2*n])
  print "  Kernel LLL: %.1f" % cputime(tk),

  for i in range(2,m/n):
    MO[:,i*n:(i+1)*n]=-(M[i*n:(i+1)*n,:2*n]*MO[:,:2*n].T).T
  #print "Total kernel computation",cputime(tk)
  tt1=cputime(t)
  tt1O=cputime(tk)
  print "  Total Step 1: %.1f" % tt1
  #L=IntegerLattice(MO)
  #Y=X.T
  #R=Y.rows()
  #print all(w in L for w in R)
  #LX=IntegerLattice(Y)
  #print all(w in LX for w in MO.rows())
  
  del M,mv
  return MO,tt1,tt10,tt1O
  #else: return MO, tt1
  


# We generate the lattice of vectors orthogonal to b modulo x0
# and also to c in the affine case
def orthoLattice(b,x0):
 m=b.length()
 M=Matrix(ZZ,m,m)

 for i in range(1,m):
      M[i,i]=1
 M[1:m,0]=-b[1:m]*inverse_mod(b[0],x0)
 M[0,0]=x0

 for i in range(1,m):
      M[i,0]=mod(M[i,0],x0)

 return M
  

def Step1_original(n,v,x0,a,X,b,m):
  
  M=orthoLattice(b,x0)

  print "Step 1",
  t=cputime()
  M2=M.LLL()
  tt10=cputime(t)
  print "LLL step1: %.1f" % cputime(t),

  assert sum([vi==0 and 1 or 0 for vi in M2*X])==m-n
  
  MOrtho=M2[:m-n]

  #print
  #for i in range(m-n+1):
  #  print i,N(log(M2[i:i+1].norm(),2)),N(log(m^(n/(2*(m-n)))*sqrt((m-n)/17),2)+iota*m+nx0/(m-n)) #N(log(sqrt((m-n)*n)*(m/2)^(m/(2*(m-n))),2)+iota*m)
  
  print "  log(Height,2)=",int(log(MOrtho.height(),2)),

  t2=cputime()
  ke=kernelLLL(MOrtho)
  tt1O=cputime(t2)
  print "  Kernel: %.1f" % cputime(t2),
  tt1=cputime(t)
  print "  Total step1: %.1f" % tt1

  return ke,tt1,tt10,tt1O


