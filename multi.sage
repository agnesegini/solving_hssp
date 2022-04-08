def ned(v):
     m=len(v)
     ones=vector([1]*m)
     return ones-v
     

def fullKernelMatrix(ke3,n,K):
  ke23=Matrix(K,ke3.nrows(),n*n)
  ind=n
  for i in range(n):
    for j in range(i,n):
      ke23[:,i*n+j]=ke3[:,ind]
      ke23[:,j*n+i]=ke3[:,ind]
      ind+=1
  return ke23   

def EigenspaceGuess(n,m,ke23,MO,kappa,verbose=True):
  li=[identity_matrix(ke23.nrows())]
 
  # We guess the j-th coordinate of the binary vector
  for j in range(n):
    ma=max(v.nrows() for v in li)
    if verbose: print "j=",j,"len(li)=",len(li),"tot=",sum(v.nrows() for v in li),"max=",ma
    if ma==1: break
    M=Matrix(GF(3),ke23.nrows(),n)
    for i in range(n):
      M[:,i]=ke23[:,i*n:(i+1)*n]*MO[:,j]
    l2=[]
    for V1 in li:
      if V1.nrows()==1: l2.append(V1)
    for e in range(2):
      if kappa>0 and j==0 and e==1: continue # only do this for constant Hamming weight
      Me=M[:,:]
      Me[:n,:]-=e*identity_matrix(n)
      V=Me.left_kernel().matrix()
      #print  V.dimensions()
      for V1 in li:
        if V1.nrows()>1:
          inter=spaceInter(V1,V)
          if inter.nrows()>0:
            l2.append(inter)
    li=l2
  li=[v[0][:n] for v in li]
  print "Depth: ", j,
  return li
  
  # We compute the intersection of the two spaces of row vectors
def spaceInter(V1,V2):
  VI=Matrix(GF(3),V1.nrows()+V2.nrows(),V1.ncols())
  VI[:V1.nrows(),:]=V1
  VI[V1.nrows():,:]=V2
  lk=VI.left_kernel().matrix()
  return lk[:,:V1.nrows()]*V1  
  
  
def bit_guessing(n,kappa,MO,x0,a,X,b,m):
  p=3
  K=GF(p)
  t2=cputime()
  xt23=Matrix(K,[(-x).list()+[x[i]*x[j]*((i==j) and 1 or 2) for i in range(n) for j in range(i,n)] for x in MO.T])

  ke3=xt23.right_kernel().matrix()
  ke23=fullKernelMatrix(ke3,n,K)
  print "dim(ker E)= ", ke3.dimensions()
  
  del xt23,ke3
  
  tei=cputime()
  li=EigenspaceGuess(n,m,ke23,MO, kappa,verbose=False)
  tef=cputime(tei)
  print "  len(li)=",len(li),"time eigen= %.1f"  %tef,
  NS=Matrix(li)*MO
  for i in range(NS.nrows()):
    if any(c==2 for c in NS[i]): NS[i]=-NS[i]

  print "  Number of recovered vectors:",NS.nrows()

  if kappa!=-1:
    ones=vector([1 for i in range(m)])
    li=NS.rows()
    for NSi in NS:
      if ones-NSi not in li:
        li.append(ones-NSi)
    NS=matrix(li)
  
  #nfound=len([True for NSi in NS if NSi in X.T])
  #print "  NFound=",nfound,"out of",n,
  ts=cputime()
  if kappa!=-1:
    e=vector([1]*n)
    if kappa<=n//2: Y=matrix(ZZ,NS[:n])
    else: Y=matrix(ZZ,NS[n:2*n])
    v0=variance(e*Y).n()
    i=0
    while True:
      j=i%n
      Y[j]=ned(Y[j])
      v1=variance(e*Y).n()
      #print v0, v1
      if v1==0: break
      if v1<v0: v0=v1
      else: Y[j]=ned(Y[j])
      i=i+1
    print "Switching rounds", i,"time= %.1f"  %cputime(ts) 
    #assert max(e*Y)==min(e*Y)==kappa, (max(e*Y),min(e*Y))
    tsf=cputime(ts)
    nyfound=len([True for y in Y if y in X.T])
    print "  NFound=",nyfound,"out of",n, 
    
    Y=Y.T
  else:
    Y=NS.T
    tsf=0
    nfound=len([True for NSi in NS if NSi in X.T])
    print "  NFound=",nfound,"out of",n,

  invNSn=matrix(Integers(x0),Y[:n]).inverse()
  ra=invNSn*b[:n]
  nrafound=len([True for rai in ra if rai in a])  
  
  print "  Total step2: %.1f" % cputime(t2),
  print "  Coefs of a found=",nrafound,"out of",n

  if kappa==n/2 and nyfound==0:
    print "\n-->Reverse"
    nyfound=len([True for y in Y.T if ones - y in X.T])
    print "  NFound=",nyfound,"out of",n,
    O=ones_matrix(ZZ,*Y.dimensions())
    Y=O-Y
    invNSn=matrix(Integers(x0),Y[:n]).inverse()
    ra=invNSn*b[:n]
    nrafound=len([True for rai in ra if rai in a])
    print "  Coefs of a found=",nrafound,"out of",n
  
  tt2=cputime(t2)
  
  return tei-t2, tef, tsf, tt2 ,nrafound
  
def eigen(n,kappa,MO,x0,a,X,b,m):
  print "Step 2",
  t2=cputime()
  K=GF(3)
  xt23=Matrix(K,[(-x).list()+[x[i]*x[j]*((i==j) and 1 or 2) for i in range(n) for j in range(i,n)] for x in MO.T])
  ke3=xt23.right_kernel().matrix()
  print "  Kernel: %.1f" % cputime(t2),

  assert xt23.nrows()==m
  assert xt23.ncols()==n*(n+1)/2+n

  ke23=fullKernelMatrix(ke3,n,K)
  tei=cputime()
  # We will compute the list of eigenvectors
  # We start with the full space.
  # We loop over the coordinates. This will split the eigenspaces.
  li=[Matrix(K,identity_matrix(n))]    
  for j in range(n):       # We loop over the coordinates of the wi vectors.
    #print "j=",j
    M=ke23[:,j*n:(j+1)*n]   # We select the submatrix corresponding to coordinate j
    li2=[]                 # We initialize the next list
    for v in li:
      if v.nrows()==1:     # We are done with this eigenvector 
        li2.append(v)
      else:     # eigenspace of dimension >1
        #print "eigenspace of dim:",v.nrows()
        A=v.solve_left(v*M)  # v*M=A*v. When we apply M on the right, this is equivalent to applying the matrix A.
                              # The eigenvalues of matrix A correspond to the jth coordinates of the wi vectors in that
                              # eigenspace
        for e,v2 in A.eigenspaces_left():    # We split the eigenspace according to the eigenvalues of A.
          vv2=v2.matrix()
          #print "  eigenspace of dim:",(vv2*v).nrows()
          li2.append(vv2*v)                   # The new eigenspaces 

    li=li2
  tef=cputime(tei)
  print "Eigenvectors computation", tef

  NS=Matrix([v[0] for v in li])*MO
  for i in range(n):
    if any(c==2 for c in NS[i]): NS[i]=-NS[i]

  print "  Number of recovered vectors:",NS.nrows(),

  nfound=len([True for NSi in NS if NSi in X.T])
  print "  NFound=",nfound,"out of",n,

  NS=NS.T
 
  # b=X*a=NS*ra
  invNSn=matrix(Integers(x0),NS[:n]).inverse()
  ra=invNSn*b[:n]
  nrafound=len([True for rai in ra if rai in a])
  
  tt2=cputime(t2)
  print "  Coefs of a found=",nrafound,"out of",n,
  print "  Total step2: %.1f" % tt2, 
  print

  return tei-t2, tef, tt2,nrafound 

  


