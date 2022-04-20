
def gen_list(n, v):
  l=0
  L=[]
  while l<v:
      p=ZZ.random_element(n)
      if p not in L: 
        L+=[p]
        l+=1
  return L
  
  
def genpseudoprime(eta,etamin=211):
  if eta<=(2*etamin):
    return random_prime(2^eta,False,2^(eta-1))
  else:
    return random_prime(2^etamin,False,2^(etamin-1))*genpseudoprime(eta-etamin)

  
def orthoLatticeMod(b,n,x0):
  m=b.length()
  assert m>=3*n
  assert m % n==0
  M=Matrix(ZZ,m,3*n)
  M[:2*n,:2*n]=identity_matrix(2*n)
  for i in range(2,m/n):
    M[i*n:(i+1)*n,2*n:3*n]=identity_matrix(n)
  
  M[1:,0]=-b[1:]*inverse_mod(b[0],x0)
  M[0,0]=x0

  for i in range(1,m):
    M[i,0]=mod(M[i,0],x0)
  return M  
  
 

def choose_m(n,B):
    if n<=150:
    	m=n^2*(int(log(B,2))+1)
    else: 
    	m=n^2*B
    return m
  
def genParams(n=10,m=20,nx0=100,v=5,B=1):
  assert v<=n, "v>n!"
  #print "Generation of x0",
  x0=genpseudoprime(nx0)


  # We generate the alpha_i's
  a=vector(ZZ,n)
  for i in xrange(n):
    a[i]=mod(ZZ.random_element(x0),x0)

  # The matrix X has m rows and must be of rank n
  if v>0:
    while True:
      X=Matrix(ZZ,m,n)
      for i in xrange(m):
        Lin=gen_list(n, v) 
        for j in Lin:
          X[i,j]=ZZ.random_element(1,B+1)
        del Lin
      print X.rank()
      if X.rank()==n: break
  else:
    while True:
      X=Matrix(ZZ,m,n)
      for i in xrange(m):
        for j in xrange(n):
          X[i,j]=ZZ.random_element(B+1)
      print X.rank()    
      if X.rank()==n: break
  print X.density().n()  

  # We generate an instance of the HSSP: b=X*a
  b=X*a
  for i in range(m):
    b[i]=mod(b[i],x0)

  return x0,a,X,b
  
