import numpy as np
from fica import ICA

def Sage2np(MO,n,m):
  MOn=np.matrix(MO)
  return MOn

def runICA(MOn,B=1):
  t1=cputime()
  A_,S=ICA(MOn,B)
  S_=np.dot(np.linalg.inv(A_),MOn)
  print "time Ica_S: ", cputime(t1),
  S2=matrix(ZZ,MOn.shape[0],MOn.shape[1], round(S_)) 
  return S2

def runICA_A(MOn,B=1,n=16,kappa=-1):
  t1=cputime()
  A_,S=ICA(MOn,B,n,kappa)
  print "time Ica_A: ", cputime(t1),
  A2=matrix(ZZ,MOn.shape[0],MOn.shape[0],round(A_))
  return A2


def statistical(MO,n,m,x0,X,a,b,kappa,B=1,variant=None):
  
  if variant==None:
    if n<=170: variant='roundA'
    else: variant='roundX'

  print "Step 2-ICA: ", variant
  
  t2=cputime()
  #print "matNbits=",matNbits(MO),
  tlll=cputime()
  MO=MO.LLL()
  print " time LLL=",cputime(tlll),"mathNbits=",matNbits(MO),

  MOn=Sage2np(MO,n,m)

  if variant=="roundA":
    A2=runICA_A(MOn,B,n,kappa)
    try:
      S2=A2.inverse()*MO
      print "mathNbits A=",matNbits(A2),
    except: 
      return 0,0
  elif variant=="roundX":
    S2=runICA(MOn,B)
    print "mathNbits X=",matNbits(X),
  else:
    raise NameError('Variant algorithm non acceptable')
  tica=cputime(t2)
  print " cputime ICA %.2f" %tica,

  tc=cputime()
  Y=X.T
  nfound=0
  for i in xrange(n):
    for j in xrange(n):
      if S2[i,:n]==Y[j,:n] and S2[i]==Y[j]:
        nfound+=1
  t=cputime(tc)      
  print "  NFound=",nfound,"out of",n,"check= %.2f" %t

  if nfound<n: 
    print
    print
    return tica,nfound

  NS=S2.T
  
  tcoff=cputime()
  # b=X*a=NS*ra
  invNSn=matrix(Integers(x0),NS[:n]).inverse()
  ra=invNSn*b[:n]
  tcf= cputime(tcoff)
  
  nrafound=len([True for rai in ra if rai in a])   
  print "  Coefs of a found=",nrafound,"out of",n, " time= %.2f" %tcf
  
  tS2=tcf+tica
  print "  Total step2: %.1f" % tS2,
  
  return tica,tS2, nfound

