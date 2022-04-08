#!/usr/bin/python
load("building.sage")
load("ortho_attack.sage")
load("multi.sage")
load("ns.sage")
load("statistical.sage")

class hssp:
  def __init__(self,n,kappa):
    self.n=n
    self.kappa=kappa
    
    
  def gen_instance(self,m=0):
    if m==0 and self.n % 2==1:
      m=self.n*(self.n+3)/2 # n is odd
    elif m==0 and self.n % 2 ==0:
      m=self.n*(self.n+4)/2 # n is even
    self.m=m
    
    print "n=",self.n,"m=",m,
    if self.kappa>-1: print "kappa=",self.kappa,
    iota=0.035
    self.nx0=int(2*iota*self.n^2+self.n*log(self.n,2))
    print "nx0=",self.nx0
    self.x0,self.a,self.X,self.b=genParams(self.n, self.m,self.nx0,self.kappa)
    
    
def hssp_attack(H,alg='default'):
  kappa=H.kappa
  n=H.n
  
  if kappa!=-1:
    print "Random instance of HSSP_n^kappa(m,Q)"
    print "n=",H.n, "kappa=",H.kappa,"m=",H.m, "log(Q)=", H.nx0 
    print
  else:
    print "Random instance of HSSP_n(m,Q)"
    print "n=",H.n,"m=",H.m, "log(Q)=", H.nx0 
    print
  if alg in ['default','multi']:
    assert H.m>(n^2+n)/2, 'm too small'
    MO,tt1= Step1(H.n,H.kappa,H.x0,H.a,H.X,H.b,H.m,BKZ=False)
    print "\nMultivariate Attack"
    if kappa>0:
      tei, tef, tsf,tt2, nrafound=bit_guessing(H.n,H.kappa,MO,H.x0,H.a,H.X,H.b,H.m) 
      return  tt1, tei, tef, tsf, tt2,nrafound, H
    else:
      tei, tef, tt2,nrafound=eigen(H.n,H.kappa,MO,H.x0,H.a,H.X,H.b,H.m) 
      return  tt1, tei, tef, tt2,nrafound, H
  
  if alg=='ns_original' or (alg=='ns' and H.m==2*n):
      print "Nguyen-Stern (Original) Attack"
      MO,tt1,tt10,tt1O= Step1_original(H.n,H.kappa,H.x0,H.a,H.X,H.b,H.m)
      beta,tt2, nrafound=ns(H,MO)
      return tt1,tt10,tt1O,beta,tt2, nrafound, H   
      
  if alg=='ns':
      MO,tt1,tt10,tt1O= Step1(H.n,H.kappa,H.x0,H.a,H.X,H.b,H.m,BKZ=True)
      print "\nNguyen-Stern (Improved) Attack"
      beta,tt2, nrafound=ns(H,MO)
      return tt1,tt10,tt1O,beta,tt2, nrafound, H
     
      
  if alg=='statistical':
    assert kappa==-1, 'The statistical attack does not work for HSSP_n^kappa(m,Q)'
    MO,tt1,tt10,tt1O= Step1(H.n,H.kappa,H.x0,H.a,H.X,H.b,H.m,BKZ=True)
    print "\nStatistical Attack"      
    tica, tt2, nrafound=statistical(MO,H.n,H.m,H.x0,H.X,H.a,H.b,H.kappa)
    return  tt1, tica, tt2, nrafound, H
  
  return None
