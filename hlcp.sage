#!/usr/bin/python
load("hssp.sage")


class hlcp:
  def __init__(self,n,B,kappa):
    self.n=n
    self.kappa=kappa
    self.B=B
    
  def gen_instance(self,m=0):
    if m<H.n:
      m=choose_m(H.n,H.B)  
    self.m=m

    print "n=",self.n,"m=",m,
    if self.kappa>-1: print "kappa=",self.kappa,
    iota=0.035
    self.nx0=int(2*iota*H.n^2+H.n*log(H.n,2)+2*H.n*log(H.B,2))
    print "nx0=",self.nx0
    self.x0,self.a,self.X,self.b=genParams(self.n, self.m,self.nx0,self.kappa,self.B)
    
    
def hlcp_attack(H,alg='default'):
  
  if H.B==1: return  hssp_attack(H,alg)
  
  print "Random instance of HLCP_(n,B)^kappa(m,Q)"
  print "n=",H.n, "B=",H.B, "kappa=",H.kappa,"m=",H.m, "log(Q)=", H.nx0 
  print
  kappa=H.kappa
  n=H.n
  B=H.B
  
        
  if alg in ['default','staistical']:
      MO,tt1,tt10,tt1O= Step1(H.n,H.kappa,H.x0,H.a,H.X,H.b,H.m,BKZ=True)      
      tica, tt2, nrafound=statistical(MO,H.n,H.m,H.x0,H.X,H.a,H.b,H.kappa,H.B)
      return  tt1, tica, tt2, nrafound, H
      
  if alg=='ns_original' or (alg=='ns' and H.m==2*n):
      print "Nguyen-Stern (Original) Attack"
      MO,tt1,tt10,tt1O= Step1_original(H.n,H.kappa,H.x0,H.a,H.X,H.b,H.m)
      beta,tt2, nrafound=ns(H,MO,B)
      return tt1,tt10,tt1O,beta,tt2, nrafound, H       
      
  if alg=='ns':
      MO,tt1,tt10,tt1O= Step1(H.n,H.kappa,H.x0,H.a,H.X,H.b,H.m,BKZ=True)
      beta,tt2,nrafound=ns(H,MO,B)
      return tt1,tt10,tt1O,beta,tt2,nrafound, H

  return  
