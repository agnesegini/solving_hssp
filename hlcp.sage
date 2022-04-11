#!/usr/bin/python
"""
This is an implementation of the attacks for solving HLCP. The class hlcp can be used to generate  instances both of HLCP_{n,B} and HSSP_{n,B}^k, either specifing k or not, respectilvely.

The function hlcp_attack performs the attacks.

To run experiments with HLCP_{n,B} and HSSP_{n,B}^k, statHLCP use statHLCP() and statHLCPk(), respectively.
""" 


load("hssp.sage")

#Class hssp to generate instances both of HLCP_{n,B} and HSSP_{n,B}^k, either specifing k or not, respectilvely.
#For HLCP_n use:
#H=hlcp(n) 
#H.gen_instance()
#For HLCP_n^kappa use:
#H=hssp(n,kappa)
#H.gen_instance()
#This generate the modulus H.x0, the weights vector H.a, the matrix H.x, the  sample vector H.b
#kappa=-1 is HLCP_n by construction
class hlcp:
  def __init__(self,n,B,kappa=-1):
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
    
 # This is the the function to perform the attacks.    
#H is the instance to be attacked and alg is the algorithm to use :
#       **if alg='default' or alg='statistica' runs the statistical attack
#       **if alg='ns_original' runs the original Nguyen-Stern attack 
#       **if alg='ns' runs the Nguyen-Stern attack with the improved orthogonal lattice attack
 
def hlcp_attack(H,alg='default'):
  # If B=1, it calls hssp_attack so actually in this case.
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
      beta,tt2,nrafound=ns(H,MO,B)
      return tt1,tt10,tt1O,beta,tt2, nrafound, H       
      
  if alg=='ns':
      MO,tt1,tt10,tt1O= Step1(H.n,H.kappa,H.x0,H.a,H.X,H.b,H.m,BKZ=True)
      beta,tt2,nrafound=ns(H,MO,B)
      return tt1,tt10,tt1O,beta,tt2,nrafound, H

  return  
  
  
def statHLCP(B=2,mi=70, ma=150, M=0,A='default'):
  L=[]
  for n in range(mi,ma,20):
    print
    H=hlcp(n,B,-1)
    H.gen_instance(m=M*n)
    L+=[hlcp_attack(H,A)]
    print 
    print L
  return L 
  
def statHLCPk(sigma,B=2,mi=70, ma=150, M=0,A='default'):
  L=[]
  for n in range(mi,ma,20):
    print
    kappa=int(sigma*n)
    H=hlcp(n,B,kappa)
    H.gen_instance(m=M*n)
    L+=[hlcp_attack(H,A)]
    print 
    print L
  return L 
  
def print_stat_hlcp(L):
  Lp=[]
  for x in L:
    H=x[-1]   
    Lp+=[[H.n,H.kappa,H.m, H.nx0]+list(x[:-1])]   
  print table_s(Lp)    
  
