#!/usr/bin/python
"""
This is an implementation of the attacks for solving HSSP. The class hssp can be used to generate  instances both of HSSP_n and HSSP_n^k, either specifing k or not, respectilvely.

The function hssp_attack performs the attacks.

To run experiments with HSSP_n and HSSP_n^k statHSSP use statHSSP() and statHSSPk(), respectively.
""" 

load("building.sage") # loads the functions for building lattices and instances
load("ortho_attack.sage")# loads the functions for performing the orthogonal lattice attack
load("multi.sage")# loads the functions for performing the multivariate attack
load("ns.sage")# loads the functions for performing the Nguyen-Stern attack
load("statistical.sage")# loads the functions for performing the statistical attack
load("print_tables.sage")# loads the functions for printing in latex tables format 


#Class hssp to generate instances both of HSSP_n and HSSP_n^k, either specifing k or not, respectilvely.
#For HSSP_n use:
#H=hssp(n) 
#H.gen_instance()
#For HSSP_n^kappa use:
#H=hssp(n,kappa)
#H.gen_instance()
#This generate the modulus H.x0, the weights vector H.a, the matrix H.x, the  sample vector H.b
#kappa=-1 is HSSP_n by construction
#
class hssp:
  def __init__(self,n,kappa=-1):
    self.n=n
    self.kappa=kappa
    
    
  def gen_instance(self,m=0): # m use to specify dimension of the sample 
    if m==0 and self.n % 2==1:
      m=self.n*(self.n+3)/2 # n is odd
    elif m==0 and self.n % 2 ==0:
      m=self.n*(self.n+4)/2 # n is even
    self.m=int(m)
    
    print "n=",self.n,"m=",m,
    if self.kappa>-1: print "kappa=",self.kappa,
    iota=0.035
    self.nx0=int(2*iota*self.n^2+self.n*log(self.n,2))# this is lower bound for log(Q) 
    print "nx0=",self.nx0
    # genParams returns the modulus x0, the weights vector a, the matrix x, the  sample vector b
    self.x0,self.a,self.X,self.b=genParams(self.n, self.m,self.nx0,self.kappa)

# This is the the function to perform the attacks.    
#H is the instance to be attacked and alg is the algorithm to use :
#       **if alg='default' or alg='multi' runs the multivariate attack
#       **if alg='ns_original' runs the original Nguyen-Stern attack 
#       **if alg='ns' runs the Nguyen-Stern attack with the improved orthogonal lattice attack
#       **if alg='statistical' runs the heuristic statistical attack 
    
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
    t=cputime()
    MO,tt1,tt10,tt1O= Step1(H.n,H.kappa,H.x0,H.a,H.X,H.b,H.m,BKZ=False)
    print "\nMultivariate Attack"
    if kappa>0:
      tei, tef, tsf,tt2, nrafound=bit_guessing(H.n,H.kappa,MO,H.x0,H.a,H.X,H.b,H.m) 
      tttot=cputime(t)
      return  tt1,tt10,tt1O, tei, tef, tsf, tt2, tttot, nrafound, H
    else:
      tei, tef, tt2,nrafound=eigen(H.n,H.kappa,MO,H.x0,H.a,H.X,H.b,H.m) 
      tttot=cputime(t) 
      return  tt1,tt10,tt1O, tei, tef, tt2,tttot,nrafound, H
  
  if alg=='ns_original' or (alg=='ns' and H.m==2*n):
      print "Nguyen-Stern (Original) Attack"
      t=cputime()
      MO,tt1,tt10,tt1O= Step1_original(H.n,H.kappa,H.x0,H.a,H.X,H.b,H.m)
      beta,tt2, nrafound,textra=ns(H,MO)
      tttot=cputime(t)
      return tt1,tt10,tt1O,beta,tt2,textra, tttot,nrafound, H   
      
  if alg=='ns':
      t=cputime()
      print "\nNguyen-Stern (Improved) Attack"
      MO,tt1,tt10,tt1O= Step1(H.n,H.kappa,H.x0,H.a,H.X,H.b,H.m,BKZ=True)
      beta,tt2, nrafound,textra=ns(H,MO)
      tttot=cputime(t)
      return tt1,tt10,tt1O,beta,tt2, textra ,tttot,nrafound, H
     
      
  if alg=='statistical':
    assert kappa==-1, 'The statistical attack does not work for HSSP_n^kappa(m,Q)'
    t=cputime()
    MO,tt1,tt10,tt1O= Step1(H.n,H.kappa,H.x0,H.a,H.X,H.b,H.m,BKZ=True)
    print "\nStatistical Attack"      
    tica, tt2, nrafound=statistical(MO,H.n,H.m,H.x0,H.X,H.a,H.b,H.kappa)
    tttot=cputime(t)
    return  tt1, tica, tt2,tttot, nrafound, H
  
  return None
  
  
def statHSSP(mi=70, ma=150, M=0,A='default',tablepr=True):
  L=[]
  for n in range(mi,ma,20):
    print
    H=hssp(n,-1)
    H.gen_instance(m=M*n)
    L+=[hssp_attack(H,A)]
    print 
    if tablepr: print print_stat(L)
  return L 
  
def statHSSPk(sigma,mi=70, ma=150, M=0,A='default',tablepr=True):
  L=[]
  for n in range(mi,ma,20):
    print
    kappa=int(sigma*n)
    H=hssp(n,kappa)
    H.gen_instance(m=M*n)
    L+=[hssp_attack(H,A)]
    print 
    if tablepr: print print_stat(L)
  return L 
  
def print_stat(L):
  Lp=[]
  for x in L:
    H=x[-1]   
    Lp+=[[H.n,H.kappa,H.m, H.nx0]+list(x[:-1])]   
  print table_s(Lp)   
