#!/usr/bin/python

import numpy as np
from sklearn.decomposition import FastICA
import math
from time import time

def ICA(X,B=1,n=16,kappa=-1,exn=0):
  t=time()
  ncomp,sample=X.shape
  ica=FastICA(n_components=ncomp,fun='cube')
  S_=ica.fit_transform(X.T)
  A_ = ica.mixing_
  assert np.allclose(X.T, np.dot(S_, A_.T) + ica.mean_)

  # we want to remove the mean
  s0=np.dot(ica.mean_,np.linalg.inv(A_.T))
  S_+=s0

  assert np.allclose(X, np.dot(A_,S_.T))  # now the mean is 0

  # we want that the components of S_ are positive
  # if the average of a column of S_ is negative, we take the opposite

  Sm=np.sign(np.average(S_,axis=0))
  S_*=Sm
  A_*=Sm

  assert np.allclose(X, np.dot(A_,S_.T))

  # we want that the components of S_ are 0 or 1
  # If all the components of S_ are 0 or 1, the standard deviation will be 1/2
  # If the components are uniformly distributed between 0 and B, the standard deviation is
  # sqrt(B(B+2)/12)
  # So we divide the components of S_ by std/sqrt(B(B+2)/12), where std is computed column by column
  if kappa==-1 and exn==0:  exn=math.sqrt(B*(B+2)/12.)
  elif exn==0: exn=math.sqrt((2*B+1)*(B+1)*kappa/(n*6.)-(kappa*(B+1)/(2*n))**2)  # exn=math.sqrt((2*B+1)*(B+1)*kappa/(n*6.))
  
  st=S_.std(axis=0)/exn
  S_/=st
  A_*=st
  
  #print " ICA.py %.2f" %(time()-t),
  #assert np.allclose(X, np.dot(A_,S_.T))
  return A_,S_.T

