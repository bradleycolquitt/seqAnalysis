#!/usr/bin/env python

import numpy as np
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
import itertools
from scipy import stats
from scipy.stats import gaussian_kde
from math import *
import argparse

DATA_NAME = "genes.fpkm_tracking"
#SUMMARY_PATH_NORM = "/media/storage2/analysis/features/norm/summaries"
#SUMMARY_PATH_RAW = "/media/storage2/analysis/features/norm/rpkm/mean/summaries"
SUMMARY_PATH_RAW = "/media/storage2/analysis/rna/summaries"
INF = float('inf')

# from MISO
# I think this is jjust a hack to get around a bug in gaussian_kde, but not sure
class gaussian_kde_covfact(gaussian_kde):
  def __init__(self, dataset, covfact = 'scotts'):
    self.covfact = covfact
    gaussian_kde.__init__(self, dataset)

  def _compute_covariance_(self):
    self.inv_cov = np.linalg.inv(self.convariance)
    self._norm_factor = sqrt(np.linalg.det(2*np.pi*self.covariance)) * self.n

  def covariance_factor(self):
    if self.covfact in ['sc', 'scotts']:
      return self.scotts_factor()
    if self.covfact in ['si', 'silverman']:
        return self.silverman_factor()
    elif self.covfact:
      return float(self.covfact)

  def reset_covfact(self, covfact):
    self.covfact = covfact
    self.covariance_factor()
    self._compute_covariance()

class NullPeakedDensity:
  def __init__(self, data):
    self.data = data
  def evaluate(self, point):
    if point[0] == 0:
      return INF
    else:
      return 0

class ExprData:
  def __init__(self, d):
    num_data = len(d.keys())
    self.wt_expr = np.zeros([num_data, 1])
    self.ko_expr = np.zeros([num_data, 1])
    i = 0
    for name, expr in d.iteritems():
      self.wt_expr[i] = expr[0]
      self.ko_expr[i] = expr[1]
      i = i+1

  def plot(self):
    wtlog2 = np.log2(self.wt_expr[np.nonzero(self.wt_expr)])
    kolog2 = np.log2(self.ko_expr[np.nonzero(self.ko_expr)])
    plt.subplot(2, 2, 1)
    plt.hist(wtlog2,100)
    plt.subplot(2, 2, 2)
    plt.hist(kolog2, 100)
    plt.subplot(2, 2, 3)
    plt.scatter(np.log2(self.wt_expr + 1),np.log2(self.ko_expr + 1), s=1)
    plt.subplot(2,2,4)
#plt.scatter(foldchange(self.wt_expr, self.ko_expr), self.compute_bayes_factor(),100)
    plt.show()

  # from the MISO paper: 
  # BF = p(D|H1)p(H1) / p(D|H0)p(H0)
  # BF approx= P(d = 0|H1) / p(d = 0|D, H1)
  # (P(d = 0|H1) is the uniform prior density, and p(d = 0|D,H1) is the KDE estimated posterior density)
  # aka Savage-Dickey density ratio
  # use a uniform prior over wt and ko
  # 
  # By the conjugacy of normal(u,s) with u known (in this case u = 0), and Gamma(a,b) prior on s
  # the posterior distribution of s is Gamma(a + n/2, (n-1)S^2), where S^2 is the sample variance.
  #
  # As such, the MAP estimate of s is 
  #
  # 
  def compute_bayes_factor(self):
    # testing the hypothesis that delta = 0
    # vs delta != 0

    # use a uniform prior on the posterior density
    maxbf = 1e12
    delta = self.wt_expr - self.ko_expr
    delta = delta - np.mean(delta)
    priorfun = lambda x: 1 + x if x <= 0 else 1 - x
    analytic_prior_density = map(priorfun, np.arange(-1,1,0.001))
    if np.mean(abs(delta)) <= 0.00009:
      #print "null peaked"
      postfun = NullPeakedDensity(delta)
    else:
      #print "not null peaked"
      postfun = self._fit_posterior()
    bf = np.zeros(delta.size)
    posterior = postfun.evaluate(np.transpose(delta))
    for i in range(bf.size):
      sign = 1
      if delta[i] < 0: sign = -1
      if posterior[i] == 0:
        bf[i] = maxbf * sign
      elif posterior[i] == INF:
        bf[i] = 0
      else:
        bf[i] = 1.0 / posterior[i] * sign
    return bf

  def plot_bayes_factors(self):
    bf = self.compute_bayes_factor()
    print sum([1 for x in bf if x > 500])
    print bf.size
    plt.hist(bf, 100)
    plt.show()

  def _fit_posterior(self):
    delta = self.wt_expr - self.ko_expr
    delta = delta - np.mean(delta)
#delta = np.zeros(self.wt_expr.size)
#    for i in range(self.wt_expr.size):
#      delta[i] = self.wt_expr[i] - self.ko_expr[i]
    return gaussian_kde_covfact(np.transpose(delta), 0.3)

def foldchange(d1, d2):
  fc = np.zeros(d1.size)
  for i in xrange(d1.size):
    if  d1[i] == 0 or d2[i] == 0:
      fc[i] = 0
    elif d1[i] > d2[i]:
      fc[i] = log(d1[i] / d2[i])
    elif d1[i] < d2[i]:
      fc[i] = - log(d2[i] / d1[i])
  return fc
 
def import_data(dname):
  f = open(dname, 'r')
  d = {}
  for line in f:
    if line[0] == "#":
      pass
    else:
      fields = line.split()
      if fields[0] == "tracking_id":  
        pass
      else: 
      # first field of tuple is WT, second field is KO
        d[fields[0]] =  (float(fields[10]), float(fields[13])) 
  return d

def import_features(set, feature, samples, data_type):
  d = {}
  sample_name = ""
  path = ""
  if data_type == "norm":
    path = SUMMARY_PATH_NORM
  elif data_type == "raw":
    path = SUMMARY_PATH_RAW
  sample_name = "/".join([path, "_".join([set, feature])])
  
  sample_data = open(sample_name)
  cnames = sample_data.readline().split()

  indices = [cnames.index(sample) + 1 for sample in samples]
  for line in sample_data:
    line = line.split()    
    d[line[0]] = (float(line[indices[0]]), float(line[indices[1]]))
  
  #for s1, s2 in itertools.izip(sample_data[0], sample_data[1]):
  #  s1 = s1.split()
  #  s2 = s2.split()
  #  d[s1[3]] = (float(s1[4]), float(s2[4]))

  return d

def worker(d):
  #n = d.keys()
  e = ExprData(d)
  bf = e.compute_bayes_factor()
  return bf
  #for i in range(len(bf)):
  #  if args.cutoff:
  #    if bf[i] > args.cutoff:
  #      print "%s\t%d" % (n[i], bf[i])
  #  else:  
  #    print "%s\t%d" % (n[i], bf[i])
  
def main():
  parser = argparse.ArgumentParser(description='Determine differentially expressed genes.')
  parser.add_argument('--fpkm', metavar='file.fpkm_tracking', type=str, help='name of gene fpkm file')
  parser.add_argument('--set')
  parser.add_argument('--feature', help='feature file', required=False)
  parser.add_argument('--data_type', required=False)
  parser.add_argument('--samples', nargs=2, help='samples to be compared', required=False)
  parser.add_argument('--cutoff', type=float, help='Cutoff for differential expression')
  args = parser.parse_args()
  d = {}
  if args.feature and args.samples:
    d = import_features(args.set, args.feature, args.samples, args.data_type)
  elif args.fpkm:  
    d = import_data(args.fpkm)
  n = d.keys()
  e = ExprData(d)
  bf = e.compute_bayes_factor()
  for i in range(len(bf)):
    if args.cutoff:
      if bf[i] > args.cutoff:
        print "%s\t%d" % (n[i], bf[i])
    else:  
      print "%s\t%d" % (n[i], bf[i])
        
if __name__ == '__main__':
  main()

