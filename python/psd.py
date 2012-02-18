#!/usr/bin/env python 
from scipy import fftpack
import numpy as np
from matplotlib import pyplot as plt
import argparse

TEST = "icam_hmc_omp.txt"

def psd(x):
  a = fftpack.fft(x)
  a = abs(a)/len(a)
  return (a**2)

def main():
  parser = argparse.ArgumentParser(description='Plot psd of frequency file')
  parser.add_argument('-i', dest='inp', type=str, help='File for input')
  args = parser.parse_args()
  x = np.genfromtxt(args.inp)
  plt.plot(psd(x))
  plt.show()
  
if __name__ == '__main__':
  main()
