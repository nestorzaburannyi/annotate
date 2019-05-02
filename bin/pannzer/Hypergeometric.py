#!/usr/bin/env python
#######################################################################
#@version: 1.0
#@date: 12.07.2010
#@author: jussi.nokso-koivisto@helsinki.fi, patrik.koskinen@helsinki.fi
#@organization: Liisa Holm's Bioinformatics Group/University of Helsinki
#@license: This file is part of PANNZER.
#
#PANNZER is free software: you can redistribute it and/or modify
#it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or
#(at your option) any later version.
#
#PANNZER is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#GNU General Public License for more details.
#
#You should have received a copy of the GNU General Public License
#along with PANNZER.  If not, see <http://www.gnu.org/licenses/>.
#######################################################################

#Standard python modules
from __future__ import print_function
from operator import truediv

#Numpy/scipy
#import numpy as np
from numpy import exp, log, cumsum, max
from scipy.special import gammaln
from scipy.stats import hypergeom
from math import sqrt
import numpy as np

#import matplotlib.pyplot as plt

def calculate_p_value_for_hypergeometric(x, m, k, n, debug=False, hyper_p_cache={}):
    """
    Calculates p-value for hypergeometric distribution.

    :param x: represents the observed number of class members in the sample

    :param m: represents the size of whole dataset

    :param k: represents the size of the sample taken from the dataset

    :param n: represents the size of class in the dataset

    :param debug: (boolean) if True, the debug information will be printed

    :param hyper_p_cache: cache for speeding up the searches

    :rtype: p-value [0,1]
    """
    try: #check if can be found from cache
        return hyper_p_cache[x, m, k, n]
    except KeyError:
        pass

    if x > n:
        print("ERROR Hypergeometric.calculate_p_value_for_hypergeometric x > n", x, n, k, m)
        return None

    Hyper_mean = expectation_value(n,k,m)
    Hyper_STD  = sqrt( variance( n,k,m, Hyper_mean))
    if x < Hyper_mean + Hyper_STD:  # arbitrary cutoff. mean + STD to speed up things
       result = 1
    else:
      cum_weights = calculate_cumulative_sum_vector(x, m, k, n, debug)
      sum_of_exp = 0.0
      constant = max(cum_weights)

#      for value in cum_weights:
#          sum_of_exp += exp(value - constant)
      sum_of_exp = np.sum(np.exp(np.array(cum_weights) - constant))
      try:
          result = exp(log(sum_of_exp)+ constant)
      except OverflowError as e:
          print("Hypergeometric.calculate_p_value_for_hypergeometric OverflowError", e, "x=", x, "m=", m, "k=", k, "n=", n)
    
    # IF ELSE ends
    hyper_p_cache[x, m, k, n] = result

    if debug:
        print("sum_of_exp", sum_of_exp)
        print("log(sum_of_exp)", log(sum_of_exp))
        print("exp(log(sum_of_exp)+ constant)", exp(log(sum_of_exp)+ constant))

    return result


def calculate_cumulative_sum_vector(x, m, k, n, debug=False):
    """
    Calculates cumulative sum vector for hypergeometric distribution.

    :param x: represents the observed number of class members in the sample

    :param m: represents the size of whole dataset

    :param k: represents the size of the sample taken from the dataset

    :param n: represents the size of class in the dataset

    :rtype: cumulative sum vector (numpy array)
    """

    maximum_value = hypergeometric_max_value(m,k,n) #Calculate maximum value for hypergeometric distribution
    if debug:
        print("maximum_value", maximum_value)

    if n < k:
        vector = range(x, n)[::-1]
    else:
        vector = range(x, k)[::-1]

    if debug:
        print("vector", vector.__str__())

    weights = count_weights(vector, m, k, n, maximum_value)
    if debug:
        print("weights", weights)
        print("len(weights)", len(weights))

    cum_weights = cumsum(weights)
    if debug:
        print("cum_weights", cum_weights)
        print("len(cum_weights)", len(cum_weights))

    return cum_weights

def hypergeometric_max_value(m,k,n,hyper_max_cache={}):
    """
    Function used to calculate the probability for the largest outcome
    from the hypergeometric distribution. This is used to calculate the
    remaining probabilities as a chain. This diminishes dramatically the
    need of gamma functions from the calculus and therefore speed the calculus.

    :param m: represents the size of whole dataset

    :param k: represents the size of the sample taken from the dataset

    :param n: represents the size of class in the dataset

    :rtype: largest outcome from the hypergeometric distribution
    """

    try: #check if can be found from cache
        return hyper_max_cache[m,k,n]
    except KeyError:
        pass

    if k > n:
        aa = m - n + 1
        bb = k - n + 1
        cc = k + 1
    else:
        aa = m - k + 1
        bb = n - k + 1
        cc = n + 1

    result = gammaln(cc) + gammaln(aa) - gammaln(bb) - gammaln(m + 1)
    hyper_max_cache[m,k,n] = result
    return result

def count_weights(vector, m, k, n, maximum_value):
    """
    Calculates the weight for each outcome of the hypergeometric distribution.
    Weighting is used to obtain probability of each bin by multiplying the
    result for the bin with the larger outcome with the weight.

    :param vector: tail vector of given hypergeometric distribution

    :param m: represents the size of whole dataset

    :param k: represents the size of the sample taken from the dataset

    :param n: represents the size of class in the dataset

    :param maximum_value: maximum value of given hypergeometric distribution

    :rtype: array of weights
    """
    vector = np.array(vector)
    result = np.log(vector + 1) + np.log(m - k - n + vector + 1) - np.log(n - vector) - np.log(k - vector)
    return [maximum_value] + list(result)

#    result = [maximum_value]
#    for value in vector:
#        result.append(log(value + 1) + log(m - k - n + value + 1) - log(n - value) - log(k - value))
#
#    return result

def expectation_value(n, k, m, hyper_e_cache={}):
    """
    Returns expectation value of class count in given hypergeometric distribution

    :param m: represents the size of whole dataset

    :param k: represents the size of the sample taken from the dataset

    :param n: represents the size of class in the dataset

    :param hyper_e_cache: cache to speed up the calculations

    :rtype: expectation value of class count
    """

    try: #check if can be found from cache
        return hyper_e_cache[m,k,n]
    except KeyError:
        pass

    result = truediv((k*n), m)
    hyper_e_cache[m,k,n] = result

    return result

def variance(n, k, m, E, hyper_v_cache={}):
    """
    Returns variance of class count in given hypergeometric distribution

    :param m: represents the size of whole dataset

    :param k: represents the size of the sample taken from the dataset

    :param n: represents the size of class in the dataset

    :param E: represents the expectation value of the hypergeometric distribution

    :rtype: variance of class count
    """
    try: #check if can be found from cache
        return hyper_v_cache[m,k,n,E]
    except KeyError:
        pass

    tmp1 = truediv((m - k), m)
    tmp2 = truediv((m - n), (m - 1))
    variance = E*tmp1*tmp2
    hyper_v_cache[m,k,n,E] = variance

    return variance

def hypergeometric_pmf(m, k, n, debug=False):
    """
    Calculates and returns probability mass function for given hypergeometric distribution.
    This function also returns the corresponding number of positive cases in the sample.

    :param m: represents the size of whole dataset

    :param k: represents the size of the sample taken from the dataset

    :param n: represents the size of class in the dataset

    :rtype: tuple (values, columns) where values is a list of the probabilities and columns
            is a list of corresponding number of positive cases in the sample.
    """
    rv = hypergeom(m, n, k)
    values = []
    columns = []

    if n < k:
        vector = reversed(range(0, n + 1))
    else:
        vector = reversed(range(0, k + 1))

    for i in vector:
        columns.append(i)
        values.append(rv.pmf(i))

    #if debug and sum(values) != 1.0:
    #    print "Hypergeometric.hypergeometric_pmf sum of the probabilities not 1.0"

    return values, columns

def calculate_pmf_value_for_hypergeometric(x, m, k, n, debug=False, hyper_p_cache={}):
    """
    Calculates p-value for hypergeometric distribution.

    :param x: represents the observed number of class members in the sample

    :param m: represents the size of whole dataset

    :param k: represents the size of the sample taken from the dataset

    :param n: represents the size of class in the dataset

    :param debug: (boolean) if True, the debug information will be printed

    :param hyper_p_cache: cache for speeding up the searches

    :rtype: p-value [0,1]
    """
    try: #check if can be found from cache
        return hyper_p_cache[x, m, k, n]
    except KeyError:
        pass

    if x > n:
        print("ERROR Hypergeometric.calculate_p_value_for_hypergeometric x > n", x, n, k, m)
        return None

    maximum_value = hypergeometric_max_value(m,k,n) #Calculate maximum value for hypergeometric distribution
    if debug:
        print("maximum_value", maximum_value)

    if n < k:
        vector = range(x, n)[::-1]
    else:
        vector = range(x, k)[::-1]

    if debug:
        print("vector", vector.__str__())

    weights = count_weights(vector, m, k, n, maximum_value)
    print("weights", weights)
    sum_of_exp = 0.0
    constant = max(weights)

    for value in weights:
        sum_of_exp += exp(value - constant)

    try:
        result = exp(log(sum_of_exp)+ constant)
        hyper_p_cache[x, m, k, n] = result
    except OverflowError as e:
        print("Hypergeometric.calculate_p_value_for_hypergeometric OverflowError", e, "x=", x, "m=", m, "k=", k, "n=", n)

    if debug:
        print("sum_of_exp", sum_of_exp)
        print("log(sum_of_exp)", log(sum_of_exp))
        print("exp(log(sum_of_exp)+ constant)", exp(log(sum_of_exp)+ constant))

    return result


