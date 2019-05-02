from __future__ import print_function
#################################
# This is a collection of
# GO distance related functions
#
# These were not coded for
# SansPanz-project originally
# Therefore, the suitability
# for that project might vary
#################################

import sys, re
import random
# import Util, GO, DB
from math import log, exp
from os import path
import cPickle as pickle
import numpy as np
import copy
############################
# Stuff Debug Commands belo#
############################
#from python_diagnostics import * # debug commands. Exclude when not needed
#import pprint

def mod_parents(GO_parents):

    parents_out = {}
    for ID in GO_parents:
       parents_out[ID] = GO_parents[ID].split(",")
    #print_variables2(vars())
    #print_some_dict2(parents_out)
    return parents_out

def get_jaccard( GO_set1, GO_set2, GO_parents,
                 add_penalty = 0, debug = False):

  def process_input(GO_set):
    # If input is a string (single GO ID
    # turn it into a list with single item
    # Filter empty items from the lists
    type_name = type(GO_set).__name__
    if not type_name == 'list':  #Case GO_set is 1 ID
      GO_set = [GO_set]
    GO_set = filter(None, GO_set)
    return GO_set

  def pool_parents(GO_set, GO_parents):

    #set = []
    set_out = copy.copy( GO_set)
    #set_out = set(GO_set)
    #print "Pool parents alkaa"
    #print_variables2(vars())
    #print GO_set
    for ID in GO_set:
      try:
        set_out.extend( GO_parents[ID] )
      except KeyError:
        pass
    #FOR
    #print "SUBFUNCTION"
    #print_variables2(vars())
    return set_out
  # SUBFUNCTIONS ENDED
  # MAIN CODE STARTS

  if debug:
    print("Before")
    print_variables2(vars())
    print(GO_set1)
    print(GO_set2)

  GO_set1 = process_input(GO_set1)
  set1 = pool_parents(GO_set1, GO_parents)

  GO_set2 = process_input(GO_set2)
  set2 = pool_parents(GO_set2, GO_parents)

  if debug:
    print("after")
    print_variables2(vars())
    print(GO_set1)
    print(GO_set2)

  set1 = filter(None, set1)
  set2 = filter(None, set2)

  set1 = set(set1)
  set2 = set(set2)
  if debug:
    print(set1)
    print(set2)

  len_S1 = len(set1)
  len_S2 = len(set2)
  len_inter = len( list( set1 & set2))
  jacc_out = len_inter/float(max(len_S1 + len_S2 - len_inter + add_penalty, 1))
  return jacc_out

def calc_max_of_jacc( GO_class, GO_list, GO_parents, add_prior = 0, debug = False):

   """
   Calculate Jacc of GO class against all classes in list
   REport the best
   Stored_Jacc is updated, when new pair is observed
   """

   #print_variables2(vars())
   jacc_out = -1 # Unrealistic value. This signals some error if gets out
   for test_GO in GO_list:
      #print test_GO
      #print GO_class
      tmp = get_jaccard(GO_class, test_GO, GO_parents, add_prior, debug = debug)
      jacc_out = max(jacc_out, tmp)
   # FOR
   if debug:
     debug_vec = []
     for test_GO in GO_list:
       debug_vec.append( get_jaccard(GO_class, test_GO, GO_parents, add_prior))
     print(debug_vec)
   return jacc_out

def Create_ClassLogScores(Data_in, RefCount = [], Conv2float = False):

   """
   Create Information Content scores
   """

   if not RefCount:
     RefCount = 0
     for key, val in Data_in.iteritems():
        RefCount = max(RefCount,  float(val))
   # IF
   AvoidZeroScore = 0.000000001  # Avoid zero-division in later analysis
   #log_refcount = log(RefCount + AvoidZeroScore)
   log_refcount = log(RefCount)
   print("Max arvot")
   print(RefCount, log_refcount)
   Data_out1 = {}
   Data_out2 = {}
   for key, val in Data_in.iteritems():
      if Conv2float:
         val = float(val)
      Data_out1[key] = -log(float(val + AvoidZeroScore))
      Data_out2[key] = log_refcount + Data_out1[key]
   # FOR
   return (Data_out1, Data_out2)


def get_GO_dists( GO_set1, GO_set2, GO_parents, InfContScores, DB_size,
                  add_penalty = 0, error_value = [],  debug = False):

  def process_input(GO_set):
    # If input is a string (single GO ID
    # turn it into a list with single item
    # Filter empty items from the lists
    type_name = type(GO_set).__name__
    if not type_name == 'list':  #Case GO_set is 1 ID
      GO_set = [GO_set]
    GO_set = filter(None, GO_set)
    return GO_set

  def pool_parents(GO_set, GO_parents):

    set = []
    for ID in GO_set:
      try:
        set.extend( GO_parents[ID] )
      except KeyError:
        pass
    #FOR
    return set

  def GetDictScore(ID, Dict, error_value):

     try:
       out = Dict[ID]
     except KeyError:
       out = error_value
     return out

  def MaxOfDictSubset(dict_in, key_list):

     out = max([ dict_in[i] for i in key_list if i in dict_in])
     return out

  def MinOfDictSubset(dict_in, key_list):

     out = min([ dict_in[i] for i in key_list if i in dict_in])
     return out

  def MaxOfDictVals(dict_in):

     out = max(dict_in.values())
     return out

  def MinOfDictVals(dict_in):

     out = min(dict_in.values())
     return out

  def PrintDictSubset(key_set, Dict_in):

    for key in key_set:
       try:
         print(key, Dict_in[key])
       except KeyError:
         print(key, "Not Found")
    #FOR


  def SumOfDictVals(key_set, Dict_in, val_if_err = [], debug = False):

    if not val_if_err:
       val_if_err = MaxOfDictVals(Dict_in)

    val_out = 0
    for key in key_set:
       try:
         val_out += Dict_in[key]
         if debug:
           print(Dict_in[key])
       except KeyError:
         val_out += val_if_err
    #FOR
    return val_out

  #######################
  # SUBFUNCTIONS ENDED
  # MAIN CODE STARTS
  #######################
  # Code is modified so that
  # it processes only single GO IDs

  GO_set1 = [ GO_set1 ]
  GO_set2 = [ GO_set2 ]

  set1 = pool_parents( GO_set1, GO_parents)
  #set1 = pool_parents( GO_set1, GO_parents)
  set2 = pool_parents( GO_set2, GO_parents)

  if debug:
    print_variables2(vars())
    print(GO_set1)
    print(GO_set2)

  set1 = filter(None, set1)
  set2 = filter(None, set2)

  set1 = set(set1)
  set2 = set(set2)
  if debug:
    print(set1)
    print(set2)

  if list(set1 & set2):  #CAse: Union set exists.
    if not error_value:
      error_value = MaxOfDictVals(InfContScores)

    Score_GO1 = GetDictScore(GO_set1[0], InfContScores, error_value)
    Score_GO2 = GetDictScore(GO_set2[0], InfContScores, error_value)

    sum_S1_only = SumOfDictVals( set1 -set2 | set(GO_set1), InfContScores,  val_if_err = error_value)
    sum_S2_only = SumOfDictVals( set2 -set1 | set(GO_set2), InfContScores,  val_if_err = error_value)

    #print "Test prints"
    #PrintDictSubset( set1 -set2 | set(GO_set1), InfContScores)
    #print "second"
    #PrintDictSubset( set1 -set2 | set(GO_set2), InfContScores)
    #print "Kolmas"
    #PrintDictSubset( set1 & set2, InfContScores)
    #print list( set1 & set2)
    #print MaxOfDictSubset(InfContScores, list( set1 & set2) )
    #print "Siina oli"
    Max_of_inters = float( MaxOfDictSubset(InfContScores, list( set1 & set2) ) )   #This is Resnik Score
    if ( Max_of_inters >= 0):
      weight = (1 - exp(-Max_of_inters)/DB_size)
    else:
       weight = (1 - exp(-Max_of_inters)/DB_size)
    #print weight
    #if Score_GO1 == 0 or Score_GO2 == 0:
    #  print "Error case"
    #  print_variables2(vars())
    #  print GO_set1
    #  print GO_set2
    #  print Score_GO1
    #  print Score_GO2
    if Max_of_inters > min( [Score_GO1, Score_GO2]):    #THIS IS BAD CORRECTION. PARENT-CHILD LINKS CORRUPTED??
       Max_of_inters = min( [Score_GO1, Score_GO2])

    Lin_Score = 2*Max_of_inters/(Score_GO1 + Score_GO2 )
    Path_Lin  = 2*Max_of_inters/(sum_S1_only + sum_S2_only )
    Wght_Lin_Score = weight*Lin_Score
    Wght_Path_Lin  = weight*Path_Lin

    if debug:                              # Following is the debug print for failure cases
      if Lin_Score > 1.5 or Path_Lin > 1.5:
         print("Test prints")
         PrintDictSubset( set1 -set2 | set(GO_set1), InfContScores)
         print("second")
         PrintDictSubset( set1 -set2 | set(GO_set2), InfContScores)
         print("Kolmas")
         PrintDictSubset( set1 & set2, InfContScores)
         print(list( set1 & set2))
         print(MaxOfDictSubset(InfContScores, list( set1 & set2) ))
         print("Siina oli")

      #print "This step"
      #print Max_of_inters
      #print sum_S1_only, sum_S2_only
      #print Score_GO1, Score_GO2
      #pprint.pprint()

  else:
    Lin_Score = 0
    Path_Lin  = 0
    Wght_Lin_Score = 0
    Wght_Path_Lin  = 0
  return (Lin_Score, Path_Lin, Wght_Lin_Score, Wght_Path_Lin)

def calc_max_go_dists( GO_class, GO_list, GO_parents, Inf_Cont, db_size, error_case_value = [], debug = False ):

   """
   Calculate GO dists of GO class against all classes in list
   REport the best
   Stored_Dist is updated, when new pair is observed
   """

   #print_variables2(vars())
   jacc_out = [-1, -1, -1, -1] # Unrealistic value. This signals some error if gets out
   for test_GO in GO_list:
      tmp = get_GO_dists(GO_class, test_GO, GO_parents, Inf_Cont, db_size, error_value = error_case_value, debug = debug)
      #print tmp
      jacc_out = [ max(jacc_out[i], tmp[i]) for i in range(4) ]
   # FOR
   if debug:
     debug_vec = []
     #for test_GO in GO_list:
     #
     #  #debug_vec.append( get_jaccard(GO_class, test_GO, GO_parents, add_prior))
     print(debug_vec)
   return jacc_out
