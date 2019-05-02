from __future__ import absolute_import
from __future__ import print_function
from operator import truediv
from numpy import sqrt

import Hypergeometric

def calculateGSZscore(sampleSize, score, sampleScoreMean, sampleScoreVariance, sizeInBackGroud, backGroundSize):
    """
    Math is represented in Toronen, Ojala, Marttinen, Holm (2009): Robust extraction of functional signals ...
    http://www.biomedcentral.com/1471-2105/10/307/
    See eq.1 - eq.4

    :param sampleSize: size of the sample taken
    :param score: element's score to be normalized with GSZ
    :param sampleScoreMean: mean of the score inside of the sample
    :param sampleScoreVariance: variance of the score inside of the sample
    :param sizeInBackGroud: size of the element in the background
    :param backGroundSize: size of the total background
    :rtype: GSZ score
    """
    if sizeInBackGroud < 2:
        return 0
    if sizeInBackGroud < 10:
        sizeInBackGroud = 10

    E = Hypergeometric.expectation_value(sizeInBackGroud, sampleSize, backGroundSize)
    V = Hypergeometric.variance(sizeInBackGroud, sampleSize, backGroundSize, E)

    if sampleSize > 1:
        E_GSZ = sampleScoreMean*E
        D2_GSZ = truediv((sampleScoreVariance*(E*(sampleSize - E) - V)), (sampleSize - 1)) + sampleScoreMean**2*V
        GSZ = truediv((score - E_GSZ), sqrt(D2_GSZ))
        return GSZ
    else:
        GSZ = E*V
        return GSZ


def testCalculateGSZscore():
    gsz = calculateGSZscore(249, 1955.5385392, 3.01751184858, 0.311070695203, 666, 11384898)

    diff = abs(gsz-float(5280.29408996))

    if diff < 0.000001:
        print("TEST OK")
    else:
        print("TEST FAILED Should be 5280.29408996, is", gsz,diff)

if __name__ == '__main__':
    testCalculateGSZscore()
    print(calculateGSZscore(52,38.6855121241,0.743952156233,0.00103566165115,92672205,92672206))
