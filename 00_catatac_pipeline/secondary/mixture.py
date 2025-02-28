import pandas as pd
import numpy as np
import csv
import os
import matplotlib.pyplot as plt
from matplotlib import gridspec
from pomegranate import *

def pgmix_mini(data):
    data = np.array(data)
    data = data.reshape(-1,1) 
    model = GeneralMixtureModel.from_samples([PoissonDistribution,NormalDistribution], 2, data)
    return model.predict_proba(data)


def ggmix(data):
    data = np.array(data)
    data = data.reshape(-1,1) 
    model = GeneralMixtureModel.from_samples([NormalDistribution,NormalDistribution], 2, data)
    return model.predict_proba(data)

def test_func(data):
    return data

def pgmix(data):
    data = np.log2(np.array(data))
    data = data.reshape(-1,1) 
    
    ##this loop checks that the model converges 
    ##and that the Poisson distribution is the lower component 
    i=0
    gmm_x = np.linspace(-2,max(data)+2,1000)
    while i==0:
        model = GeneralMixtureModel.from_samples([PoissonDistribution,NormalDistribution],2,data)
        if numpy.isnan(model.probability(gmm_x)).any():
            i=0
        else:
            if model.distributions[0].name == 'PoissonDistribution':
                if model.distributions[0].parameters[0]<model.distributions[1].parameters[0]:
                    i=1    
                else:
                    i=0
            elif model.distributions[0].parameters[0]>model.distributions[1].parameters[0]:
                    i=0
    
    ##append positive calls to table
    if model.distributions[0].name == 'PoissonDistribution':
        return model.predict_proba(data)
    else:
        return 1 - model.predict_proba(data)
