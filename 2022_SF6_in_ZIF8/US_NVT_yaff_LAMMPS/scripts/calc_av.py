#!/usr/bin/env python
import os
import numpy as np
from molmod.units import angstrom

fn_name = '1SF6'

print('Simulations with K = 25 kJ/mol')
CV_all = np.loadtxt('../results/CVs_%s.csv' %fn_name, delimiter=',')
CV_av = np.mean(CV_all,axis=1)
for i in range(len(CV_av)):
    print('Target %.1f A, actual average %.1f A' %(-11+i*0.5, CV_av[i]/angstrom))

print('\n\nSimulations with K = 50 kJ/mol')
CV_all = np.loadtxt('../results/CVs_%s_50.csv' %fn_name, delimiter=',')
CV_av = np.mean(CV_all,axis=1)
for i in range(len(CV_av)):
    print('Target %.1f A, actual average %.1f A' %(-3+i*0.5, CV_av[i]/angstrom))

print('\n\nSimulations with K = 75 kJ/mol')
CV_all = np.loadtxt('../results/CVs_%s_75.csv' %fn_name, delimiter=',')
CV_av = np.mean(CV_all,axis=1)
for i in range(len(CV_av)):
    print('Target %.1f A, actual average %.1f A' %(-3+i*0.5, CV_av[i]/angstrom))

print('\n\nSimulations with K = 100 kJ/mol')
CV_all = np.loadtxt('../results/CVs_%s_100.csv' %fn_name, delimiter=',')
CV_av = np.mean(CV_all,axis=1)
for i in range(len(CV_av)):
    print('Target %.1f A, actual average %.1f A' %(-3+i*0.5, CV_av[i]/angstrom))

print('\n\nSimulations with K = 150 kJ/mol')
CV_all = np.loadtxt('../results/CVs_%s_150.csv' %fn_name, delimiter=',')
CV_av = np.mean(CV_all,axis=1)
for i in range(len(CV_av)):
    print('Target %.1f A, actual average %.1f A' %(-3+i*0.5, CV_av[i]/angstrom))

print('\n\nSimulations with K = 200 kJ/mol')
CV_all = np.loadtxt('../results/CVs_%s_200.csv' %fn_name, delimiter=',')
CV_av = np.mean(CV_all,axis=1)
for i in range(len(CV_av)):
    print('Target %.1f A, actual average %.1f A' %(-3+i*0.5, CV_av[i]/angstrom))

print('\n\nSimulations with K = 300 kJ/mol')
CV_all = np.loadtxt('../results/CVs_%s_300.csv' %fn_name, delimiter=',')
CV_av = np.mean(CV_all,axis=1)
for i in range(len(CV_av)):
    print('Target %.1f A, actual average %.1f A' %(-3+i*0.5, CV_av[i]/angstrom))

print('\n\nSimulations with K = 400 kJ/mol')
CV_all = np.loadtxt('../results/CVs_%s_400.csv' %fn_name, delimiter=',')
CV_av = np.mean(CV_all,axis=1)
for i in range(len(CV_av)):
    print('Target %.1f A, actual average %.1f A' %(-3+i*0.5, CV_av[i]/angstrom))

