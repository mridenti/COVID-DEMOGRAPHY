# -*- coding: utf-8 -*-
import sympy as sym
import numpy as np
import sys
import os
from scipy import linalg
import pandas as pd


def compute_r0(input_folder, contact_matrix, flag_where):
    if flag_where:
        curr_dir = os.path.abspath(os.curdir)
        print(u'Diretório atual ' + curr_dir)
        print(u'Movendo para o diretório de entrada (input) ')
        os.chdir("..")
        os.chdir(os.path.join('input', 'cenarios', input_folder))

        curr_dir = os.path.abspath(os.curdir)
        print(u'Diretório de entrada (input) ' + curr_dir)

    # Infected columns
    # Exposed: de 0*12 a 1*12
    # X: de x*12 a (x+1)*12
    E = 0
    A = 1
    I = 2
    Qi = 3
    Qa = 4
    H = 5

    age_strata = 16

    # Read parameters
    df = pd.read_csv('parameters.csv', sep=',')
    mu_eq = df.values[1, 1:]
    mu_cov = df.values[2, 1:]
    theta = 0
    #alphas = df.values[3, 1:]
    rho = df.values[4, 1:]
    #phi = df.values[5, 1:]
    #etas = df.values[6, 1:]
    a = df.values[7, 1:]
    gamma_H = df.values[8, 1:]
    gamma_HR = df.values[9, 1:]
    gamma_RI = df.values[10, 1:]
    gamma_RA = df.values[11, 1:]
    gamma_RQI = df.values[12, 1:]
    gamma_RQA = df.values[13, 1:]
    gamma_HQI = df.values[14, 1:]

    # Read demography
    df = pd.read_csv('demographic_data.csv', sep=',')
    pop_strata = df.values[0, 1:]
    pop_strata = pop_strata/float(np.sum(pop_strata))

    # Read beta gamma
    df = pd.read_csv('beta_gama.csv', sep=',')
    #xI = df.values[0, 17:17 + age_strata]
    #xA = df.values[0, 17 + age_strata:17 + 2 * age_strata]
    #betas = df.values[0:age_strata, 17 + 2 * age_strata:17 + 3 * age_strata]
    #gamma_QI = np.multiply(xI, gamma_H + gamma_RI) / (1 - xI)
    #gamma_QA = np.multiply(xA, gamma_RA) / (1 - xA)
    betas = contact_matrix
    xI = np.zeros([age_strata, ], dtype=np.float64)
    xA = np.zeros([age_strata, ], dtype=np.float64)
    gamma_QI = np.multiply(xI, gamma_H + gamma_RI) / (1 - xI)
    gamma_QA = np.multiply(xA, gamma_RA) / (1 - xA)
    # Estado inicial
    N = sym.Symbol('N')
    S0 = N

    # A unica entrada de novas infecções é em E, vindo de A, I ou E
    F = np.zeros((6 * age_strata, 6 * age_strata))

    # Some matrix manipulation
    xi = 0.5
    alpha = 0.75
    B_1 = xi * np.multiply(betas, np.tile(pop_strata.reshape(age_strata, 1), (1, age_strata)))
    B_2 = alpha * np.multiply(betas, np.tile(pop_strata.reshape(age_strata, 1), (1, age_strata)))
    B_3 = np.multiply(betas, np.tile(pop_strata.reshape(age_strata, 1), (1, age_strata)))
    F[E * age_strata: (E + 1) * age_strata, E * age_strata: (E + 1) * age_strata] = B_1
    F[E * age_strata: (E + 1) * age_strata, A * age_strata: (A + 1) * age_strata] = B_2
    F[E * age_strata: (E + 1) * age_strata, I * age_strata: (I + 1) * age_strata] = B_3

    # Transições entre compartimentos de infectados
    V = np.zeros((6 * age_strata, 6 * age_strata))
    # Exposto -> Exposto
    for i in range(age_strata):
        j = E + i
        V[j, j] = mu_eq[i] + a[i]
    # Assintomático
    for i in range(age_strata):
        j = A * age_strata + i
        V[j, j] = mu_eq[i] + gamma_RA[i] + gamma_QA[i]
        V[j, E * age_strata + i] = -a[i] * (1 - rho[i])
    # Infectado sintomático
    for i in range(age_strata):
        j = I * age_strata + i
        V[j, j] = mu_eq[i] + theta * mu_cov[i] + gamma_H[i] + gamma_RI[i] + gamma_QI[i]
        V[j, E * age_strata + i] = -a[i] * rho[i]
    # Quarentenado sintomático
    for i in range(age_strata):
        j = Qi * age_strata + i
        V[j, j] = mu_eq[i] + gamma_HQI[i] + gamma_RQI[i]
        V[j, I * age_strata + i] = -gamma_QI[i]
    # Quarentenado assintomático
    for i in range(age_strata):
        j = Qa * age_strata + i
        V[j, j] = mu_eq[i] + gamma_RQA[i]
        V[j, A * age_strata + i] = -gamma_QA[i]
    # Hospitalizado
    for i in range(age_strata):
        j = H * age_strata + i
        V[j, j] = mu_eq[i] + gamma_HR[i] + mu_cov[i]
        V[j, I * age_strata + i] = -gamma_H[i]
        V[j, Qi * age_strata + i] = -gamma_HQI[i]

    # print(V)
    v_1 = np.linalg.inv(V)
    np.set_printoptions(threshold=sys.maxsize)
    np.savetxt("V.csv", V, delimiter=",")
    np.savetxt("F.csv", F, delimiter=",")
    next_gen_matrix = np.matmul(F, v_1)
    w, _ = linalg.eig(next_gen_matrix)
    R0 = np.max(np.abs(w))
    print("R0 obtido: " + str(R0))

    if flag_where:
        print(u'Voltando para o diretório de script')
        os.chdir("..")
        os.chdir("..")
        os.chdir("..")
        os.chdir("scripts")

    return R0
