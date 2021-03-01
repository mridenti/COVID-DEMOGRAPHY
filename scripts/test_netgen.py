# -*- coding: utf-8 -*-
from netgen import compute_r0
import numpy as np
import csv
import os
import sys
import matplotlib.pyplot as plt
import subprocess


def read_output(out_file):
    s = np.zeros([t_days, age_strata], dtype=np.float64)
    e = np.zeros([t_days, age_strata], dtype=np.float64)
    y = np.zeros([t_days, age_strata], dtype=np.float64)
    r = np.zeros([t_days, age_strata], dtype=np.float64)
    n = np.zeros([t_days, age_strata], dtype=np.float64)
    a = np.zeros([t_days, age_strata], dtype=np.float64)
    c = np.zeros([t_days, age_strata], dtype=np.float64)
    h = np.zeros([t_days, age_strata], dtype=np.float64)
    ll = np.zeros([t_days, age_strata], dtype=np.float64)
    ri = np.zeros([t_days, age_strata], dtype=np.float64)
    with open(out_file, "r") as csvfile:
        spamreader = csv.reader(csvfile, delimiter=',')
        next(spamreader, None)
        j = 0
        for row in spamreader:
            for i in range(age_strata):
                s[j, i] = row[compartments * (i + 1)]
                e[j, i] = row[compartments * (i + 1) + 1]
                y[j, i] = row[compartments * (i + 1) + 2]
                r[j, i] = row[compartments * (i + 1) + 3]
                n[j, i] = row[compartments * (i + 1) + 4]
                a[j, i] = row[compartments * (i + 1) + 5]
                c[j, i] = row[compartments * (i + 1) + 6]
                h[j, i] = row[compartments * (i + 1) + 7]
                ll[j, i] = row[compartments * (i + 1) + 8]
            for ii in range(age_strata):
                ri[j, ii] = row[compartments * (age_strata + 1) + ii + 1]
            j = j + 1
        return s, e, y, r, n, a, c, h, ll, ri


def compute_sum(y, c, h, ll, ri):
    # Soma todas as faixas etárias

    #  S = np.sum(s, axis=1)
    #  E = np.sum(e, axis=1)
    I = np.sum(y, axis=1)
    #  R = np.sum(r, axis=1)
    #  N = np.sum(n, axis=1)
    #  A = np.sum(a, axis=1)
    C = np.sum(c, axis=1)
    H = np.sum(h, axis=1)
    L = np.sum(ll, axis=1)
    RI = np.sum(ri, axis=1)
    Ac = C + RI + I + H

    # Criam-se grupos de faixa etária: 1: 0 - 20 / 2: 20 - 55 / 3: 55 ou mais

    I_A = np.sum(y[:, 0:4], axis=1)
    I_B = np.sum(y[:, 4:11], axis=1)
    I_C = np.sum(y[:, 11:16], axis=1)

    RI_A = np.sum(ri[:, 0:4], axis=1)
    RI_B = np.sum(ri[:, 4:11], axis=1)
    RI_C = np.sum(ri[:, 11:16], axis=1)

    H_A = np.sum(h[:, 0:4], axis=1)
    H_B = np.sum(h[:, 4:11], axis=1)
    H_C = np.sum(h[:, 11:16], axis=1)

    L_A = np.sum(ll[:, 0:4], axis=1)
    L_B = np.sum(ll[:, 4:11], axis=1)
    L_C = np.sum(ll[:, 11:16], axis=1)

    C_A = np.sum(c[:, 0:4], axis=1)
    C_B = np.sum(c[:, 4:11], axis=1)
    C_C = np.sum(c[:, 11:16], axis=1)

    Ac_A = C_A + RI_A + I_A + H_A
    Ac_B = C_B + RI_B + I_B + H_B
    Ac_C = C_C + RI_C + I_C + H_C

    return I, H, L, C, Ac, I_A, I_B, I_C, H_A, H_B, H_C, L_A, L_B, L_C, C_A, C_B, C_C, Ac_A, Ac_B, Ac_C


def plot_scenario(y, c, h, ll, ac, hA, hB, hC, lA, lB, lC, cA, cB, cC, country_name, scenario_name):
    t_dates = np.linspace(0, t_days - 1, t_days, dtype=int)
    plt.figure(1)
    _ = plt.subplot(121)
    plt.plot(t_dates, y, '-', label='I=actual infected')
    plt.plot(t_dates, c, '-', label='C=casualties')
    plt.plot(t_dates, h, '-', label='H=actual hospitalized')
    plt.plot(t_dates, ll, '-', label='L=actual ICU occupied')
    # plt.plot(t_dates, ac, '-', label='Ac=accumulated infected')
    plt.xlabel(u'days')
    plt.ylabel(u'individuals')
    plt.xlim([0, last_day])
    plt.legend(loc='upper left', shadow=False, fontsize='small', framealpha=0)

    _ = plt.subplot(122)
    plt.semilogy(t_dates, y, label='I=actual infected')
    plt.semilogy(t_dates, c, label='C=casualties')
    plt.semilogy(t_dates, h, label='H=actual hospitalized')
    plt.semilogy(t_dates, ll, label='L=actual ICU occupied')
    plt.semilogy(t_dates, ac, label='Ac=accumulated infected')
    plt.legend(loc='upper left', shadow=False, fontsize='small', framealpha=0)
    plt.xlim([0, last_day])
    plt.ylim([1, 2.0 * np.max(ac)])
    plt.xlabel(u'days')
    max_H = np.max(h)
    max_L = np.max(ll)
    max_I = np.max(y)
    max_C = np.max(c)
    t_max_I = t_dates[np.where(y == max_I)]
    t_max_L = t_dates[np.where(ll == max_L)]
    textstr = '\n'.join([r'$Max(H)=%.2e$' % (max_H,), r'$Max(L)=%.2e$' % (max_L,), r'$Max(I)=%.2e$' % (max_I,),
                         r'$t(max(I))=%.f$ days' % (t_max_I,),
                         r'$t(max(L))=%.f$ days' % (t_max_L,),
                         r'Casualties $=%.2e$' % (max_C,)])
    props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
    # place a text box in upper left in axes coords
    plt.text(0.68, 0.15, textstr, transform=plt.gca().transAxes, fontsize='small', verticalalignment='center',
             bbox=props)
    plt.suptitle(country_name[7:])

    name_fig = country_name[:] + scenario_name[:] + '_fitting.png'
    save_fig(name_fig)

    plt.figure(2)
    _ = plt.subplot(121)
    plt.plot(t_dates, hA, '-b', label='H: 0 - 20')
    plt.plot(t_dates, hB, '-g', label='H: 20 - 55')
    plt.plot(t_dates, hC, '-r', label='H: 55 or more')
    plt.plot(t_dates, lA, '--b', label='L: 0 - 20')
    plt.plot(t_dates, lB, '--g', label='L: 20 - 55')
    plt.plot(t_dates, lC, '--r', label='L: 55 or more')
    plt.xlim([0, last_day])
    plt.xlabel(u'days')
    plt.ylabel(u'individuals')
    plt.legend(loc='upper left', shadow=False, fontsize='small', framealpha=0)

    _ = plt.subplot(122)
    plt.semilogy(t_dates, hA, '-b', label='H: 0 - 20')
    plt.semilogy(t_dates, hB, '-g', label='H: 20 - 55')
    plt.semilogy(t_dates, hC, '-r', label='H: 55 or more')
    plt.semilogy(t_dates, lA, '--b', label='L: 0 - 20')
    plt.semilogy(t_dates, lB, '--g', label='L: 20 - 55')
    plt.semilogy(t_dates, lC, '--r', label='L: 55 or more')
    plt.semilogy(t_dates, cA, '-.b', label='C: 0 - 20')
    plt.semilogy(t_dates, cB, '-.g', label='C: 20 - 55')
    plt.semilogy(t_dates, cC, '-.r', label='C: 55 or more')
    plt.xlabel('days')
    plt.xlim([0, last_day])
    plt.ylim([1, 2.0 * np.max([hA, hB, hC])])
    plt.legend(loc='upper left', shadow=False, fontsize='small', framealpha=0.5)
    plt.suptitle(country_name[7:])

    plt.close('all')


def save_fig(fig_name):
    os.chdir("..")
    my_path = os.getcwd()
    dirs = os.path.join(my_path, 'figures')
    try:
        os.makedirs(dirs)
    except OSError:
        pass
    plt.savefig(os.path.join(dirs, fig_name), format='png',
                dpi=300, bbox_inches='tight')
    os.chdir("scripts")


input_folder = 'CenarioCH'
cenario_folder = 'CenarioCH'
output_file = 'result_dataCH.csv'
contact_matrix_all_file = 'contact_matrix_all.csv'

# Number of days in simulation
t_days = 400

# Last day of the year
last_day = 360

# Number of compartments in the output file
compartments = 11

# Initial number of infected
ni = 1

age_strata = 16
day_init = 1
day_next_1 = 100
day_next_2 = 300  # this data is arbitrary
day_next_3 = 350  # this data is arbitrary

r0_post = 2.00
r_0 = 2.00
g_e = 0

# Correction factor for interventions - all set to one
g_e0 = 1.0
g_e1 = 1.0
g_e2 = 1.0
g_e3 = 1.0

interv= ['0', '9', '9', '9']

matrix_all = np.zeros([age_strata, age_strata], dtype=np.float64)

curr_dir = os.path.abspath(os.curdir)
print(u'Diretório atual ' + curr_dir)
print(u'Movendo para o diretório de entrada (input) ')
os.chdir("..")
os.chdir(os.path.join('input', 'cenarios', input_folder))
curr_dir = os.path.abspath(os.curdir)
print(u'Diretório de entrada (input) ' + curr_dir)

# Read contact matrices
# All
with open(contact_matrix_all_file, "r") as csvfile:
    spamreader = csv.reader(csvfile, delimiter=',')
    next(spamreader, None)
    j = 0
    for row in spamreader:
        for i in range(age_strata):
            matrix_all[j, i] = row[i]
        j = j + 1

print(u'Voltando para o diretório de script')
os.chdir("..")
os.chdir("..")
os.chdir("..")
os.chdir("scripts")

R0 = compute_r0(input_folder, matrix_all, True)

subprocess.call(
        ['python', 'cenario_generator.py', '-i', cenario_folder, '-d', str(day_init), str(day_next_1),
         str(day_next_2), str(day_next_3), '-m', '3',
         '-I0', str(ni), '-R0', str(r_0), '-Rp', str(r0_post), '-p', str(g_e0), str(g_e1), str(g_e2), str(g_e3),
         '-itv', interv[0], interv[1], interv[2], interv[3], '-pa', '1.0'])

os.chdir("..")
subprocess.call(['bin/csv_to_input', cenario_folder], stdout=open(os.devnull, 'wb'))
subprocess.call(['bin/spatial_covid0d_estrat.exe', 'input/generated-input.txt', '/'.join(['output',
                    output_file]), '3'], stdout=open(os.devnull, 'wb'))
os.chdir("output")
s, e, y, r, n, a, c, h, l, ri = read_output(output_file)
os.chdir("..")
os.chdir("scripts")
I, H, L, C, Ac, I_A, I_B, I_C, H_A, H_B, H_C, L_A, L_B, L_C, \
    C_A, C_B, C_C, Ac_A, Ac_B, Ac_C = compute_sum(y, c, h, l, ri)
plot_scenario(I, C, H, L, Ac, H_A, H_B, H_C, L_A,
              L_B, L_C, C_A, C_B, C_C, 'test_net_gen', 'scenario1')