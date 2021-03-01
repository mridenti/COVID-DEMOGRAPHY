# coding: utf-8
import matplotlib.pyplot as plt
import numpy as np
import csv
import os
import subprocess

# CONSTANTS AND PARAMETERS

# Output file - this should be fixed
output_file = ['result_dataBR.csv', 'result_dataGR.csv', 'result_dataNG.csv', 'result_dataUS.csv']

# Cenario folder
cenario_folder = ['cenarioBR', 'cenarioGR', 'cenarioNG', 'cenarioUS']

# Number of age groups
age_strata = 16

# Number of days in simulation
t_days = 400

# Last day of the year
last_day = 360

# Number of compartments in the output file
compartments = 11

# Initial number of infected
ni = 1

# Population size

# population = [211291345, 83517045, 200963603, 329064917]
population = [10000000, 10000000, 10000000, 10000000]

# Correction factor for interventions - all set to one
g_e0 = 1.0
g_e1 = 1.0
g_e2 = 1.0
g_e3 = 1.0

# Define some parameters which are usually fitted

r0_post = 1.5
r_0 = 2.00
g_e = 0

day_init = 1
day_next_1 = [41, 48, 24, 48]  # to be better defined later
day_next_2 = 60  # this data is arbitrary
day_next_3 = 120  # this data is arbitrary


# FUNCTION DEFINITIONS


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


def run_scenario(index, interv, scenario_tag):
    subprocess.call(
        ['python', 'cenario_generator.py', '-i', cenario_folder[index], '-d', str(day_init), str(day_next_1[index]),
         str(day_next_2), str(day_next_3), '-m', '3',
         '-I0', str(ni), '-R0', str(r_0), '-Rp', str(r0_post), '-p', str(g_e0), str(g_e1), str(g_e2), str(g_e3),
         '-itv', interv[0], interv[1], interv[2], interv[3]], '-pa', '1.0')

    os.chdir("..")
    subprocess.call(['bin/csv_to_input', cenario_folder[index]], stdout=open(os.devnull, 'wb'))
    subprocess.call(['bin/spatial_covid0d_estrat.exe', 'input/generated-input.txt', '/'.join(['output',
                                                                                              output_file[index]]),
                     '3'],
                    stdout=open(os.devnull, 'wb'))
    os.chdir("output")
    sBR, eBR, yBR, rBR, nBR, aBR, cBR, hBR, lBR, riBR = read_output(output_file[index])
    os.chdir("..")
    os.chdir("scripts")
    I_BR, H_BR, L_BR, C_BR, Ac_BR, I_A_BR, I_B_BR, I_C_BR, H_A_BR, H_B_BR, H_C_BR, L_A_BR, L_B_BR, L_C_BR, \
    C_A_BR, C_B_BR, C_C_BR, Ac_A_BR, Ac_B_BR, Ac_C_BR = compute_sum(yBR, cBR, hBR, lBR, riBR)
    plot_scenario(I_BR, C_BR, H_BR, L_BR, Ac_BR, H_A_BR, H_B_BR, H_C_BR, L_A_BR,
                  L_B_BR, L_C_BR, C_A_BR, C_B_BR, C_C_BR, cenario_folder[index], scenario_tag)
    return sBR, eBR, yBR, rBR, nBR, aBR, cBR, hBR, lBR, riBR, I_BR, H_BR, L_BR, C_BR, Ac_BR, I_A_BR, I_B_BR, I_C_BR, \
           H_A_BR, H_B_BR, H_C_BR, L_A_BR, L_B_BR, L_C_BR, C_A_BR, C_B_BR, C_C_BR, Ac_A_BR, Ac_B_BR, Ac_C_BR


# Gráfico de barras - cenário fixo
def plot_same_scenario(an_preval_a, an_preval_b, an_preval_c, an_preval_d,
                       an_preval_a_A, an_preval_b_A, an_preval_c_A, an_preval_d_A,
                       an_preval_a_B, an_preval_b_B, an_preval_c_B, an_preval_d_B,
                       an_preval_a_C, an_preval_b_C, an_preval_c_C, an_preval_d_C,
                       hosp_max_a, hosp_max_b, hosp_max_c, hosp_max_d,
                       hosp_max_a_A, hosp_max_b_A, hosp_max_c_A, hosp_max_d_A,
                       hosp_max_a_B, hosp_max_b_B, hosp_max_c_B, hosp_max_d_B,
                       hosp_max_a_C, hosp_max_b_C, hosp_max_c_C, hosp_max_d_C,
                       fat_rate_a, fat_rate_b, fat_rate_c, fat_rate_d,
                       fat_rate_a_A, fat_rate_b_A, fat_rate_c_A, fat_rate_d_A,
                       fat_rate_a_B, fat_rate_b_B, fat_rate_c_B, fat_rate_d_B,
                       fat_rate_a_C, fat_rate_b_C, fat_rate_c_C, fat_rate_d_C,
                       cenario_string):
    countries = ['a', 'b', 'c', 'd']
    x = np.arange(len(countries))  # the label locations
    width = 0.2

    # plot period prevalence
    _, ax1 = plt.subplots()
    an_preval = [an_preval_a, an_preval_b, an_preval_c, an_preval_d]
    an_preval_A = [an_preval_a_A, an_preval_b_A, an_preval_c_A, an_preval_d_A]
    an_preval_B = [an_preval_a_B, an_preval_b_B, an_preval_c_B, an_preval_d_B]
    an_preval_C = [an_preval_a_C, an_preval_b_C, an_preval_c_C, an_preval_d_C]
    ax1.bar(x - 1.5 * width, an_preval, width, label='All')
    ax1.bar(x - 0.5 * width, an_preval_A, width, label='0-20')
    ax1.bar(x + 0.5 * width, an_preval_B, width, label='20-55')
    ax1.bar(x + 1.5 * width, an_preval_C, width, label='>55')
    ax1.set_xticks(x)
    ax1.set_xticklabels(countries)
    ax1.set_ylabel('Annual prevalence (%)')
    ax1.set_title('Annual prevalence (%) by country - ' + cenario_string)
    ax1.legend()

    name_fig = cenario_string[:] + '_prevalence.png'
    save_fig(name_fig)

    # plot max hosp
    _, ax2 = plt.subplots()
    hosp_max = [hosp_max_a, hosp_max_b, hosp_max_c, hosp_max_d]
    hosp_max_A = [hosp_max_a_A, hosp_max_b_A, hosp_max_c_A, hosp_max_d_A]
    hosp_max_B = [hosp_max_a_B, hosp_max_b_B, hosp_max_c_B, hosp_max_d_B]
    hosp_max_C = [hosp_max_a_C, hosp_max_b_C, hosp_max_c_C, hosp_max_d_C]
    ax2.bar(x - 1.5 * width, hosp_max, width, label='All')
    ax2.bar(x - 0.5 * width, hosp_max_A, width, label='0-20')
    ax2.bar(x + 0.5 * width, hosp_max_B, width, label='20-55')
    ax2.bar(x + 1.5 * width, hosp_max_C, width, label='>55')
    ax2.set_xticks(x)
    ax2.set_xticklabels(countries)
    ax2.set_ylabel('Maximum hospitalization demand')
    ax2.set_title('Maximum hospitalization demand by country - ' + cenario_string)
    ax2.legend()

    name_fig = cenario_string[:] + '_hospitalization.png'
    save_fig(name_fig)

    # plot fatality rate
    _, ax3 = plt.subplots()
    fat_rate = [fat_rate_a, fat_rate_b, fat_rate_c, fat_rate_d]
    fat_rate_A = [fat_rate_a_A, fat_rate_b_A, fat_rate_c_A, fat_rate_d_A]
    fat_rate_B = [fat_rate_a_B, fat_rate_b_B, fat_rate_c_B, fat_rate_d_B]
    fat_rate_C = [fat_rate_a_C, fat_rate_b_C, fat_rate_c_C, fat_rate_d_C]
    ax3.bar(x - 1.5 * width, fat_rate, width, label='All')
    ax3.bar(x - 0.5 * width, fat_rate_A, width, label='0-20')
    ax3.bar(x + 0.5 * width, fat_rate_B, width, label='20-55')
    ax3.bar(x + 1.5 * width, fat_rate_C, width, label='>55')
    ax3.set_xticks(x)
    ax3.set_xticklabels(countries)
    ax3.set_ylabel('Annual fatality rate (deaths/100k)')
    ax3.set_title('Annual fatality rate (deaths/100k) - ' + cenario_string)
    ax3.legend()

    name_fig = cenario_string[:] + '_fatality_per_100k.png'
    save_fig(name_fig)

    plt.show()


# Gráfico de barras - cenário fixo
def save_sheet_same_scenario(an_preval_a, an_preval_b, an_preval_c, an_preval_d,
                             an_preval_a_A, an_preval_b_A, an_preval_c_A, an_preval_d_A,
                             an_preval_a_B, an_preval_b_B, an_preval_c_B, an_preval_d_B,
                             an_preval_a_C, an_preval_b_C, an_preval_c_C, an_preval_d_C,
                             hosp_max_a, hosp_max_b, hosp_max_c, hosp_max_d,
                             hosp_max_a_A, hosp_max_b_A, hosp_max_c_A, hosp_max_d_A,
                             hosp_max_a_B, hosp_max_b_B, hosp_max_c_B, hosp_max_d_B,
                             hosp_max_a_C, hosp_max_b_C, hosp_max_c_C, hosp_max_d_C,
                             fat_rate_a, fat_rate_b, fat_rate_c, fat_rate_d,
                             fat_rate_a_A, fat_rate_b_A, fat_rate_c_A, fat_rate_d_A,
                             fat_rate_a_B, fat_rate_b_B, fat_rate_c_B, fat_rate_d_B,
                             fat_rate_a_C, fat_rate_b_C, fat_rate_c_C, fat_rate_d_C,
                             cenario_string):
    name_sheet = cenario_string[:] + '_datasheet.csv'
    countries = ['a', 'b', 'c', 'd']

    os.chdir("..")
    my_path = os.getcwd()
    dirs = os.path.join(my_path, 'output')

    sheet_path = os.path.join(dirs, name_sheet)

    title1 = ['prevalence_(%)']
    title2 = ['_maximum_demand_for_hospitalization']
    title3 = ['_fatality_per_100k']

    an_preval = [an_preval_a, an_preval_b, an_preval_c, an_preval_d]
    an_preval_A = [an_preval_a_A, an_preval_b_A, an_preval_c_A, an_preval_d_A]
    an_preval_B = [an_preval_a_B, an_preval_b_B, an_preval_c_B, an_preval_d_B]
    an_preval_C = [an_preval_a_C, an_preval_b_C, an_preval_c_C, an_preval_d_C]

    hosp_max = [hosp_max_a, hosp_max_b, hosp_max_c, hosp_max_d]
    hosp_max_A = [hosp_max_a_A, hosp_max_b_A, hosp_max_c_A, hosp_max_d_A]
    hosp_max_B = [hosp_max_a_B, hosp_max_b_B, hosp_max_c_B, hosp_max_d_B]
    hosp_max_C = [hosp_max_a_C, hosp_max_b_C, hosp_max_c_C, hosp_max_d_C]

    fat_rate = [fat_rate_a, fat_rate_b, fat_rate_c, fat_rate_d]
    fat_rate_A = [fat_rate_a_A, fat_rate_b_A, fat_rate_c_A, fat_rate_d_A]
    fat_rate_B = [fat_rate_a_B, fat_rate_b_B, fat_rate_c_B, fat_rate_d_B]
    fat_rate_C = [fat_rate_a_C, fat_rate_b_C, fat_rate_c_C, fat_rate_d_C]

    with open(sheet_path, 'wb') as csvfile:
        spamwriter = csv.writer(csvfile)
        spamwriter.writerow(title1)
        spamwriter.writerow(np.concatenate((['Country'], countries)))
        spamwriter.writerow(np.concatenate((['All'], an_preval)))
        spamwriter.writerow(np.concatenate((['0-20'], an_preval_A)))
        spamwriter.writerow(np.concatenate((['20-55'], an_preval_B)))
        spamwriter.writerow(np.concatenate((['55 or more'], an_preval_C)))
        spamwriter.writerow([''])
        spamwriter.writerow(title2)
        spamwriter.writerow(np.concatenate((['Country'], countries)))
        spamwriter.writerow(np.concatenate((['All'], hosp_max)))
        spamwriter.writerow(np.concatenate((['0-20'], hosp_max_A)))
        spamwriter.writerow(np.concatenate((['20-55'], hosp_max_B)))
        spamwriter.writerow(np.concatenate((['55 or more'], hosp_max_C)))
        spamwriter.writerow([''])
        spamwriter.writerow(title3)
        spamwriter.writerow(np.concatenate((['Country'], countries)))
        spamwriter.writerow(np.concatenate((['All'], fat_rate)))
        spamwriter.writerow(np.concatenate((['0-20'], fat_rate_A)))
        spamwriter.writerow(np.concatenate((['20-55'], fat_rate_B)))
        spamwriter.writerow(np.concatenate((['55 or more'], fat_rate_C)))

    os.chdir("scripts")


def plot_same_country(an_preval_1, an_preval_A_1, an_preval_B_1, an_preval_C_1,
                      an_preval_2, an_preval_A_2, an_preval_B_2, an_preval_C_2,
                      an_preval_3, an_preval_A_3, an_preval_B_3, an_preval_C_3,
                      an_preval_4, an_preval_A_4, an_preval_B_4, an_preval_C_4,
                      an_preval_5, an_preval_A_5, an_preval_B_5, an_preval_C_5,
                      an_preval_6, an_preval_A_6, an_preval_B_6, an_preval_C_6,
                      hosp_max_1, hosp_max_A_1, hosp_max_B_1, hosp_max_C_1,
                      hosp_max_2, hosp_max_A_2, hosp_max_B_2, hosp_max_C_2,
                      hosp_max_3, hosp_max_A_3, hosp_max_B_3, hosp_max_C_3,
                      hosp_max_4, hosp_max_A_4, hosp_max_B_4, hosp_max_C_4,
                      hosp_max_5, hosp_max_A_5, hosp_max_B_5, hosp_max_C_5,
                      hosp_max_6, hosp_max_A_6, hosp_max_B_6, hosp_max_C_6,
                      fat_rate_1, fat_rate_A_1, fat_rate_B_1, fat_rate_C_1,
                      fat_rate_2, fat_rate_A_2, fat_rate_B_2, fat_rate_C_2,
                      fat_rate_3, fat_rate_A_3, fat_rate_B_3, fat_rate_C_3,
                      fat_rate_4, fat_rate_A_4, fat_rate_B_4, fat_rate_C_4,
                      fat_rate_5, fat_rate_A_5, fat_rate_B_5, fat_rate_C_5,
                      fat_rate_6, fat_rate_A_6, fat_rate_B_6, fat_rate_C_6,
                      country_tag):
    scenarios = ['1', '2', '3', '4', '5', '6']
    x = np.arange(len(scenarios))  # the label locations
    width = 0.20

    # plot period prevalence
    _, ax1 = plt.subplots()
    an_preval = [an_preval_1, an_preval_2, an_preval_3, an_preval_4, an_preval_5, an_preval_6]
    an_preval_A = [an_preval_A_1, an_preval_A_2, an_preval_A_3, an_preval_A_4, an_preval_A_5, an_preval_A_6]
    an_preval_B = [an_preval_B_1, an_preval_B_2, an_preval_B_3, an_preval_B_4, an_preval_B_5, an_preval_B_6]
    an_preval_C = [an_preval_C_1, an_preval_C_2, an_preval_C_3, an_preval_C_4, an_preval_C_5, an_preval_C_6]
    ax1.bar(x - 1.5 * width, an_preval, width, label='All')
    ax1.bar(x - 0.5 * width, an_preval_A, width, label='0-20')
    ax1.bar(x + 0.5 * width, an_preval_B, width, label='20-55')
    ax1.bar(x + 1.5 * width, an_preval_C, width, label='>55')
    ax1.set_xticks(x)
    ax1.set_xticklabels(scenarios)
    ax1.set_ylabel('Annual prevalence (%)')
    ax1.set_title('Annual prevalence (%) for country - ' + country_tag)
    ax1.legend()

    name_fig = country_tag[:] + '_prevalence.png'
    save_fig(name_fig)

    # plot max hosp
    _, ax2 = plt.subplots()
    hosp_max = [hosp_max_1, hosp_max_2, hosp_max_3, hosp_max_4, hosp_max_5, hosp_max_6]
    hosp_max_A = [hosp_max_A_1, hosp_max_A_2, hosp_max_A_3, hosp_max_A_4, hosp_max_A_5, hosp_max_A_6]
    hosp_max_B = [hosp_max_B_1, hosp_max_B_2, hosp_max_B_3, hosp_max_B_4, hosp_max_B_5, hosp_max_B_6]
    hosp_max_C = [hosp_max_C_1, hosp_max_C_2, hosp_max_C_3, hosp_max_C_4, hosp_max_C_5, hosp_max_C_6]
    ax2.bar(x - 1.5 * width, hosp_max, width, label='All')
    ax2.bar(x - 0.5 * width, hosp_max_A, width, label='0-20')
    ax2.bar(x + 0.5 * width, hosp_max_B, width, label='20-55')
    ax2.bar(x + 1.5 * width, hosp_max_C, width, label='>55')
    ax2.set_xticks(x)
    ax2.set_xticklabels(scenarios)
    ax2.set_ylabel('Maximum hospitalization demand')
    ax2.set_title('Maximum hospitalization demand for country - ' + country_tag)
    ax2.legend()

    name_fig = country_tag[:] + '_hospitalization.png'
    save_fig(name_fig)

    # plot fatality rate
    _, ax3 = plt.subplots()
    fat_rate = [fat_rate_1, fat_rate_2, fat_rate_3, fat_rate_4, fat_rate_5, fat_rate_6]
    fat_rate_A = [fat_rate_A_1, fat_rate_A_2, fat_rate_A_3, fat_rate_A_4, fat_rate_A_5, fat_rate_A_6]
    fat_rate_B = [fat_rate_B_1, fat_rate_B_2, fat_rate_B_3, fat_rate_B_4, fat_rate_B_5, fat_rate_B_6]
    fat_rate_C = [fat_rate_C_1, fat_rate_C_2, fat_rate_C_3, fat_rate_C_4, fat_rate_C_5, fat_rate_C_6]
    ax3.bar(x - 1.5 * width, fat_rate, width, label='All')
    ax3.bar(x - 0.5 * width, fat_rate_A, width, label='0-20')
    ax3.bar(x + 0.5 * width, fat_rate_B, width, label='20-55')
    ax3.bar(x + 1.5 * width, fat_rate_C, width, label='>55')
    ax3.set_xticks(x)
    ax3.set_xticklabels(scenarios)
    ax3.set_ylabel('Annual fatality rate (deaths/100k)')
    ax3.set_title('Annual fatality rate (deaths/100k) - ' + country_tag)
    ax3.legend()

    name_fig = country_tag[:] + '_fatality_per_100k.png'
    save_fig(name_fig)

    plt.show()


def save_sheet_same_country(an_preval_1, an_preval_A_1, an_preval_B_1, an_preval_C_1,
                            an_preval_2, an_preval_A_2, an_preval_B_2, an_preval_C_2,
                            an_preval_3, an_preval_A_3, an_preval_B_3, an_preval_C_3,
                            an_preval_4, an_preval_A_4, an_preval_B_4, an_preval_C_4,
                            an_preval_5, an_preval_A_5, an_preval_B_5, an_preval_C_5,
                            an_preval_6, an_preval_A_6, an_preval_B_6, an_preval_C_6,
                            hosp_max_1, hosp_max_A_1, hosp_max_B_1, hosp_max_C_1,
                            hosp_max_2, hosp_max_A_2, hosp_max_B_2, hosp_max_C_2,
                            hosp_max_3, hosp_max_A_3, hosp_max_B_3, hosp_max_C_3,
                            hosp_max_4, hosp_max_A_4, hosp_max_B_4, hosp_max_C_4,
                            hosp_max_5, hosp_max_A_5, hosp_max_B_5, hosp_max_C_5,
                            hosp_max_6, hosp_max_A_6, hosp_max_B_6, hosp_max_C_6,
                            fat_rate_1, fat_rate_A_1, fat_rate_B_1, fat_rate_C_1,
                            fat_rate_2, fat_rate_A_2, fat_rate_B_2, fat_rate_C_2,
                            fat_rate_3, fat_rate_A_3, fat_rate_B_3, fat_rate_C_3,
                            fat_rate_4, fat_rate_A_4, fat_rate_B_4, fat_rate_C_4,
                            fat_rate_5, fat_rate_A_5, fat_rate_B_5, fat_rate_C_5,
                            fat_rate_6, fat_rate_A_6, fat_rate_B_6, fat_rate_C_6,
                            country_tag):

    name_sheet = country_tag[:] + '_datasheet.csv'
    scenarios = ['1', '2', '3', '4', '5', '6']

    os.chdir("..")
    my_path = os.getcwd()
    dirs = os.path.join(my_path, 'output')

    sheet_path = os.path.join(dirs, name_sheet)

    title1 = ['prevalence_(%)']
    title2 = ['_maximum_demand_for_hospitalization']
    title3 = ['_fatality_per_100k']

    an_preval = [an_preval_1, an_preval_2, an_preval_3, an_preval_4, an_preval_5, an_preval_6]
    an_preval_A = [an_preval_A_1, an_preval_A_2, an_preval_A_3, an_preval_A_4, an_preval_A_5, an_preval_A_6]
    an_preval_B = [an_preval_B_1, an_preval_B_2, an_preval_B_3, an_preval_B_4, an_preval_B_5, an_preval_B_6]
    an_preval_C = [an_preval_C_1, an_preval_C_2, an_preval_C_3, an_preval_C_4, an_preval_C_5, an_preval_C_6]

    hosp_max = [hosp_max_1, hosp_max_2, hosp_max_3, hosp_max_4, hosp_max_5, hosp_max_6]
    hosp_max_A = [hosp_max_A_1, hosp_max_A_2, hosp_max_A_3, hosp_max_A_4, hosp_max_A_5, hosp_max_A_6]
    hosp_max_B = [hosp_max_B_1, hosp_max_B_2, hosp_max_B_3, hosp_max_B_4, hosp_max_B_5, hosp_max_B_6]
    hosp_max_C = [hosp_max_C_1, hosp_max_C_2, hosp_max_C_3, hosp_max_C_4, hosp_max_C_5, hosp_max_C_6]

    fat_rate = [fat_rate_1, fat_rate_2, fat_rate_3, fat_rate_4, fat_rate_5, fat_rate_6]
    fat_rate_A = [fat_rate_A_1, fat_rate_A_2, fat_rate_A_3, fat_rate_A_4, fat_rate_A_5, fat_rate_A_6]
    fat_rate_B = [fat_rate_B_1, fat_rate_B_2, fat_rate_B_3, fat_rate_B_4, fat_rate_B_5, fat_rate_B_6]
    fat_rate_C = [fat_rate_C_1, fat_rate_C_2, fat_rate_C_3, fat_rate_C_4, fat_rate_C_5, fat_rate_C_6]

    with open(sheet_path, 'wb') as csvfile:
        spamwriter = csv.writer(csvfile)
        spamwriter.writerow(title1)
        spamwriter.writerow(np.concatenate((['Scenarios'], scenarios)))
        spamwriter.writerow(np.concatenate((['All'], an_preval)))
        spamwriter.writerow(np.concatenate((['0-20'], an_preval_A)))
        spamwriter.writerow(np.concatenate((['20-55'], an_preval_B)))
        spamwriter.writerow(np.concatenate((['55 or more'], an_preval_C)))
        spamwriter.writerow([''])
        spamwriter.writerow(title2)
        spamwriter.writerow(np.concatenate((['Scenarios'], scenarios)))
        spamwriter.writerow(np.concatenate((['All'], hosp_max)))
        spamwriter.writerow(np.concatenate((['0-20'], hosp_max_A)))
        spamwriter.writerow(np.concatenate((['20-55'], hosp_max_B)))
        spamwriter.writerow(np.concatenate((['55 or more'], hosp_max_C)))
        spamwriter.writerow([''])
        spamwriter.writerow(title3)
        spamwriter.writerow(np.concatenate((['Scenarios'], scenarios)))
        spamwriter.writerow(np.concatenate((['All'], fat_rate)))
        spamwriter.writerow(np.concatenate((['0-20'], fat_rate_A)))
        spamwriter.writerow(np.concatenate((['20-55'], fat_rate_B)))
        spamwriter.writerow(np.concatenate((['55 or more'], fat_rate_C)))

    os.chdir("scripts")

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


# MAIN CALLS AND RUNS

intervention_string_1 = ['0', '1', '1', '1']  # Scenario 1 - no intervention
intervention_string_2 = ['0', '9', '9', '9']  # Scenario 2 - Fech. Esc./Isol. Idos./Dist. Soc./Trabal.
intervention_string_3 = ['0', '8', '8', '8']  # Scenario 3 - Fech. Esc./Isol. Idos./Dist. Soc./
intervention_string_4 = ['0', '3', '3', '3']  # Scenario 4 - Fech. Esc./Dist. Soc.
intervention_string_5 = ['0', '13', '13', '13']  # Scenario 5 - Isol. Idos./Dist. Soc.
intervention_string_6 = ['0', '12', '12', '12']  # Scenario 6 - Dist. Soc.

# BR simulation
index = 0

# Scenario 1 - no intervention
sBR1, eBR1, yBR1, rBR1, nBR1, aBR1, cBR1, hBR1, lBR1, riBR1, I1_BR, H1_BR, L1_BR, C1_BR, Ac1_BR, I1_A_BR, I1_B_BR, \
I1_C_BR, H1_A_BR, H1_B_BR, H1_C_BR, L1_A_BR, L1_B_BR, L1_C_BR, C1_A_BR, C1_B_BR, C1_C_BR, Ac1_A_BR, Ac1_B_BR, \
Ac1_C_BR = run_scenario(index, intervention_string_1, 'scenario1')

# Scenario 2 - Fech. Esc./Isol. Idos./Dist. Soc./Trabal.
sBR2, eBR2, yBR2, rBR2, nBR2, aBR2, cBR2, hBR2, lBR2, riBR2, I2_BR, H2_BR, L2_BR, C2_BR, Ac2_BR, I2_A_BR, I2_B_BR, \
I2_C_BR, H2_A_BR, H2_B_BR, H2_C_BR, L2_A_BR, L2_B_BR, L2_C_BR, C2_A_BR, C2_B_BR, C2_C_BR, Ac2_A_BR, Ac2_B_BR, \
Ac2_C_BR = run_scenario(index, intervention_string_2, 'scenario2')

# Scenario 3 - Fech. Esc./Isol. Idos./Dist. Soc./
sBR3, eBR3, yBR3, rBR3, nBR3, aBR3, cBR3, hBR3, lBR3, riBR3, I3_BR, H3_BR, L3_BR, C3_BR, Ac3_BR, I3_A_BR, I3_B_BR, \
I3_C_BR, H3_A_BR, H3_B_BR, H3_C_BR, L3_A_BR, L3_B_BR, L3_C_BR, C3_A_BR, C3_B_BR, C3_C_BR, Ac3_A_BR, Ac3_B_BR, \
Ac3_C_BR = run_scenario(index, intervention_string_3, 'scenario3')

# Scenario 4 - Fech. Esc./Dist. Soc.
sBR4, eBR4, yBR4, rBR4, nBR4, aBR4, cBR4, hBR4, lBR4, riBR4, I4_BR, H4_BR, L4_BR, C4_BR, Ac4_BR, I4_A_BR, I4_B_BR, \
I4_C_BR, H4_A_BR, H4_B_BR, H4_C_BR, L4_A_BR, L4_B_BR, L4_C_BR, C4_A_BR, C4_B_BR, C4_C_BR, Ac4_A_BR, Ac4_B_BR, \
Ac4_C_BR = run_scenario(index, intervention_string_4, 'scenario4')

# Scenario 5 - Iso. Idos./Dist. Soc.
sBR5, eBR5, yBR5, rBR5, nBR5, aBR5, cBR5, hBR5, lBR5, riBR5, I5_BR, H5_BR, L5_BR, C5_BR, Ac5_BR, I5_A_BR, I5_B_BR, \
I5_C_BR, H5_A_BR, H5_B_BR, H5_C_BR, L5_A_BR, L5_B_BR, L5_C_BR, C5_A_BR, C5_B_BR, C5_C_BR, Ac5_A_BR, Ac5_B_BR, \
Ac5_C_BR = run_scenario(index, intervention_string_5, 'scenario5')

# Scenario 6 - Dist. Soc.
sBR6, eBR6, yBR6, rBR6, nBR6, aBR6, cBR6, hBR6, lBR6, riBR6, I6_BR, H6_BR, L6_BR, C6_BR, Ac6_BR, I6_A_BR, I6_B_BR, \
I6_C_BR, H6_A_BR, H6_B_BR, H6_C_BR, L6_A_BR, L6_B_BR, L6_C_BR, C6_A_BR, C6_B_BR, C6_C_BR, Ac6_A_BR, Ac6_B_BR, \
Ac6_C_BR = run_scenario(index, intervention_string_6, 'scenario6')

# GR simulation
index = 1

# Scenario 1 - no intervention
sGR1, eGR1, yGR1, rGR1, nGR1, aGR1, cGR1, hGR1, lGR1, riGR1, I1_GR, H1_GR, L1_GR, C1_GR, Ac1_GR, I1_A_GR, I1_B_GR, \
I1_C_GR, H1_A_GR, H1_B_GR, H1_C_GR, L1_A_GR, L1_B_GR, L1_C_GR, C1_A_GR, C1_B_GR, C1_C_GR, Ac1_A_GR, Ac1_B_GR, \
Ac1_C_GR = run_scenario(index, intervention_string_1, 'scenario1')

# Scenario 2 - Fech. Esc./Isol. Idos./Dist. Soc./Trabal.
sGR2, eGR2, yGR2, rGR2, nGR2, aGR2, cGR2, hGR2, lGR2, riGR2, I2_GR, H2_GR, L2_GR, C2_GR, Ac2_GR, I2_A_GR, I2_B_GR, \
I2_C_GR, H2_A_GR, H2_B_GR, H2_C_GR, L2_A_GR, L2_B_GR, L2_C_GR, C2_A_GR, C2_B_GR, C2_C_GR, Ac2_A_GR, Ac2_B_GR, \
Ac2_C_GR = run_scenario(index, intervention_string_2, 'scenario2')

# Scenario 3 - Fech. Esc./Isol. Idos./Dist. Soc./
sGR3, eGR3, yGR3, rGR3, nGR3, aGR3, cGR3, hGR3, lGR3, riGR3, I3_GR, H3_GR, L3_GR, C3_GR, Ac3_GR, I3_A_GR, I3_B_GR, \
I3_C_GR, H3_A_GR, H3_B_GR, H3_C_GR, L3_A_GR, L3_B_GR, L3_C_GR, C3_A_GR, C3_B_GR, C3_C_GR, Ac3_A_GR, Ac3_B_GR, \
Ac3_C_GR = run_scenario(index, intervention_string_3, 'scenario3')

# Scenario 4 - Fech. Esc./Dist. Soc.
sGR4, eGR4, yGR4, rGR4, nGR4, aGR4, cGR4, hGR4, lGR4, riGR4, I4_GR, H4_GR, L4_GR, C4_GR, Ac4_GR, I4_A_GR, I4_B_GR, \
I4_C_GR, H4_A_GR, H4_B_GR, H4_C_GR, L4_A_GR, L4_B_GR, L4_C_GR, C4_A_GR, C4_B_GR, C4_C_GR, Ac4_A_GR, Ac4_B_GR, \
Ac4_C_GR = run_scenario(index, intervention_string_4, 'scenario4')

# Scenario 5 - Iso. Idos./Dist. Soc.
sGR5, eGR5, yGR5, rGR5, nGR5, aGR5, cGR5, hGR5, lGR5, riGR5, I5_GR, H5_GR, L5_GR, C5_GR, Ac5_GR, I5_A_GR, I5_B_GR, \
I5_C_GR, H5_A_GR, H5_B_GR, H5_C_GR, L5_A_GR, L5_B_GR, L5_C_GR, C5_A_GR, C5_B_GR, C5_C_GR, Ac5_A_GR, Ac5_B_GR, \
Ac5_C_GR = run_scenario(index, intervention_string_5, 'scenario5')

# Scenario 6 - Dist. Soc.
sGR6, eGR6, yGR6, rGR6, nGR6, aGR6, cGR6, hGR6, lGR6, riGR6, I6_GR, H6_GR, L6_GR, C6_GR, Ac6_GR, I6_A_GR, I6_B_GR, \
I6_C_GR, H6_A_GR, H6_B_GR, H6_C_GR, L6_A_GR, L6_B_GR, L6_C_GR, C6_A_GR, C6_B_GR, C6_C_GR, Ac6_A_GR, Ac6_B_GR, \
Ac6_C_GR = run_scenario(index, intervention_string_6, 'scenario6')

# NG simulation
index = 2

# Scenario 1 - no intervention
sNG1, eNG1, yNG1, rNG1, nNG1, aNG1, cNG1, hNG1, lNG1, riNG1, I1_NG, H1_NG, L1_NG, C1_NG, Ac1_NG, I1_A_NG, I1_B_NG, \
I1_C_NG, H1_A_NG, H1_B_NG, H1_C_NG, L1_A_NG, L1_B_NG, L1_C_NG, C1_A_NG, C1_B_NG, C1_C_NG, Ac1_A_NG, Ac1_B_NG, \
Ac1_C_NG = run_scenario(index, intervention_string_1, 'scenario1')

# Scenario 2 - Fech. Esc./Isol. Idos./Dist. Soc./Trabal.
sNG2, eNG2, yNG2, rNG2, nNG2, aNG2, cNG2, hNG2, lNG2, riNG2, I2_NG, H2_NG, L2_NG, C2_NG, Ac2_NG, I2_A_NG, I2_B_NG, \
I2_C_NG, H2_A_NG, H2_B_NG, H2_C_NG, L2_A_NG, L2_B_NG, L2_C_NG, C2_A_NG, C2_B_NG, C2_C_NG, Ac2_A_NG, Ac2_B_NG, \
Ac2_C_NG = run_scenario(index, intervention_string_2, 'scenario2')

# Scenario 3 - Fech. Esc./Isol. Idos./Dist. Soc./
sNG3, eNG3, yNG3, rNG3, nNG3, aNG3, cNG3, hNG3, lNG3, riNG3, I3_NG, H3_NG, L3_NG, C3_NG, Ac3_NG, I3_A_NG, I3_B_NG, \
I3_C_NG, H3_A_NG, H3_B_NG, H3_C_NG, L3_A_NG, L3_B_NG, L3_C_NG, C3_A_NG, C3_B_NG, C3_C_NG, Ac3_A_NG, Ac3_B_NG, \
Ac3_C_NG = run_scenario(index, intervention_string_3, 'scenario3')

# Scenario 4 - Fech. Esc./Dist. Soc.
sNG4, eNG4, yNG4, rNG4, nNG4, aNG4, cNG4, hNG4, lNG4, riNG4, I4_NG, H4_NG, L4_NG, C4_NG, Ac4_NG, I4_A_NG, I4_B_NG, \
I4_C_NG, H4_A_NG, H4_B_NG, H4_C_NG, L4_A_NG, L4_B_NG, L4_C_NG, C4_A_NG, C4_B_NG, C4_C_NG, Ac4_A_NG, Ac4_B_NG, \
Ac4_C_NG = run_scenario(index, intervention_string_4, 'scenario4')

# Scenario 5 - Iso. Idos./Dist. Soc.
sNG5, eNG5, yNG5, rNG5, nNG5, aNG5, cNG5, hNG5, lNG5, riNG5, I5_NG, H5_NG, L5_NG, C5_NG, Ac5_NG, I5_A_NG, I5_B_NG, \
I5_C_NG, H5_A_NG, H5_B_NG, H5_C_NG, L5_A_NG, L5_B_NG, L5_C_NG, C5_A_NG, C5_B_NG, C5_C_NG, Ac5_A_NG, Ac5_B_NG, \
Ac5_C_NG = run_scenario(index, intervention_string_5, 'scenario5')

# Scenario 6 - Dist. Soc.
sNG6, eNG6, yNG6, rNG6, nNG6, aNG6, cNG6, hNG6, lNG6, riNG6, I6_NG, H6_NG, L6_NG, C6_NG, Ac6_NG, I6_A_NG, I6_B_NG, \
I6_C_NG, H6_A_NG, H6_B_NG, H6_C_NG, L6_A_NG, L6_B_NG, L6_C_NG, C6_A_NG, C6_B_NG, C6_C_NG, Ac6_A_NG, Ac6_B_NG, \
Ac6_C_NG = run_scenario(index, intervention_string_6, 'scenario6')

# US simulation
index = 3

# Scenario 1 - no intervention
sUS1, eUS1, yUS1, rUS1, nUS1, aUS1, cUS1, hUS1, lUS1, riUS1, I1_US, H1_US, L1_US, C1_US, Ac1_US, I1_A_US, I1_B_US, \
I1_C_US, H1_A_US, H1_B_US, H1_C_US, L1_A_US, L1_B_US, L1_C_US, C1_A_US, C1_B_US, C1_C_US, Ac1_A_US, Ac1_B_US, \
Ac1_C_US = run_scenario(index, intervention_string_1, 'scenario1')

# Scenario 2 - Fech. Esc./Isol. Idos./Dist. Soc./Trabal.
sUS2, eUS2, yUS2, rUS2, nUS2, aUS2, cUS2, hUS2, lUS2, riUS2, I2_US, H2_US, L2_US, C2_US, Ac2_US, I2_A_US, I2_B_US, \
I2_C_US, H2_A_US, H2_B_US, H2_C_US, L2_A_US, L2_B_US, L2_C_US, C2_A_US, C2_B_US, C2_C_US, Ac2_A_US, Ac2_B_US, \
Ac2_C_US = run_scenario(index, intervention_string_2, 'scenario2')

# Scenario 3 - Fech. Esc./Isol. Idos./Dist. Soc./
sUS3, eUS3, yUS3, rUS3, nUS3, aUS3, cUS3, hUS3, lUS3, riUS3, I3_US, H3_US, L3_US, C3_US, Ac3_US, I3_A_US, I3_B_US, \
I3_C_US, H3_A_US, H3_B_US, H3_C_US, L3_A_US, L3_B_US, L3_C_US, C3_A_US, C3_B_US, C3_C_US, Ac3_A_US, Ac3_B_US, \
Ac3_C_US = run_scenario(index, intervention_string_3, 'scenario3')

# Scenario 4 - Fech. Esc./Dist. Soc.
sUS4, eUS4, yUS4, rUS4, nUS4, aUS4, cUS4, hUS4, lUS4, riUS4, I4_US, H4_US, L4_US, C4_US, Ac4_US, I4_A_US, I4_B_US, \
I4_C_US, H4_A_US, H4_B_US, H4_C_US, L4_A_US, L4_B_US, L4_C_US, C4_A_US, C4_B_US, C4_C_US, Ac4_A_US, Ac4_B_US, \
Ac4_C_US = run_scenario(index, intervention_string_4, 'scenario4')

# Scenario 5 - Iso. Idos./Dist. Soc.
sUS5, eUS5, yUS5, rUS5, nUS5, aUS5, cUS5, hUS5, lUS5, riUS5, I5_US, H5_US, L5_US, C5_US, Ac5_US, I5_A_US, I5_B_US, \
I5_C_US, H5_A_US, H5_B_US, H5_C_US, L5_A_US, L5_B_US, L5_C_US, C5_A_US, C5_B_US, C5_C_US, Ac5_A_US, Ac5_B_US, \
Ac5_C_US = run_scenario(index, intervention_string_5, 'scenario5')

# Scenario 6 - Dist. Soc.
sUS6, eUS6, yUS6, rUS6, nUS6, aUS6, cUS6, hUS6, lUS6, riUS6, I6_US, H6_US, L6_US, C6_US, Ac6_US, I6_A_US, I6_B_US, \
I6_C_US, H6_A_US, H6_B_US, H6_C_US, L6_A_US, L6_B_US, L6_C_US, C6_A_US, C6_B_US, C6_C_US, Ac6_A_US, Ac6_B_US, \
Ac6_C_US = run_scenario(index, intervention_string_6, 'scenario6')

# Prevalência anual (%) - an_preval

# Demanda máxima de hospitalizações (número absoluto) - hosp_max

# Letalidade (fatality rate), por milhão de habitantes - fat_rate

# Cenário 1

an_preval_BR1 = 100 * Ac1_BR[last_day] / population[0]
an_preval_GR1 = 100 * Ac1_GR[last_day] / population[1]
an_preval_NG1 = 100 * Ac1_NG[last_day] / population[2]
an_preval_US1 = 100 * Ac1_US[last_day] / population[3]

an_preval_A_BR1 = 100 * Ac1_A_BR[last_day] / population[0]
an_preval_A_GR1 = 100 * Ac1_A_GR[last_day] / population[1]
an_preval_A_NG1 = 100 * Ac1_A_NG[last_day] / population[2]
an_preval_A_US1 = 100 * Ac1_A_US[last_day] / population[3]

an_preval_B_BR1 = 100 * Ac1_B_BR[last_day] / population[0]
an_preval_B_GR1 = 100 * Ac1_B_GR[last_day] / population[1]
an_preval_B_NG1 = 100 * Ac1_B_NG[last_day] / population[2]
an_preval_B_US1 = 100 * Ac1_B_US[last_day] / population[3]

an_preval_C_BR1 = 100 * Ac1_C_BR[last_day] / population[0]
an_preval_C_GR1 = 100 * Ac1_C_GR[last_day] / population[1]
an_preval_C_NG1 = 100 * Ac1_C_NG[last_day] / population[2]
an_preval_C_US1 = 100 * Ac1_C_US[last_day] / population[3]

hosp_max_BR1 = np.max(H1_BR)
hosp_max_GR1 = np.max(H1_GR)
hosp_max_NG1 = np.max(H1_NG)
hosp_max_US1 = np.max(H1_US)

hosp_max_A_BR1 = np.max(H1_A_BR)
hosp_max_A_GR1 = np.max(H1_A_GR)
hosp_max_A_NG1 = np.max(H1_A_NG)
hosp_max_A_US1 = np.max(H1_A_US)

hosp_max_B_BR1 = np.max(H1_B_BR)
hosp_max_B_GR1 = np.max(H1_B_GR)
hosp_max_B_NG1 = np.max(H1_B_NG)
hosp_max_B_US1 = np.max(H1_B_US)

hosp_max_C_BR1 = np.max(H1_C_BR)
hosp_max_C_GR1 = np.max(H1_C_GR)
hosp_max_C_NG1 = np.max(H1_C_NG)
hosp_max_C_US1 = np.max(H1_C_US)

fat_rate_BR1 = 1E5 * C1_BR[last_day] / population[0]
fat_rate_GR1 = 1E5 * C1_GR[last_day] / population[1]
fat_rate_NG1 = 1E5 * C1_NG[last_day] / population[2]
fat_rate_US1 = 1E5 * C1_US[last_day] / population[3]

fat_rate_A_BR1 = 1E5 * C1_A_BR[last_day] / population[0]
fat_rate_A_GR1 = 1E5 * C1_A_GR[last_day] / population[1]
fat_rate_A_NG1 = 1E5 * C1_A_NG[last_day] / population[2]
fat_rate_A_US1 = 1E5 * C1_A_US[last_day] / population[3]

fat_rate_B_BR1 = 1E5 * C1_B_BR[last_day] / population[0]
fat_rate_B_GR1 = 1E5 * C1_B_GR[last_day] / population[1]
fat_rate_B_NG1 = 1E5 * C1_B_NG[last_day] / population[2]
fat_rate_B_US1 = 1E5 * C1_B_US[last_day] / population[3]

fat_rate_C_BR1 = 1E5 * C1_C_BR[last_day] / population[0]
fat_rate_C_GR1 = 1E5 * C1_C_GR[last_day] / population[1]
fat_rate_C_NG1 = 1E5 * C1_C_NG[last_day] / population[2]
fat_rate_C_US1 = 1E5 * C1_C_US[last_day] / population[3]

# Cenário 2

an_preval_BR2 = 100 * Ac2_BR[last_day] / population[0]
an_preval_GR2 = 100 * Ac2_GR[last_day] / population[1]
an_preval_NG2 = 100 * Ac2_NG[last_day] / population[2]
an_preval_US2 = 100 * Ac2_US[last_day] / population[3]

an_preval_A_BR2 = 100 * Ac2_A_BR[last_day] / population[0]
an_preval_A_GR2 = 100 * Ac2_A_GR[last_day] / population[1]
an_preval_A_NG2 = 100 * Ac2_A_NG[last_day] / population[2]
an_preval_A_US2 = 100 * Ac2_A_US[last_day] / population[3]

an_preval_B_BR2 = 100 * Ac2_B_BR[last_day] / population[0]
an_preval_B_GR2 = 100 * Ac2_B_GR[last_day] / population[1]
an_preval_B_NG2 = 100 * Ac2_B_NG[last_day] / population[2]
an_preval_B_US2 = 100 * Ac2_B_US[last_day] / population[3]

an_preval_C_BR2 = 100 * Ac2_C_BR[last_day] / population[0]
an_preval_C_GR2 = 100 * Ac2_C_GR[last_day] / population[1]
an_preval_C_NG2 = 100 * Ac2_C_NG[last_day] / population[2]
an_preval_C_US2 = 100 * Ac2_C_US[last_day] / population[3]

hosp_max_BR2 = np.max(H2_BR)
hosp_max_GR2 = np.max(H2_GR)
hosp_max_NG2 = np.max(H2_NG)
hosp_max_US2 = np.max(H2_US)

hosp_max_A_BR2 = np.max(H2_A_BR)
hosp_max_A_GR2 = np.max(H2_A_GR)
hosp_max_A_NG2 = np.max(H2_A_NG)
hosp_max_A_US2 = np.max(H2_A_US)

hosp_max_B_BR2 = np.max(H2_B_BR)
hosp_max_B_GR2 = np.max(H2_B_GR)
hosp_max_B_NG2 = np.max(H2_B_NG)
hosp_max_B_US2 = np.max(H2_B_US)

hosp_max_C_BR2 = np.max(H2_C_BR)
hosp_max_C_GR2 = np.max(H2_C_GR)
hosp_max_C_NG2 = np.max(H2_C_NG)
hosp_max_C_US2 = np.max(H2_C_US)

fat_rate_BR2 = 1E5 * C2_BR[last_day] / population[0]
fat_rate_GR2 = 1E5 * C2_GR[last_day] / population[1]
fat_rate_NG2 = 1E5 * C2_NG[last_day] / population[2]
fat_rate_US2 = 1E5 * C2_US[last_day] / population[3]

fat_rate_A_BR2 = 1E5 * C2_A_BR[last_day] / population[0]
fat_rate_A_GR2 = 1E5 * C2_A_GR[last_day] / population[1]
fat_rate_A_NG2 = 1E5 * C2_A_NG[last_day] / population[2]
fat_rate_A_US2 = 1E5 * C2_A_US[last_day] / population[3]

fat_rate_B_BR2 = 1E5 * C2_B_BR[last_day] / population[0]
fat_rate_B_GR2 = 1E5 * C2_B_GR[last_day] / population[1]
fat_rate_B_NG2 = 1E5 * C2_B_NG[last_day] / population[2]
fat_rate_B_US2 = 1E5 * C2_B_US[last_day] / population[3]

fat_rate_C_BR2 = 1E5 * C2_C_BR[last_day] / population[0]
fat_rate_C_GR2 = 1E5 * C2_C_GR[last_day] / population[1]
fat_rate_C_NG2 = 1E5 * C2_C_NG[last_day] / population[2]
fat_rate_C_US2 = 1E5 * C2_C_US[last_day] / population[3]

# Scenario 3

an_preval_BR3 = 100 * Ac3_BR[last_day] / population[0]
an_preval_GR3 = 100 * Ac3_GR[last_day] / population[1]
an_preval_NG3 = 100 * Ac3_NG[last_day] / population[2]
an_preval_US3 = 100 * Ac3_US[last_day] / population[3]

an_preval_A_BR3 = 100 * Ac3_A_BR[last_day] / population[0]
an_preval_A_GR3 = 100 * Ac3_A_GR[last_day] / population[1]
an_preval_A_NG3 = 100 * Ac3_A_NG[last_day] / population[2]
an_preval_A_US3 = 100 * Ac3_A_US[last_day] / population[3]

an_preval_B_BR3 = 100 * Ac3_B_BR[last_day] / population[0]
an_preval_B_GR3 = 100 * Ac3_B_GR[last_day] / population[1]
an_preval_B_NG3 = 100 * Ac3_B_NG[last_day] / population[2]
an_preval_B_US3 = 100 * Ac3_B_US[last_day] / population[3]

an_preval_C_BR3 = 100 * Ac3_C_BR[last_day] / population[0]
an_preval_C_GR3 = 100 * Ac3_C_GR[last_day] / population[1]
an_preval_C_NG3 = 100 * Ac3_C_NG[last_day] / population[2]
an_preval_C_US3 = 100 * Ac3_C_US[last_day] / population[3]

hosp_max_BR3 = np.max(H3_BR)
hosp_max_GR3 = np.max(H3_GR)
hosp_max_NG3 = np.max(H3_NG)
hosp_max_US3 = np.max(H3_US)

hosp_max_A_BR3 = np.max(H3_A_BR)
hosp_max_A_GR3 = np.max(H3_A_GR)
hosp_max_A_NG3 = np.max(H3_A_NG)
hosp_max_A_US3 = np.max(H3_A_US)

hosp_max_B_BR3 = np.max(H3_B_BR)
hosp_max_B_GR3 = np.max(H3_B_GR)
hosp_max_B_NG3 = np.max(H3_B_NG)
hosp_max_B_US3 = np.max(H3_B_US)

hosp_max_C_BR3 = np.max(H3_C_BR)
hosp_max_C_GR3 = np.max(H3_C_GR)
hosp_max_C_NG3 = np.max(H3_C_NG)
hosp_max_C_US3 = np.max(H3_C_US)

fat_rate_BR3 = 1E5 * C3_BR[last_day] / population[0]
fat_rate_GR3 = 1E5 * C3_GR[last_day] / population[1]
fat_rate_NG3 = 1E5 * C3_NG[last_day] / population[2]
fat_rate_US3 = 1E5 * C3_US[last_day] / population[3]

fat_rate_A_BR3 = 1E5 * C3_A_BR[last_day] / population[0]
fat_rate_A_GR3 = 1E5 * C3_A_GR[last_day] / population[1]
fat_rate_A_NG3 = 1E5 * C3_A_NG[last_day] / population[2]
fat_rate_A_US3 = 1E5 * C3_A_US[last_day] / population[3]

fat_rate_B_BR3 = 1E5 * C3_B_BR[last_day] / population[0]
fat_rate_B_GR3 = 1E5 * C3_B_GR[last_day] / population[1]
fat_rate_B_NG3 = 1E5 * C3_B_NG[last_day] / population[2]
fat_rate_B_US3 = 1E5 * C3_B_US[last_day] / population[3]

fat_rate_C_BR3 = 1E5 * C3_C_BR[last_day] / population[0]
fat_rate_C_GR3 = 1E5 * C3_C_GR[last_day] / population[1]
fat_rate_C_NG3 = 1E5 * C3_C_NG[last_day] / population[2]
fat_rate_C_US3 = 1E5 * C3_C_US[last_day] / population[3]

# Scenario 4

an_preval_BR4 = 100 * Ac4_BR[last_day] / population[0]
an_preval_GR4 = 100 * Ac4_GR[last_day] / population[1]
an_preval_NG4 = 100 * Ac4_NG[last_day] / population[2]
an_preval_US4 = 100 * Ac4_US[last_day] / population[3]

an_preval_A_BR4 = 100 * Ac4_A_BR[last_day] / population[0]
an_preval_A_GR4 = 100 * Ac4_A_GR[last_day] / population[1]
an_preval_A_NG4 = 100 * Ac4_A_NG[last_day] / population[2]
an_preval_A_US4 = 100 * Ac4_A_US[last_day] / population[3]

an_preval_B_BR4 = 100 * Ac4_B_BR[last_day] / population[0]
an_preval_B_GR4 = 100 * Ac4_B_GR[last_day] / population[1]
an_preval_B_NG4 = 100 * Ac4_B_NG[last_day] / population[2]
an_preval_B_US4 = 100 * Ac4_B_US[last_day] / population[3]

an_preval_C_BR4 = 100 * Ac4_C_BR[last_day] / population[0]
an_preval_C_GR4 = 100 * Ac4_C_GR[last_day] / population[1]
an_preval_C_NG4 = 100 * Ac4_C_NG[last_day] / population[2]
an_preval_C_US4 = 100 * Ac4_C_US[last_day] / population[3]

hosp_max_BR4 = np.max(H4_BR)
hosp_max_GR4 = np.max(H4_GR)
hosp_max_NG4 = np.max(H4_NG)
hosp_max_US4 = np.max(H4_US)

hosp_max_A_BR4 = np.max(H4_A_BR)
hosp_max_A_GR4 = np.max(H4_A_GR)
hosp_max_A_NG4 = np.max(H4_A_NG)
hosp_max_A_US4 = np.max(H4_A_US)

hosp_max_B_BR4 = np.max(H4_B_BR)
hosp_max_B_GR4 = np.max(H4_B_GR)
hosp_max_B_NG4 = np.max(H4_B_NG)
hosp_max_B_US4 = np.max(H4_B_US)

hosp_max_C_BR4 = np.max(H4_C_BR)
hosp_max_C_GR4 = np.max(H4_C_GR)
hosp_max_C_NG4 = np.max(H4_C_NG)
hosp_max_C_US4 = np.max(H4_C_US)

fat_rate_BR4 = 1E5 * C4_BR[last_day] / population[0]
fat_rate_GR4 = 1E5 * C4_GR[last_day] / population[1]
fat_rate_NG4 = 1E5 * C4_NG[last_day] / population[2]
fat_rate_US4 = 1E5 * C4_US[last_day] / population[3]

fat_rate_A_BR4 = 1E5 * C4_A_BR[last_day] / population[0]
fat_rate_A_GR4 = 1E5 * C4_A_GR[last_day] / population[1]
fat_rate_A_NG4 = 1E5 * C4_A_NG[last_day] / population[2]
fat_rate_A_US4 = 1E5 * C4_A_US[last_day] / population[3]

fat_rate_B_BR4 = 1E5 * C4_B_BR[last_day] / population[0]
fat_rate_B_GR4 = 1E5 * C4_B_GR[last_day] / population[1]
fat_rate_B_NG4 = 1E5 * C4_B_NG[last_day] / population[2]
fat_rate_B_US4 = 1E5 * C4_B_US[last_day] / population[3]

fat_rate_C_BR4 = 1E5 * C4_C_BR[last_day] / population[0]
fat_rate_C_GR4 = 1E5 * C4_C_GR[last_day] / population[1]
fat_rate_C_NG4 = 1E5 * C4_C_NG[last_day] / population[2]
fat_rate_C_US4 = 1E5 * C4_C_US[last_day] / population[3]

# Scenario 5

an_preval_BR5 = 100 * Ac5_BR[last_day] / population[0]
an_preval_GR5 = 100 * Ac5_GR[last_day] / population[1]
an_preval_NG5 = 100 * Ac5_NG[last_day] / population[2]
an_preval_US5 = 100 * Ac5_US[last_day] / population[3]

an_preval_A_BR5 = 100 * Ac5_A_BR[last_day] / population[0]
an_preval_A_GR5 = 100 * Ac5_A_GR[last_day] / population[1]
an_preval_A_NG5 = 100 * Ac5_A_NG[last_day] / population[2]
an_preval_A_US5 = 100 * Ac5_A_US[last_day] / population[3]

an_preval_B_BR5 = 100 * Ac5_B_BR[last_day] / population[0]
an_preval_B_GR5 = 100 * Ac5_B_GR[last_day] / population[1]
an_preval_B_NG5 = 100 * Ac5_B_NG[last_day] / population[2]
an_preval_B_US5 = 100 * Ac5_B_US[last_day] / population[3]

an_preval_C_BR5 = 100 * Ac5_C_BR[last_day] / population[0]
an_preval_C_GR5 = 100 * Ac5_C_GR[last_day] / population[1]
an_preval_C_NG5 = 100 * Ac5_C_NG[last_day] / population[2]
an_preval_C_US5 = 100 * Ac5_C_US[last_day] / population[3]

hosp_max_BR5 = np.max(H5_BR)
hosp_max_GR5 = np.max(H5_GR)
hosp_max_NG5 = np.max(H5_NG)
hosp_max_US5 = np.max(H5_US)

hosp_max_A_BR5 = np.max(H5_A_BR)
hosp_max_A_GR5 = np.max(H5_A_GR)
hosp_max_A_NG5 = np.max(H5_A_NG)
hosp_max_A_US5 = np.max(H5_A_US)

hosp_max_B_BR5 = np.max(H5_B_BR)
hosp_max_B_GR5 = np.max(H5_B_GR)
hosp_max_B_NG5 = np.max(H5_B_NG)
hosp_max_B_US5 = np.max(H5_B_US)

hosp_max_C_BR5 = np.max(H5_C_BR)
hosp_max_C_GR5 = np.max(H5_C_GR)
hosp_max_C_NG5 = np.max(H5_C_NG)
hosp_max_C_US5 = np.max(H5_C_US)

fat_rate_BR5 = 1E5 * C5_BR[last_day] / population[0]
fat_rate_GR5 = 1E5 * C5_GR[last_day] / population[1]
fat_rate_NG5 = 1E5 * C5_NG[last_day] / population[2]
fat_rate_US5 = 1E5 * C5_US[last_day] / population[3]

fat_rate_A_BR5 = 1E5 * C5_A_BR[last_day] / population[0]
fat_rate_A_GR5 = 1E5 * C5_A_GR[last_day] / population[1]
fat_rate_A_NG5 = 1E5 * C5_A_NG[last_day] / population[2]
fat_rate_A_US5 = 1E5 * C5_A_US[last_day] / population[3]

fat_rate_B_BR5 = 1E5 * C5_B_BR[last_day] / population[0]
fat_rate_B_GR5 = 1E5 * C5_B_GR[last_day] / population[1]
fat_rate_B_NG5 = 1E5 * C5_B_NG[last_day] / population[2]
fat_rate_B_US5 = 1E5 * C5_B_US[last_day] / population[3]

fat_rate_C_BR5 = 1E5 * C5_C_BR[last_day] / population[0]
fat_rate_C_GR5 = 1E5 * C5_C_GR[last_day] / population[1]
fat_rate_C_NG5 = 1E5 * C5_C_NG[last_day] / population[2]
fat_rate_C_US5 = 1E5 * C5_C_US[last_day] / population[3]

# Scenario 6

an_preval_BR6 = 100 * Ac6_BR[last_day] / population[0]
an_preval_GR6 = 100 * Ac6_GR[last_day] / population[1]
an_preval_NG6 = 100 * Ac6_NG[last_day] / population[2]
an_preval_US6 = 100 * Ac6_US[last_day] / population[3]

an_preval_A_BR6 = 100 * Ac6_A_BR[last_day] / population[0]
an_preval_A_GR6 = 100 * Ac6_A_GR[last_day] / population[1]
an_preval_A_NG6 = 100 * Ac6_A_NG[last_day] / population[2]
an_preval_A_US6 = 100 * Ac6_A_US[last_day] / population[3]

an_preval_B_BR6 = 100 * Ac6_B_BR[last_day] / population[0]
an_preval_B_GR6 = 100 * Ac6_B_GR[last_day] / population[1]
an_preval_B_NG6 = 100 * Ac6_B_NG[last_day] / population[2]
an_preval_B_US6 = 100 * Ac6_B_US[last_day] / population[3]

an_preval_C_BR6 = 100 * Ac6_C_BR[last_day] / population[0]
an_preval_C_GR6 = 100 * Ac6_C_GR[last_day] / population[1]
an_preval_C_NG6 = 100 * Ac6_C_NG[last_day] / population[2]
an_preval_C_US6 = 100 * Ac6_C_US[last_day] / population[3]

hosp_max_BR6 = np.max(H6_BR)
hosp_max_GR6 = np.max(H6_GR)
hosp_max_NG6 = np.max(H6_NG)
hosp_max_US6 = np.max(H6_US)

hosp_max_A_BR6 = np.max(H6_A_BR)
hosp_max_A_GR6 = np.max(H6_A_GR)
hosp_max_A_NG6 = np.max(H6_A_NG)
hosp_max_A_US6 = np.max(H6_A_US)

hosp_max_B_BR6 = np.max(H6_B_BR)
hosp_max_B_GR6 = np.max(H6_B_GR)
hosp_max_B_NG6 = np.max(H6_B_NG)
hosp_max_B_US6 = np.max(H6_B_US)

hosp_max_C_BR6 = np.max(H6_C_BR)
hosp_max_C_GR6 = np.max(H6_C_GR)
hosp_max_C_NG6 = np.max(H6_C_NG)
hosp_max_C_US6 = np.max(H6_C_US)

fat_rate_BR6 = 1E5 * C6_BR[last_day] / population[0]
fat_rate_GR6 = 1E5 * C6_GR[last_day] / population[1]
fat_rate_NG6 = 1E5 * C6_NG[last_day] / population[2]
fat_rate_US6 = 1E5 * C6_US[last_day] / population[3]

fat_rate_A_BR6 = 1E5 * C6_A_BR[last_day] / population[0]
fat_rate_A_GR6 = 1E5 * C6_A_GR[last_day] / population[1]
fat_rate_A_NG6 = 1E5 * C6_A_NG[last_day] / population[2]
fat_rate_A_US6 = 1E5 * C6_A_US[last_day] / population[3]

fat_rate_B_BR6 = 1E5 * C6_B_BR[last_day] / population[0]
fat_rate_B_GR6 = 1E5 * C6_B_GR[last_day] / population[1]
fat_rate_B_NG6 = 1E5 * C6_B_NG[last_day] / population[2]
fat_rate_B_US6 = 1E5 * C6_B_US[last_day] / population[3]

fat_rate_C_BR6 = 1E5 * C6_C_BR[last_day] / population[0]
fat_rate_C_GR6 = 1E5 * C6_C_GR[last_day] / population[1]
fat_rate_C_NG6 = 1E5 * C6_C_NG[last_day] / population[2]
fat_rate_C_US6 = 1E5 * C6_C_US[last_day] / population[3]

# Comandos para plotagem

cenario_tag = 'Scenario 1'

plot_same_scenario(an_preval_BR1, an_preval_GR1, an_preval_NG1, an_preval_US1,
                   an_preval_A_BR1, an_preval_A_GR1, an_preval_A_NG1, an_preval_A_US1,
                   an_preval_B_BR1, an_preval_B_GR1, an_preval_B_NG1, an_preval_B_US1,
                   an_preval_C_BR1, an_preval_C_GR1, an_preval_C_NG1, an_preval_C_US1,
                   hosp_max_BR1, hosp_max_GR1, hosp_max_NG1, hosp_max_US1,
                   hosp_max_A_BR1, hosp_max_A_GR1, hosp_max_A_NG1, hosp_max_A_US1,
                   hosp_max_B_BR1, hosp_max_B_GR1, hosp_max_B_NG1, hosp_max_B_US1,
                   hosp_max_C_BR1, hosp_max_C_GR1, hosp_max_C_NG1, hosp_max_C_US1,
                   fat_rate_BR1, fat_rate_GR1, fat_rate_NG1, fat_rate_US1,
                   fat_rate_A_BR1, fat_rate_A_GR1, fat_rate_A_NG1, fat_rate_A_US1,
                   fat_rate_B_BR1, fat_rate_B_GR1, fat_rate_B_NG1, fat_rate_B_US1,
                   fat_rate_C_BR1, fat_rate_C_GR1, fat_rate_C_NG1, fat_rate_C_US1,
                   cenario_tag)

save_sheet_same_scenario(an_preval_BR1, an_preval_GR1, an_preval_NG1, an_preval_US1,
                         an_preval_A_BR1, an_preval_A_GR1, an_preval_A_NG1, an_preval_A_US1,
                         an_preval_B_BR1, an_preval_B_GR1, an_preval_B_NG1, an_preval_B_US1,
                         an_preval_C_BR1, an_preval_C_GR1, an_preval_C_NG1, an_preval_C_US1,
                         hosp_max_BR1, hosp_max_GR1, hosp_max_NG1, hosp_max_US1,
                         hosp_max_A_BR1, hosp_max_A_GR1, hosp_max_A_NG1, hosp_max_A_US1,
                         hosp_max_B_BR1, hosp_max_B_GR1, hosp_max_B_NG1, hosp_max_B_US1,
                         hosp_max_C_BR1, hosp_max_C_GR1, hosp_max_C_NG1, hosp_max_C_US1,
                         fat_rate_BR1, fat_rate_GR1, fat_rate_NG1, fat_rate_US1,
                         fat_rate_A_BR1, fat_rate_A_GR1, fat_rate_A_NG1, fat_rate_A_US1,
                         fat_rate_B_BR1, fat_rate_B_GR1, fat_rate_B_NG1, fat_rate_B_US1,
                         fat_rate_C_BR1, fat_rate_C_GR1, fat_rate_C_NG1, fat_rate_C_US1,
                         cenario_tag)

cenario_tag = 'Scenario 2'

plot_same_scenario(an_preval_BR2, an_preval_GR2, an_preval_NG2, an_preval_US2,
                   an_preval_A_BR2, an_preval_A_GR2, an_preval_A_NG2, an_preval_A_US2,
                   an_preval_B_BR2, an_preval_B_GR2, an_preval_B_NG2, an_preval_B_US2,
                   an_preval_C_BR2, an_preval_C_GR2, an_preval_C_NG2, an_preval_C_US2,
                   hosp_max_BR2, hosp_max_GR2, hosp_max_NG2, hosp_max_US2,
                   hosp_max_A_BR2, hosp_max_A_GR2, hosp_max_A_NG2, hosp_max_A_US2,
                   hosp_max_B_BR2, hosp_max_B_GR2, hosp_max_B_NG2, hosp_max_B_US2,
                   hosp_max_C_BR2, hosp_max_C_GR2, hosp_max_C_NG2, hosp_max_C_US2,
                   fat_rate_BR2, fat_rate_GR2, fat_rate_NG2, fat_rate_US2,
                   fat_rate_A_BR2, fat_rate_A_GR2, fat_rate_A_NG2, fat_rate_A_US2,
                   fat_rate_B_BR2, fat_rate_B_GR2, fat_rate_B_NG2, fat_rate_B_US2,
                   fat_rate_C_BR2, fat_rate_C_GR2, fat_rate_C_NG2, fat_rate_C_US2,
                   cenario_tag)

save_sheet_same_scenario(an_preval_BR2, an_preval_GR2, an_preval_NG2, an_preval_US2,
                         an_preval_A_BR2, an_preval_A_GR2, an_preval_A_NG2, an_preval_A_US2,
                         an_preval_B_BR2, an_preval_B_GR2, an_preval_B_NG2, an_preval_B_US2,
                         an_preval_C_BR2, an_preval_C_GR2, an_preval_C_NG2, an_preval_C_US2,
                         hosp_max_BR2, hosp_max_GR2, hosp_max_NG2, hosp_max_US2,
                         hosp_max_A_BR2, hosp_max_A_GR2, hosp_max_A_NG2, hosp_max_A_US2,
                         hosp_max_B_BR2, hosp_max_B_GR2, hosp_max_B_NG2, hosp_max_B_US2,
                         hosp_max_C_BR2, hosp_max_C_GR2, hosp_max_C_NG2, hosp_max_C_US2,
                         fat_rate_BR2, fat_rate_GR2, fat_rate_NG2, fat_rate_US2,
                         fat_rate_A_BR2, fat_rate_A_GR2, fat_rate_A_NG2, fat_rate_A_US2,
                         fat_rate_B_BR2, fat_rate_B_GR2, fat_rate_B_NG2, fat_rate_B_US2,
                         fat_rate_C_BR2, fat_rate_C_GR2, fat_rate_C_NG2, fat_rate_C_US2,
                         cenario_tag)

cenario_tag = 'Scenario 3'

plot_same_scenario(an_preval_BR3, an_preval_GR3, an_preval_NG3, an_preval_US3,
                   an_preval_A_BR3, an_preval_A_GR3, an_preval_A_NG3, an_preval_A_US3,
                   an_preval_B_BR3, an_preval_B_GR3, an_preval_B_NG3, an_preval_B_US3,
                   an_preval_C_BR3, an_preval_C_GR3, an_preval_C_NG3, an_preval_C_US3,
                   hosp_max_BR3, hosp_max_GR3, hosp_max_NG3, hosp_max_US3,
                   hosp_max_A_BR3, hosp_max_A_GR3, hosp_max_A_NG3, hosp_max_A_US3,
                   hosp_max_B_BR3, hosp_max_B_GR3, hosp_max_B_NG3, hosp_max_B_US3,
                   hosp_max_C_BR3, hosp_max_C_GR3, hosp_max_C_NG3, hosp_max_C_US3,
                   fat_rate_BR3, fat_rate_GR3, fat_rate_NG3, fat_rate_US3,
                   fat_rate_A_BR3, fat_rate_A_GR3, fat_rate_A_NG3, fat_rate_A_US3,
                   fat_rate_B_BR3, fat_rate_B_GR3, fat_rate_B_NG3, fat_rate_B_US3,
                   fat_rate_C_BR3, fat_rate_C_GR3, fat_rate_C_NG3, fat_rate_C_US3,
                   cenario_tag)

save_sheet_same_scenario(an_preval_BR3, an_preval_GR3, an_preval_NG3, an_preval_US3,
                         an_preval_A_BR3, an_preval_A_GR3, an_preval_A_NG3, an_preval_A_US3,
                         an_preval_B_BR3, an_preval_B_GR3, an_preval_B_NG3, an_preval_B_US3,
                         an_preval_C_BR3, an_preval_C_GR3, an_preval_C_NG3, an_preval_C_US3,
                         hosp_max_BR3, hosp_max_GR3, hosp_max_NG3, hosp_max_US3,
                         hosp_max_A_BR3, hosp_max_A_GR3, hosp_max_A_NG3, hosp_max_A_US3,
                         hosp_max_B_BR3, hosp_max_B_GR3, hosp_max_B_NG3, hosp_max_B_US3,
                         hosp_max_C_BR3, hosp_max_C_GR3, hosp_max_C_NG3, hosp_max_C_US3,
                         fat_rate_BR3, fat_rate_GR3, fat_rate_NG3, fat_rate_US3,
                         fat_rate_A_BR3, fat_rate_A_GR3, fat_rate_A_NG3, fat_rate_A_US3,
                         fat_rate_B_BR3, fat_rate_B_GR3, fat_rate_B_NG3, fat_rate_B_US3,
                         fat_rate_C_BR3, fat_rate_C_GR3, fat_rate_C_NG3, fat_rate_C_US3,
                         cenario_tag)

cenario_tag = 'Scenario 4'

plot_same_scenario(an_preval_BR4, an_preval_GR4, an_preval_NG4, an_preval_US4,
                   an_preval_A_BR4, an_preval_A_GR4, an_preval_A_NG4, an_preval_A_US4,
                   an_preval_B_BR4, an_preval_B_GR4, an_preval_B_NG4, an_preval_B_US4,
                   an_preval_C_BR4, an_preval_C_GR4, an_preval_C_NG4, an_preval_C_US4,
                   hosp_max_BR4, hosp_max_GR4, hosp_max_NG4, hosp_max_US4,
                   hosp_max_A_BR4, hosp_max_A_GR4, hosp_max_A_NG4, hosp_max_A_US4,
                   hosp_max_B_BR4, hosp_max_B_GR4, hosp_max_B_NG4, hosp_max_B_US4,
                   hosp_max_C_BR4, hosp_max_C_GR4, hosp_max_C_NG4, hosp_max_C_US4,
                   fat_rate_BR4, fat_rate_GR4, fat_rate_NG4, fat_rate_US4,
                   fat_rate_A_BR4, fat_rate_A_GR4, fat_rate_A_NG4, fat_rate_A_US4,
                   fat_rate_B_BR4, fat_rate_B_GR4, fat_rate_B_NG4, fat_rate_B_US4,
                   fat_rate_C_BR4, fat_rate_C_GR4, fat_rate_C_NG4, fat_rate_C_US4,
                   cenario_tag)

save_sheet_same_scenario(an_preval_BR4, an_preval_GR4, an_preval_NG4, an_preval_US4,
                         an_preval_A_BR4, an_preval_A_GR4, an_preval_A_NG4, an_preval_A_US4,
                         an_preval_B_BR4, an_preval_B_GR4, an_preval_B_NG4, an_preval_B_US4,
                         an_preval_C_BR4, an_preval_C_GR4, an_preval_C_NG4, an_preval_C_US4,
                         hosp_max_BR4, hosp_max_GR4, hosp_max_NG4, hosp_max_US4,
                         hosp_max_A_BR4, hosp_max_A_GR4, hosp_max_A_NG4, hosp_max_A_US4,
                         hosp_max_B_BR4, hosp_max_B_GR4, hosp_max_B_NG4, hosp_max_B_US4,
                         hosp_max_C_BR4, hosp_max_C_GR4, hosp_max_C_NG4, hosp_max_C_US4,
                         fat_rate_BR4, fat_rate_GR4, fat_rate_NG4, fat_rate_US4,
                         fat_rate_A_BR4, fat_rate_A_GR4, fat_rate_A_NG4, fat_rate_A_US4,
                         fat_rate_B_BR4, fat_rate_B_GR4, fat_rate_B_NG4, fat_rate_B_US4,
                         fat_rate_C_BR4, fat_rate_C_GR4, fat_rate_C_NG4, fat_rate_C_US4,
                         cenario_tag)

cenario_tag = 'Scenario 5'

plot_same_scenario(an_preval_BR5, an_preval_GR5, an_preval_NG5, an_preval_US5,
                   an_preval_A_BR5, an_preval_A_GR5, an_preval_A_NG5, an_preval_A_US5,
                   an_preval_B_BR5, an_preval_B_GR5, an_preval_B_NG5, an_preval_B_US5,
                   an_preval_C_BR5, an_preval_C_GR5, an_preval_C_NG5, an_preval_C_US5,
                   hosp_max_BR5, hosp_max_GR5, hosp_max_NG5, hosp_max_US5,
                   hosp_max_A_BR5, hosp_max_A_GR5, hosp_max_A_NG5, hosp_max_A_US5,
                   hosp_max_B_BR5, hosp_max_B_GR5, hosp_max_B_NG5, hosp_max_B_US5,
                   hosp_max_C_BR5, hosp_max_C_GR5, hosp_max_C_NG5, hosp_max_C_US5,
                   fat_rate_BR5, fat_rate_GR5, fat_rate_NG5, fat_rate_US5,
                   fat_rate_A_BR5, fat_rate_A_GR5, fat_rate_A_NG5, fat_rate_A_US5,
                   fat_rate_B_BR5, fat_rate_B_GR5, fat_rate_B_NG5, fat_rate_B_US5,
                   fat_rate_C_BR5, fat_rate_C_GR5, fat_rate_C_NG5, fat_rate_C_US5,
                   cenario_tag)

save_sheet_same_scenario(an_preval_BR5, an_preval_GR5, an_preval_NG5, an_preval_US5,
                         an_preval_A_BR5, an_preval_A_GR5, an_preval_A_NG5, an_preval_A_US5,
                         an_preval_B_BR5, an_preval_B_GR5, an_preval_B_NG5, an_preval_B_US5,
                         an_preval_C_BR5, an_preval_C_GR5, an_preval_C_NG5, an_preval_C_US5,
                         hosp_max_BR5, hosp_max_GR5, hosp_max_NG5, hosp_max_US5,
                         hosp_max_A_BR5, hosp_max_A_GR5, hosp_max_A_NG5, hosp_max_A_US5,
                         hosp_max_B_BR5, hosp_max_B_GR5, hosp_max_B_NG5, hosp_max_B_US5,
                         hosp_max_C_BR5, hosp_max_C_GR5, hosp_max_C_NG5, hosp_max_C_US5,
                         fat_rate_BR5, fat_rate_GR5, fat_rate_NG5, fat_rate_US5,
                         fat_rate_A_BR5, fat_rate_A_GR5, fat_rate_A_NG5, fat_rate_A_US5,
                         fat_rate_B_BR5, fat_rate_B_GR5, fat_rate_B_NG5, fat_rate_B_US5,
                         fat_rate_C_BR5, fat_rate_C_GR5, fat_rate_C_NG5, fat_rate_C_US5,
                         cenario_tag)

cenario_tag = 'Scenario 6'

plot_same_scenario(an_preval_BR6, an_preval_GR6, an_preval_NG6, an_preval_US6,
                   an_preval_A_BR6, an_preval_A_GR6, an_preval_A_NG6, an_preval_A_US6,
                   an_preval_B_BR6, an_preval_B_GR6, an_preval_B_NG6, an_preval_B_US6,
                   an_preval_C_BR6, an_preval_C_GR6, an_preval_C_NG6, an_preval_C_US6,
                   hosp_max_BR6, hosp_max_GR6, hosp_max_NG6, hosp_max_US6,
                   hosp_max_A_BR6, hosp_max_A_GR6, hosp_max_A_NG6, hosp_max_A_US6,
                   hosp_max_B_BR6, hosp_max_B_GR6, hosp_max_B_NG6, hosp_max_B_US6,
                   hosp_max_C_BR6, hosp_max_C_GR6, hosp_max_C_NG6, hosp_max_C_US6,
                   fat_rate_BR6, fat_rate_GR6, fat_rate_NG6, fat_rate_US6,
                   fat_rate_A_BR6, fat_rate_A_GR6, fat_rate_A_NG6, fat_rate_A_US6,
                   fat_rate_B_BR6, fat_rate_B_GR6, fat_rate_B_NG6, fat_rate_B_US6,
                   fat_rate_C_BR6, fat_rate_C_GR6, fat_rate_C_NG6, fat_rate_C_US6,
                   cenario_tag)

save_sheet_same_scenario(an_preval_BR6, an_preval_GR6, an_preval_NG6, an_preval_US6,
                         an_preval_A_BR6, an_preval_A_GR6, an_preval_A_NG6, an_preval_A_US6,
                         an_preval_B_BR6, an_preval_B_GR6, an_preval_B_NG6, an_preval_B_US6,
                         an_preval_C_BR6, an_preval_C_GR6, an_preval_C_NG6, an_preval_C_US6,
                         hosp_max_BR6, hosp_max_GR6, hosp_max_NG6, hosp_max_US6,
                         hosp_max_A_BR6, hosp_max_A_GR6, hosp_max_A_NG6, hosp_max_A_US6,
                         hosp_max_B_BR6, hosp_max_B_GR6, hosp_max_B_NG6, hosp_max_B_US6,
                         hosp_max_C_BR6, hosp_max_C_GR6, hosp_max_C_NG6, hosp_max_C_US6,
                         fat_rate_BR6, fat_rate_GR6, fat_rate_NG6, fat_rate_US6,
                         fat_rate_A_BR6, fat_rate_A_GR6, fat_rate_A_NG6, fat_rate_A_US6,
                         fat_rate_B_BR6, fat_rate_B_GR6, fat_rate_B_NG6, fat_rate_B_US6,
                         fat_rate_C_BR6, fat_rate_C_GR6, fat_rate_C_NG6, fat_rate_C_US6,
                         cenario_tag)

# Now plot same country, but distinct scenarios

country_tag = 'a'

plot_same_country(an_preval_BR1, an_preval_A_BR1, an_preval_B_BR1, an_preval_C_BR1,
                  an_preval_BR2, an_preval_A_BR2, an_preval_B_BR2, an_preval_C_BR2,
                  an_preval_BR3, an_preval_A_BR3, an_preval_B_BR3, an_preval_C_BR3,
                  an_preval_BR4, an_preval_A_BR4, an_preval_B_BR4, an_preval_C_BR4,
                  an_preval_BR5, an_preval_A_BR5, an_preval_B_BR5, an_preval_C_BR5,
                  an_preval_BR6, an_preval_A_BR6, an_preval_B_BR6, an_preval_C_BR6,
                  hosp_max_BR1, hosp_max_A_BR1, hosp_max_B_BR1, hosp_max_C_BR1,
                  hosp_max_BR2, hosp_max_A_BR2, hosp_max_B_BR2, hosp_max_C_BR2,
                  hosp_max_BR3, hosp_max_A_BR3, hosp_max_B_BR3, hosp_max_C_BR3,
                  hosp_max_BR4, hosp_max_A_BR4, hosp_max_B_BR4, hosp_max_C_BR4,
                  hosp_max_BR5, hosp_max_A_BR5, hosp_max_B_BR5, hosp_max_C_BR5,
                  hosp_max_BR6, hosp_max_A_BR6, hosp_max_B_BR6, hosp_max_C_BR6,
                  fat_rate_BR1, fat_rate_A_BR1, fat_rate_B_BR1, fat_rate_C_BR1,
                  fat_rate_BR2, fat_rate_A_BR2, fat_rate_B_BR2, fat_rate_C_BR2,
                  fat_rate_BR3, fat_rate_A_BR3, fat_rate_B_BR3, fat_rate_C_BR3,
                  fat_rate_BR4, fat_rate_A_BR4, fat_rate_B_BR4, fat_rate_C_BR4,
                  fat_rate_BR5, fat_rate_A_BR5, fat_rate_B_BR5, fat_rate_C_BR5,
                  fat_rate_BR6, fat_rate_A_BR6, fat_rate_B_BR6, fat_rate_C_BR6,
                  country_tag)

save_sheet_same_country(an_preval_BR1, an_preval_A_BR1, an_preval_B_BR1, an_preval_C_BR1,
                        an_preval_BR2, an_preval_A_BR2, an_preval_B_BR2, an_preval_C_BR2,
                        an_preval_BR3, an_preval_A_BR3, an_preval_B_BR3, an_preval_C_BR3,
                        an_preval_BR4, an_preval_A_BR4, an_preval_B_BR4, an_preval_C_BR4,
                        an_preval_BR5, an_preval_A_BR5, an_preval_B_BR5, an_preval_C_BR5,
                        an_preval_BR6, an_preval_A_BR6, an_preval_B_BR6, an_preval_C_BR6,
                        hosp_max_BR1, hosp_max_A_BR1, hosp_max_B_BR1, hosp_max_C_BR1,
                        hosp_max_BR2, hosp_max_A_BR2, hosp_max_B_BR2, hosp_max_C_BR2,
                        hosp_max_BR3, hosp_max_A_BR3, hosp_max_B_BR3, hosp_max_C_BR3,
                        hosp_max_BR4, hosp_max_A_BR4, hosp_max_B_BR4, hosp_max_C_BR4,
                        hosp_max_BR5, hosp_max_A_BR5, hosp_max_B_BR5, hosp_max_C_BR5,
                        hosp_max_BR6, hosp_max_A_BR6, hosp_max_B_BR6, hosp_max_C_BR6,
                        fat_rate_BR1, fat_rate_A_BR1, fat_rate_B_BR1, fat_rate_C_BR1,
                        fat_rate_BR2, fat_rate_A_BR2, fat_rate_B_BR2, fat_rate_C_BR2,
                        fat_rate_BR3, fat_rate_A_BR3, fat_rate_B_BR3, fat_rate_C_BR3,
                        fat_rate_BR4, fat_rate_A_BR4, fat_rate_B_BR4, fat_rate_C_BR4,
                        fat_rate_BR5, fat_rate_A_BR5, fat_rate_B_BR5, fat_rate_C_BR5,
                        fat_rate_BR6, fat_rate_A_BR6, fat_rate_B_BR6, fat_rate_C_BR6,
                        country_tag)

country_tag = 'b'

plot_same_country(an_preval_GR1, an_preval_A_GR1, an_preval_B_GR1, an_preval_C_GR1,
                  an_preval_GR2, an_preval_A_GR2, an_preval_B_GR2, an_preval_C_GR2,
                  an_preval_GR3, an_preval_A_GR3, an_preval_B_GR3, an_preval_C_GR3,
                  an_preval_GR4, an_preval_A_GR4, an_preval_B_GR4, an_preval_C_GR4,
                  an_preval_GR5, an_preval_A_GR5, an_preval_B_GR5, an_preval_C_GR5,
                  an_preval_GR6, an_preval_A_GR6, an_preval_B_GR6, an_preval_C_GR6,
                  hosp_max_GR1, hosp_max_A_GR1, hosp_max_B_GR1, hosp_max_C_GR1,
                  hosp_max_GR2, hosp_max_A_GR2, hosp_max_B_GR2, hosp_max_C_GR2,
                  hosp_max_GR3, hosp_max_A_GR3, hosp_max_B_GR3, hosp_max_C_GR3,
                  hosp_max_GR4, hosp_max_A_GR4, hosp_max_B_GR4, hosp_max_C_GR4,
                  hosp_max_GR5, hosp_max_A_GR5, hosp_max_B_GR5, hosp_max_C_GR5,
                  hosp_max_GR6, hosp_max_A_GR6, hosp_max_B_GR6, hosp_max_C_GR6,
                  fat_rate_GR1, fat_rate_A_GR1, fat_rate_B_GR1, fat_rate_C_GR1,
                  fat_rate_GR2, fat_rate_A_GR2, fat_rate_B_GR2, fat_rate_C_GR2,
                  fat_rate_GR3, fat_rate_A_GR3, fat_rate_B_GR3, fat_rate_C_GR3,
                  fat_rate_GR4, fat_rate_A_GR4, fat_rate_B_GR4, fat_rate_C_GR4,
                  fat_rate_GR5, fat_rate_A_GR5, fat_rate_B_GR5, fat_rate_C_GR5,
                  fat_rate_GR6, fat_rate_A_GR6, fat_rate_B_GR6, fat_rate_C_GR6,
                  country_tag)

save_sheet_same_country(an_preval_GR1, an_preval_A_GR1, an_preval_B_GR1, an_preval_C_GR1,
                        an_preval_GR2, an_preval_A_GR2, an_preval_B_GR2, an_preval_C_GR2,
                        an_preval_GR3, an_preval_A_GR3, an_preval_B_GR3, an_preval_C_GR3,
                        an_preval_GR4, an_preval_A_GR4, an_preval_B_GR4, an_preval_C_GR4,
                        an_preval_GR5, an_preval_A_GR5, an_preval_B_GR5, an_preval_C_GR5,
                        an_preval_GR6, an_preval_A_GR6, an_preval_B_GR6, an_preval_C_GR6,
                        hosp_max_GR1, hosp_max_A_GR1, hosp_max_B_GR1, hosp_max_C_GR1,
                        hosp_max_GR2, hosp_max_A_GR2, hosp_max_B_GR2, hosp_max_C_GR2,
                        hosp_max_GR3, hosp_max_A_GR3, hosp_max_B_GR3, hosp_max_C_GR3,
                        hosp_max_GR4, hosp_max_A_GR4, hosp_max_B_GR4, hosp_max_C_GR4,
                        hosp_max_GR5, hosp_max_A_GR5, hosp_max_B_GR5, hosp_max_C_GR5,
                        hosp_max_GR6, hosp_max_A_GR6, hosp_max_B_GR6, hosp_max_C_GR6,
                        fat_rate_GR1, fat_rate_A_GR1, fat_rate_B_GR1, fat_rate_C_GR1,
                        fat_rate_GR2, fat_rate_A_GR2, fat_rate_B_GR2, fat_rate_C_GR2,
                        fat_rate_GR3, fat_rate_A_GR3, fat_rate_B_GR3, fat_rate_C_GR3,
                        fat_rate_GR4, fat_rate_A_GR4, fat_rate_B_GR4, fat_rate_C_GR4,
                        fat_rate_GR5, fat_rate_A_GR5, fat_rate_B_GR5, fat_rate_C_GR5,
                        fat_rate_GR6, fat_rate_A_GR6, fat_rate_B_GR6, fat_rate_C_GR6,
                        country_tag)

country_tag = 'c'

plot_same_country(an_preval_NG1, an_preval_A_NG1, an_preval_B_NG1, an_preval_C_NG1,
                  an_preval_NG2, an_preval_A_NG2, an_preval_B_NG2, an_preval_C_NG2,
                  an_preval_NG3, an_preval_A_NG3, an_preval_B_NG3, an_preval_C_NG3,
                  an_preval_NG4, an_preval_A_NG4, an_preval_B_NG4, an_preval_C_NG4,
                  an_preval_NG5, an_preval_A_NG5, an_preval_B_NG5, an_preval_C_NG5,
                  an_preval_NG6, an_preval_A_NG6, an_preval_B_NG6, an_preval_C_NG6,
                  hosp_max_NG1, hosp_max_A_NG1, hosp_max_B_NG1, hosp_max_C_NG1,
                  hosp_max_NG2, hosp_max_A_NG2, hosp_max_B_NG2, hosp_max_C_NG2,
                  hosp_max_NG3, hosp_max_A_NG3, hosp_max_B_NG3, hosp_max_C_NG3,
                  hosp_max_NG4, hosp_max_A_NG4, hosp_max_B_NG4, hosp_max_C_NG4,
                  hosp_max_NG5, hosp_max_A_NG5, hosp_max_B_NG5, hosp_max_C_NG5,
                  hosp_max_NG6, hosp_max_A_NG6, hosp_max_B_NG6, hosp_max_C_NG6,
                  fat_rate_NG1, fat_rate_A_NG1, fat_rate_B_NG1, fat_rate_C_NG1,
                  fat_rate_NG2, fat_rate_A_NG2, fat_rate_B_NG2, fat_rate_C_NG2,
                  fat_rate_NG3, fat_rate_A_NG3, fat_rate_B_NG3, fat_rate_C_NG3,
                  fat_rate_NG4, fat_rate_A_NG4, fat_rate_B_NG4, fat_rate_C_NG4,
                  fat_rate_NG5, fat_rate_A_NG5, fat_rate_B_NG5, fat_rate_C_NG5,
                  fat_rate_NG6, fat_rate_A_NG6, fat_rate_B_NG6, fat_rate_C_NG6,
                  country_tag)

save_sheet_same_country(an_preval_NG1, an_preval_A_NG1, an_preval_B_NG1, an_preval_C_NG1,
                        an_preval_NG2, an_preval_A_NG2, an_preval_B_NG2, an_preval_C_NG2,
                        an_preval_NG3, an_preval_A_NG3, an_preval_B_NG3, an_preval_C_NG3,
                        an_preval_NG4, an_preval_A_NG4, an_preval_B_NG4, an_preval_C_NG4,
                        an_preval_NG5, an_preval_A_NG5, an_preval_B_NG5, an_preval_C_NG5,
                        an_preval_NG6, an_preval_A_NG6, an_preval_B_NG6, an_preval_C_NG6,
                        hosp_max_NG1, hosp_max_A_NG1, hosp_max_B_NG1, hosp_max_C_NG1,
                        hosp_max_NG2, hosp_max_A_NG2, hosp_max_B_NG2, hosp_max_C_NG2,
                        hosp_max_NG3, hosp_max_A_NG3, hosp_max_B_NG3, hosp_max_C_NG3,
                        hosp_max_NG4, hosp_max_A_NG4, hosp_max_B_NG4, hosp_max_C_NG4,
                        hosp_max_NG5, hosp_max_A_NG5, hosp_max_B_NG5, hosp_max_C_NG5,
                        hosp_max_NG6, hosp_max_A_NG6, hosp_max_B_NG6, hosp_max_C_NG6,
                        fat_rate_NG1, fat_rate_A_NG1, fat_rate_B_NG1, fat_rate_C_NG1,
                        fat_rate_NG2, fat_rate_A_NG2, fat_rate_B_NG2, fat_rate_C_NG2,
                        fat_rate_NG3, fat_rate_A_NG3, fat_rate_B_NG3, fat_rate_C_NG3,
                        fat_rate_NG4, fat_rate_A_NG4, fat_rate_B_NG4, fat_rate_C_NG4,
                        fat_rate_NG5, fat_rate_A_NG5, fat_rate_B_NG5, fat_rate_C_NG5,
                        fat_rate_NG6, fat_rate_A_NG6, fat_rate_B_NG6, fat_rate_C_NG6,
                        country_tag)

country_tag = 'd'

plot_same_country(an_preval_US1, an_preval_A_US1, an_preval_B_US1, an_preval_C_US1,
                  an_preval_US2, an_preval_A_US2, an_preval_B_US2, an_preval_C_US2,
                  an_preval_US3, an_preval_A_US3, an_preval_B_US3, an_preval_C_US3,
                  an_preval_US4, an_preval_A_US4, an_preval_B_US4, an_preval_C_US4,
                  an_preval_US5, an_preval_A_US5, an_preval_B_US5, an_preval_C_US5,
                  an_preval_US6, an_preval_A_US6, an_preval_B_US6, an_preval_C_US6,
                  hosp_max_US1, hosp_max_A_US1, hosp_max_B_US1, hosp_max_C_US1,
                  hosp_max_US2, hosp_max_A_US2, hosp_max_B_US2, hosp_max_C_US2,
                  hosp_max_US3, hosp_max_A_US3, hosp_max_B_US3, hosp_max_C_US3,
                  hosp_max_US4, hosp_max_A_US4, hosp_max_B_US4, hosp_max_C_US4,
                  hosp_max_US5, hosp_max_A_US5, hosp_max_B_US5, hosp_max_C_US5,
                  hosp_max_US6, hosp_max_A_US6, hosp_max_B_US6, hosp_max_C_US6,
                  fat_rate_US1, fat_rate_A_US1, fat_rate_B_US1, fat_rate_C_US1,
                  fat_rate_US2, fat_rate_A_US2, fat_rate_B_US2, fat_rate_C_US2,
                  fat_rate_US3, fat_rate_A_US3, fat_rate_B_US3, fat_rate_C_US3,
                  fat_rate_US4, fat_rate_A_US4, fat_rate_B_US4, fat_rate_C_US4,
                  fat_rate_US5, fat_rate_A_US5, fat_rate_B_US5, fat_rate_C_US5,
                  fat_rate_US6, fat_rate_A_US6, fat_rate_B_US6, fat_rate_C_US6,
                  country_tag)

save_sheet_same_country(an_preval_US1, an_preval_A_US1, an_preval_B_US1, an_preval_C_US1,
                        an_preval_US2, an_preval_A_US2, an_preval_B_US2, an_preval_C_US2,
                        an_preval_US3, an_preval_A_US3, an_preval_B_US3, an_preval_C_US3,
                        an_preval_US4, an_preval_A_US4, an_preval_B_US4, an_preval_C_US4,
                        an_preval_US5, an_preval_A_US5, an_preval_B_US5, an_preval_C_US5,
                        an_preval_US6, an_preval_A_US6, an_preval_B_US6, an_preval_C_US6,
                        hosp_max_US1, hosp_max_A_US1, hosp_max_B_US1, hosp_max_C_US1,
                        hosp_max_US2, hosp_max_A_US2, hosp_max_B_US2, hosp_max_C_US2,
                        hosp_max_US3, hosp_max_A_US3, hosp_max_B_US3, hosp_max_C_US3,
                        hosp_max_US4, hosp_max_A_US4, hosp_max_B_US4, hosp_max_C_US4,
                        hosp_max_US5, hosp_max_A_US5, hosp_max_B_US5, hosp_max_C_US5,
                        hosp_max_US6, hosp_max_A_US6, hosp_max_B_US6, hosp_max_C_US6,
                        fat_rate_US1, fat_rate_A_US1, fat_rate_B_US1, fat_rate_C_US1,
                        fat_rate_US2, fat_rate_A_US2, fat_rate_B_US2, fat_rate_C_US2,
                        fat_rate_US3, fat_rate_A_US3, fat_rate_B_US3, fat_rate_C_US3,
                        fat_rate_US4, fat_rate_A_US4, fat_rate_B_US4, fat_rate_C_US4,
                        fat_rate_US5, fat_rate_A_US5, fat_rate_B_US5, fat_rate_C_US5,
                        fat_rate_US6, fat_rate_A_US6, fat_rate_B_US6, fat_rate_C_US6,
                        country_tag)
