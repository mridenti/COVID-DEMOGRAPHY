# coding: utf-8
import numpy as np
import csv
import os
import subprocess
import argparse
from scipy import optimize
from scipy import stats
from scipy import linalg
from numpy.linalg import multi_dot

# Optimization script
# Objective function is the chi2 of death and corrected infection notification

##### Process command line options
parser = argparse.ArgumentParser(description='Optimization of notification data based on SEAHIR model. ')
parser.add_argument('-o', '--output_file', help='Output file name', required=True)
parser.add_argument('-oi', '--output_file_i', help='Dummy file to output intermediary results', required=True)
parser.add_argument('-sc', '--scenario', help='Scenario folder', required=True)
parser.add_argument('-td', '--transition_day', nargs=3, help='Transition days, 3 days required', required=True)
parser.add_argument('-ti', '--transition_intervention', nargs=3, help='Codes of intervention levels, 3 values required'
                    , required=True)
parser.add_argument('-ni', '--initial_infected', type=float, help='Initial number of infected', required=True)

args = parser.parse_args()

# Output file - this should be fixed
output_file = args.output_file

# Output file - intermediary analysis
output_file_i = args.output_file_i

# Initial number of infected
ni = args.initial_infected

# Output uncertainty file  - this should be fixed
output_unc_file = 'uncertainty_data.csv'

# Notification data file
notification_file = 'notification.csv'

# Scenario folder
cenario_folder = args.scenario

# Transition days
day_t1 = args.transition_day[0]
day_t2 = args.transition_day[1]
day_t3 = args.transition_day[2]

# Transition interventions
inter_1 = args.transition_intervention[0]
inter_2 = args.transition_intervention[1]
inter_3 = args.transition_intervention[2]

# Number of age groups
age_strata = 16

# Number of days
t_days = 400

# Number of compartments in the output file
compartments = 11

# Read experimental data
# Move to input folder
os.chdir("..")
os.chdir(os.path.join('input', 'cenarios', cenario_folder))
with open(notification_file, "r") as csvfile:
    Spamreader = csv.reader(csvfile, delimiter=',')
    j = 0
    for Row in Spamreader:
        if j == 0:
            YData = Row[1:]
            j = j + 1
        elif j == 1:
            tData = Row[1:]
            j = j + 1
        elif j == 2:
            CData = Row[1:]
            j = j + 1
        elif j == 3:
            tcData = Row[1:]
            j = j + 1
        elif j == 4:
            day_one = Row[1]
            j = j + 1
        elif j == 5:
            month_one = Row[1]
            j = j + 1
        elif j == 6:
            year_one = Row[1]
            j = j + 1
        elif j == 7:
            leitos = Row[1]

YData = np.array([xx for xx in YData if xx], dtype=np.float64)
tData = np.array([xx for xx in tData if xx], dtype=np.int64)
CData = np.array([xx for xx in CData if xx], dtype=np.float64)
tcData = np.array([xx for xx in tcData if xx], dtype=np.int64)

# Return to script directory
os.chdir("..")
os.chdir("..")
os.chdir("..")
os.chdir("scripts")


def save_parameters(r_0, r_0m, r_0p, g_s, g_sm, g_sp, g_e, g_em, g_ep, day, day_m, day_p):
    os.chdir("..")
    os.chdir(os.path.join('input', 'cenarios', cenario_folder))
    with open('optimized_parameters.csv', 'ab') as csv_file:
        spamwriter = csv.writer(csv_file)
        spamwriter.writerow(np.concatenate((['r_0'], [str(r_0)])))
        spamwriter.writerow(np.concatenate((['r_0m'], [str(r_0m)])))
        spamwriter.writerow(np.concatenate((['r_0p'], [str(r_0p)])))
        spamwriter.writerow(np.concatenate((['g_s'], [str(g_s)])))
        spamwriter.writerow(np.concatenate((['g_sm'], [str(g_sm)])))
        spamwriter.writerow(np.concatenate((['g_sp'], [str(g_sp)])))
        spamwriter.writerow(np.concatenate((['g_e'], [str(g_e)])))
        spamwriter.writerow(np.concatenate((['g_es'], [str(g_em)])))
        spamwriter.writerow(np.concatenate((['g_ep'], [str(g_ep)])))
        spamwriter.writerow(np.concatenate((['day_t'], [str(day)])))
        spamwriter.writerow(np.concatenate((['day_tm'], [str(day_m)])))
        spamwriter.writerow(np.concatenate((['day_tp'], [str(day_p)])))
    # Return to script directory
    os.chdir("..")
    os.chdir("..")
    os.chdir("..")
    os.chdir("scripts")
    return 0


def read_output(out_file):
    y = np.zeros([t_days, age_strata], dtype=np.float64)
    h = np.zeros([t_days, age_strata], dtype=np.float64)
    c = np.zeros([t_days, age_strata], dtype=np.float64)
    ri = np.zeros([t_days, age_strata], dtype=np.float64)
    ac = np.zeros(t_days, dtype=np.float64)
    with open(out_file, "r") as csvfile:
        spamreader = csv.reader(csvfile, delimiter=',')
        next(spamreader, None)
        j = 0
        for row in spamreader:
            for i in range(age_strata):
                y[j, i] = row[compartments * (i + 1) + 2]
                c[j, i] = row[compartments * (i + 1) + 6]
                h[j, i] = row[compartments * (i + 1) + 7]
            for ii in range(age_strata):
                ri[j, ii] = row[compartments * (age_strata + 1) + ii + 1]
            j = j + 1
        c_s = np.sum(c, axis=1)
        ri_s = np.sum(ri, axis=1)
        y_s = np.sum(y, axis=1)
        h_s = np.sum(h, axis=1)
        ac = c_s + ri_s + y_s + h_s
        return ac, c_s


def build_cov_matrix(data):
    length = len(data)
    s = np.eye(length)
    for i in range(0, length):
        s[i, i] = data[i]
    return s


def f_cost(vec):
    l_vec = len(vec)
    new_vec = np.zeros(l_vec, dtype=np.float64)
    for i in range(0, l_vec):
        if vec[i] >= 0:
            new_vec[i] = 1
        else:
            new_vec[i] = 1 + np.square(vec[i])
    return new_vec


def f(x):
    r_0 = x[0]
    g_s = x[1]
    g_e = x[2]
    day_t = x[3]
    print(r_0, g_s, g_e, day_t)
    subprocess.call(['python', 'cenario_generator.py', '-i', cenario_folder, '-d', '0', str(day_t), day_t2, day_t3,
                     '-m', '3', '-I0', str(ni * g_s), '-R0', str(r_0), '-Rp', str(r_0), '-p', '1', str(g_e),
                     str(g_e), str(g_e), '-itv', '0', inter_1, inter_2, inter_3], stdout=open(os.devnull, 'wb'))
    os.chdir("..")
    subprocess.call(['bin/csv_to_input', cenario_folder], stdout=open(os.devnull, 'wb'))
    subprocess.call(['bin/spatial_covid0d_estrat.exe', 'input/generated-input.txt', '/'.join(['output', output_file]),
                     '3'], stdout=open(os.devnull, 'wb'))
    os.chdir("output")
    f_ac, f_c = read_output(output_file)
    os.chdir("..")
    os.chdir("scripts")

    ac_vec = np.take(f_ac, tData - 1) / g_s - YData
    s_ac = build_cov_matrix(YData)
    term1 = multi_dot([ac_vec, linalg.inv(s_ac), ac_vec])

    c_vec = np.take(f_c, tcData - 1) - CData
    c_vec = np.multiply(c_vec, f_cost(c_vec))
    s_c = build_cov_matrix(CData)
    term2 = multi_dot([c_vec, linalg.inv(s_c), c_vec])

    f_chi2 = term1 + term2

    return f_chi2


def f_s(x, params):
    g_s = x[0]
    r_0 = params[0]
    g_e = params[1]
    day_t = params[2]
    dead_fac = params[3]
    subprocess.call(['python', 'cenario_generator.py', '-i', cenario_folder, '-d', '0', str(day_t), day_t2, day_t3,
                     '-m', '3', '-I0', str(ni * g_s), '-R0', str(r_0), '-Rp', str(r_0), '-p', '1', str(g_e),
                     str(g_e), str(g_e), '-itv', '0', inter_1, inter_2, inter_3, '-f', str(dead_fac)],
                    stdout=open(os.devnull, 'wb'))
    os.chdir("..")
    subprocess.call(['bin/csv_to_input', cenario_folder], stdout=open(os.devnull, 'wb'))
    subprocess.call(['bin/spatial_covid0d_estrat.exe', 'input/generated-input.txt',
                     '/'.join(['output', output_unc_file]), '3'], stdout=open(os.devnull, 'wb'))
    os.chdir("output")
    f_ac, f_c = read_output(output_unc_file)
    os.chdir("..")
    os.chdir("scripts")

    ac_vec = np.take(f_ac, tData - 1) / g_s - YData
    s_ac = build_cov_matrix(YData)
    term1 = multi_dot([ac_vec, linalg.inv(s_ac), ac_vec])

    c_vec = np.take(f_c, tcData - 1) - CData
    c_vec = np.multiply(c_vec, f_cost(c_vec))
    s_c = build_cov_matrix(CData)
    term2 = multi_dot([c_vec, linalg.inv(s_c), c_vec])

    f_chi2 = term1 + term2

    return f_chi2


def f_g(x, params):
    r_0 = x[0]
    g_s = x[1]
    g_e = x[2]
    day_t = x[3]
    dead_fac = params
    subprocess.call(['python', 'cenario_generator.py', '-i', cenario_folder, '-d', '0', str(day_t), day_t2, day_t3,
                     '-m', '3', '-I0', str(ni * g_s), '-R0', str(r_0), '-Rp', str(r_0), '-p', '1', str(g_e),
                     str(g_e), str(g_e), '-itv', '0', inter_1, inter_2, inter_3, '-f', str(dead_fac)],
                    stdout=open(os.devnull, 'wb'))
    os.chdir("..")
    subprocess.call(['bin/csv_to_input', cenario_folder], stdout=open(os.devnull, 'wb'))
    subprocess.call(['bin/spatial_covid0d_estrat.exe', 'input/generated-input.txt',
                     '/'.join(['output', output_unc_file]), '3'], stdout=open(os.devnull, 'wb'))
    os.chdir("output")
    f_ac, f_c = read_output(output_unc_file)
    os.chdir("..")
    os.chdir("scripts")

    ac_vec = np.take(f_ac, tData - 1) / g_s - YData
    s_ac = build_cov_matrix(YData)
    term1 = multi_dot([ac_vec, linalg.inv(s_ac), ac_vec])

    c_vec = np.take(f_c, tcData - 1) - CData
    c_vec = np.multiply(c_vec, f_cost(c_vec))
    s_c = build_cov_matrix(CData)
    term2 = multi_dot([c_vec, linalg.inv(s_c), c_vec])

    f_chi2 = term1 + term2

    return f_chi2


def pchi2_limit(chi2_0, ngl, par, d_bound, fat_fac):
    r_0 = par[0]
    g_s = par[1]
    g_e = par[2]
    day_t = par[3]
    day_bound_plus = d_bound[1]
    day_bound_minus = d_bound[0]
    # plus limit of R0
    delta_par = 0.001 * par[0]
    prob0 = stats.chi2.cdf(ngl, ngl)
    prob = prob0
    while 0.95 > prob > 0.05:
        r_0 = r_0 + delta_par
        chi2_upd = ngl * f_s([g_s], [r_0, g_e, day_t, fat_fac]) / chi2_0
        prob = stats.chi2.cdf(chi2_upd, ngl)
        err_prob = (prob - 0.95) / 0.95
        if prob > 0.95 and err_prob > 0.01:
            delta_par = 0.5 * delta_par
            prob = 0.95
    r0_plus = r_0
    r_0 = par[0]
    # minus limit of R0
    delta_par = 0.001 * par[0]
    prob = prob0
    while 0.95 > prob > 0.05:
        r_0 = r_0 - delta_par
        chi2_upd = ngl * f_s([g_s], [r_0, g_e, day_t, fat_fac]) / chi2_0
        prob = stats.chi2.cdf(chi2_upd, ngl)
        err_prob = (prob - 0.95) / 0.95
        if prob > 0.95 and err_prob > 0.01:
            delta_par = 0.5 * delta_par
            prob = 0.95
    r0_minus = r_0
    r_0 = par[0]
    # plus limit of g_e
    delta_par = 0.001 * par[2]
    prob = prob0
    while 0.95 > prob > 0.05:
        g_e = g_e + delta_par
        chi2_upd = ngl * f_s([g_s], [r_0, g_e, day_t, fat_fac]) / chi2_0
        prob = stats.chi2.cdf(chi2_upd, ngl)
        err_prob = (prob - 0.95) / 0.95
        if prob > 0.95 and err_prob > 0.01:
            delta_par = 0.5 * delta_par
            prob = 0.95
    g_e_plus = g_e
    g_e = par[2]
    # minus limit of g_e
    delta_par = 0.001 * par[2]
    prob = prob0
    while 0.95 > prob > 0.05:
        g_e = g_e - delta_par
        chi2_upd = ngl * f_s([g_s], [r_0, g_e, day_t, fat_fac]) / chi2_0
        prob = stats.chi2.cdf(chi2_upd, ngl)
        err_prob = (prob - 0.95) / 0.95
        if prob > 0.95 and err_prob > 0.01:
            delta_par = 0.5 * delta_par
            prob = 0.95
    g_e_minus = g_e
    g_e = par[2]
    # plus limit of day_t
    delta_par = 0.01 * par[3]
    prob = prob0
    while 0.95 > prob > 0.05 and day_t < day_bound_plus:
        day_t = day_t + delta_par
        chi2_upd = ngl * f_s([g_s], [r_0, g_e, day_t, fat_fac]) / chi2_0
        prob = stats.chi2.cdf(chi2_upd, ngl)
        err_prob = (prob - 0.95) / 0.95
        if prob > 0.95 and err_prob > 0.01:
            delta_par = 0.5 * delta_par
            prob = 0.95
    day_t_plus = day_t
    day_t = par[3]
    # minus limit of day_t
    delta_par = 0.01 * par[3]
    prob = prob0
    while 0.95 > prob > 0.05 and day_t > day_bound_minus:
        day_t = day_t - delta_par
        chi2_upd = ngl * f_s([g_s], [r_0, g_e, day_t, fat_fac]) / chi2_0
        prob = stats.chi2.cdf(chi2_upd, ngl)
        err_prob = (prob - 0.95) / 0.95
        if prob > 0.95 and err_prob > 0.01:
            delta_par = 0.5 * delta_par
            prob = 0.95
    day_t_minus = day_t

    return r0_plus, r0_minus, g_e_plus, g_e_minus, day_t_plus, day_t_minus


r_bound = [1.0, 20.0]
g_bound = [1.0, 60.0]
g_e_bound = [0.1, 1.0]
day_bound = [3.0, 40]
ret = optimize.differential_evolution(f, [r_bound, g_bound, g_e_bound, day_bound], disp=True,
                                      maxiter=20, popsize=15, seed=1234)

# Uncertainty estimate:
# First: estimate of the R0 and attenuation factor uncertainties
print("Uncertainty estimate of R0 and attenuation factor:")

df = len(YData) + len(CData) - 4
chi2_init = ret.fun
fat_fac = 2.0
r0_Plus, r0_Minus, g_e_Plus, g_e_Minus, day_t_Plus, day_t_Minus = pchi2_limit(chi2_init, df,
                                                                              [ret.x[0], ret.x[1], ret.x[2], ret.x[3]],
                                                                              day_bound, fat_fac)

# Second: estimate of the under notification factor uncertainty
print("Uncertainty estimate of the underreporting index:")
minus_fac = 1.0
plus_fac = 3.0
g_unc_bound = [1.1*ret.x[1], 3 * ret.x[1]]
g_unc_bound_minus = [1.0, 0.9*ret.x[1]]
ret_minus = optimize.differential_evolution(f_g, [r_bound, g_unc_bound, g_e_bound, day_bound], args=(minus_fac,),
                                            disp=True, maxiter=20, seed=1234)
ret_plus = optimize.differential_evolution(f_g, [r_bound, g_unc_bound_minus, g_e_bound, day_bound], args=(plus_fac,),
                                           disp=True, maxiter=20, seed=1234)

# r0_Plus = np.max([ret_plus.x[0], ret_minus.x[0], r0_Plus])
# r0_Minus = np.min([ret_plus.x[0], ret_minus.x[0], r0_Minus])
# g_e_Plus = np.max([ret_plus.x[2], ret_minus.x[2], g_e_Plus])
# g_e_Minus = np.min([ret_plus.x[2], ret_minus.x[2], g_e_Minus])
# day_t_Plus = np.max([ret_plus.x[3], ret_minus.x[3], day_t_Plus])
# day_t_Minus = np.min([ret_plus.x[3], ret_minus.x[3], day_t_Minus])

day_t1 = str(int(np.ceil(ret.x[3])))
g_e0 = str(1.0)
g_e1 = str(ret.x[2])
g_e2 = str(ret.x[2])
g_e3 = str(ret.x[2])

subprocess.call(['python', 'cenario_generator.py', '-i', cenario_folder, '-d', '0', day_t1, day_t2, day_t3, '-m', '3',
                 '-I0', str(ni * ret.x[1]), '-R0', str(ret.x[0]), '-Rp', str(ret.x[0]), '-p', '1',
                 str(ret.x[2]), str(ret.x[2]), str(ret.x[2]), '-itv', '0', inter_1, inter_2, inter_3, '-s'])

print("global minimum: x = [%.4f, %.4f, %.4f, %.4f], f(x0) = %.4f" % (ret.x[0], ret.x[1], ret.x[2], ret.x[3], ret.fun))
print ("Under notification index uncertainty: [%.4f, %.4f]" % (ret_minus.x[1], ret_plus.x[1]))
print ("R0 uncertainty: [%.4f, %.4f]" % (r0_Plus, r0_Minus))
print ("Attenuation factor uncertainty: [%.4f, %.4f]" % (g_e_Plus, g_e_Minus))
print ("Transition day uncertainty: [%.1f, %.1f]" % (day_t_Plus, day_t_Minus))

save_parameters(ret.x[0], r0_Minus, r0_Plus, ret.x[1], ret_plus.x[1], ret_minus.x[1], ret.x[2], g_e_Minus, g_e_Plus,
                ret.x[3], day_t_Minus, day_t_Plus)

os.chdir("..")
subprocess.call(['bin/csv_to_input', cenario_folder], stdout=open(os.devnull, 'wb'))
subprocess.call(['bin/spatial_covid0d_estrat.exe', 'input/generated-input.txt', '/'.join(['output', output_file]),
                 '3'], stdout=open(os.devnull, 'wb'))
os.chdir("scripts")

subprocess.call(['python', 'plot_output_SEAHIR.py', '-d', '0', day_t1, day_t2, day_t3, '-s', str(ret.x[1]),
                 '-o', output_file, '-sc', cenario_folder, '-i', day_one, month_one, year_one])

subprocess.call(['python', 'scenario_explore_A.py', '-d', day_t1, '-b', leitos,
                 '-o', output_file_i, '-sc', cenario_folder, '-i', day_one, month_one, year_one, '-p', g_e0, g_e1,
                 g_e2, g_e3, '-ni', str(ni)])

subprocess.call(['python', 'scenario_explore_B.py', '-d', day_t1, '-b', leitos,
                 '-o', output_file_i, '-sc', cenario_folder, '-i', day_one, month_one, year_one, '-p', g_e0, g_e1,
                 g_e2, g_e3, '-ni', str(ni)])

subprocess.call(['python', 'scenario_explore_C.py', '-d', day_t1, '-b', leitos,
                 '-o', output_file_i, '-sc', cenario_folder, '-i', day_one, month_one, year_one, '-p', g_e0, g_e1,
                 g_e2, g_e3, '-ni', str(ni)])

subprocess.call(['python', 'range_plot.py', '-d', day_t1, '-b', leitos,
                 '-o', output_file_i, '-sc', cenario_folder, '-i', day_one, month_one, year_one, '-p', g_e0, g_e1,
                 g_e2, g_e3, '-ni', str(ni)])
