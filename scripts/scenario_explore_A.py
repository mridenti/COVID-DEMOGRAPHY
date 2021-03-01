# coding: utf-8
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import numpy as np
import csv
import os
import subprocess
import datetime
import argparse

##### Process command line options
##### Variable parameters, for error estimation within reasonable bounds
parser = argparse.ArgumentParser(description=u'This script plots the projection under set A of interventions.')
parser.add_argument('-d', '--day', type=int, help='Day of first intervention',
                    required=True)
parser.add_argument('-b', '--beds', type=int, help='Available ICUs beds',
                    required=False)
parser.add_argument('-o', '--output_file', help='Dummy file to output intermediary results', required=True)
parser.add_argument('-sc', '--scenario', help='Scenario folder', required=True)
parser.add_argument('-i', '--init_date', type=int, nargs=3, help='Date of the first day of notification [dat], [month],'
                                                                 '[year]', required=True)
parser.add_argument('-p', '--prob', type=float, nargs=4, help='Transmission probability attenuation factor',
                    required=True)
parser.add_argument('-ni', '--initial_infected', type=float, help='Initial number of infected', required=True)

args = parser.parse_args()

# Output file - this should be fixed
output_file = args.output_file

# Cenario folder
cenario_folder = args.scenario

# Number of age groups
age_strata = 16

# Number of days
t_days = 400

# Number of compartments in the output file
compartments = 11

# Initial number of infected
ni = args.initial_infected

# Time specifications for plot
months = mdates.MonthLocator()  # every month
weeks = mdates.WeekdayLocator()  # every week
month_fmt = mdates.DateFormatter('%b')

# Days of transition
day_init = 0
day_next_1 = int(args.day)
day = int(args.init_date[0])
month = int(args.init_date[1])
year = int(args.init_date[2])

temp_date = datetime.date(2020, 05, 10) - datetime.date(year, month, day)
day_next_2 = temp_date.days
temp_date = datetime.date(2020, 07, 31) - datetime.date(year, month, day)
day_next_3 = temp_date.days
temp_date = datetime.date(2020, 06, 30) - datetime.date(year, month, day)
day_mid = temp_date.days

g_e0 = float(args.prob[0])
g_e1 = float(args.prob[1])
g_e2 = float(args.prob[2])
g_e3 = float(args.prob[3])

# ICU beds
leitos = int(args.beds)

# Read optimal parameters and uncertainties
os.chdir("..")
os.chdir(os.path.join('input', 'cenarios', cenario_folder))
with open('optimized_parameters.csv', 'r') as csvfile:
    Spamreader = csv.reader(csvfile, delimiter=',')
    j = 0
    for Row in Spamreader:
        if j == 0:
            r0_post = float(Row[1])
            j = j + 1
        elif j == 1:
            r_0 = float(Row[1])
            j = j + 1
        elif j == 2:
            r_0m = float(Row[1])
            j = j + 1
        elif j == 3:
            r_0p = float(Row[1])
            j = j + 1
        elif j == 4:
            g_s = float(Row[1])
            j = j + 1
        elif j == 5:
            g_sm = float(Row[1])
            j = j + 1
        elif j == 6:
            g_sp = float(Row[1])
            j = j + 1
        elif j == 7:
            g_e = float(Row[1])
            j = j + 1
        elif j == 8:
            g_em = float(Row[1])
            j = j + 1
        elif j == 9:
            g_ep = float(Row[1])
            j = j + 1

if r0_post < 2.0:
    r0_post = 2.0

# Return to script directory
os.chdir("..")
os.chdir("..")
os.chdir("..")
os.chdir("scripts")


def save_fig(suffix_fig):
    os.chdir("..")
    my_path = os.getcwd()
    dirs = os.path.join(my_path, 'figures', cenario_folder)
    try:
        os.makedirs(dirs)
    except OSError:
        pass
    plt.savefig(os.path.join(dirs, cenario_folder[:] + suffix_fig), format='svg',
                dpi=300, bbox_inches='tight')
    os.chdir("scripts")
    return 0


def read_output(out_file):
    s = np.zeros([t_days, age_strata], dtype=np.float64)
    e = np.zeros([t_days, age_strata], dtype=np.float64)
    y = np.zeros([t_days, age_strata], dtype=np.float64)
    r = np.zeros([t_days, age_strata], dtype=np.float64)
    n = np.zeros([t_days, age_strata], dtype=np.float64)
    a = np.zeros([t_days, age_strata], dtype=np.float64)
    c = np.zeros([t_days, age_strata], dtype=np.float64)
    h = np.zeros([t_days, age_strata], dtype=np.float64)
    l = np.zeros([t_days, age_strata], dtype=np.float64)
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
                l[j, i] = row[compartments * (i + 1) + 8]
            for ii in range(age_strata):
                ri[j, ii] = row[compartments * (age_strata + 1) + ii + 1]
            j = j + 1
        return s, e, y, r, n, a, c, h, l, ri


# First scenario
# Lock, never relax until the 31 of July
subprocess.call(['python', 'cenario_generator.py', '-i', cenario_folder, '-d', str(day_init), str(day_next_1),
                 str(day_next_2), str(day_next_3), '-m', '3',
                 '-I0', str(ni*g_s), '-R0', str(r_0), '-Rp', str(r_0), '-p', str(g_e0), str(g_e1), str(g_e1), str(g_e0),
                 '-itv', '0', '10', '10', '0'])
os.chdir("..")
subprocess.call(['bin/csv_to_input', cenario_folder], stdout=open(os.devnull, 'wb'))
subprocess.call(['bin/spatial_covid0d_estrat.exe', 'input/generated-input.txt', '/'.join(['output', output_file]),
                 '3'], stdout=open(os.devnull, 'wb'))
os.chdir("output")
s1, e1, y1, r1, n1, a1, c1, h1, l1, ri1 = read_output(output_file)
os.chdir("..")
os.chdir("scripts")

# Lock, then relax to level 1 on 10 of May, then open on 31st of July
subprocess.call(['python', 'cenario_generator.py', '-i', cenario_folder, '-d', str(day_init), str(day_next_1),
                 str(day_next_2), str(day_next_3), '-m', '3',
                 '-I0', str(ni*g_s), '-R0', str(r_0), '-Rp', str(r_0), '-p', str(g_e0), str(g_e1), str(g_e1), str(g_e0),
                 '-itv', '0', '10', '9', '0'])
os.chdir("..")
subprocess.call(['bin/csv_to_input', cenario_folder], stdout=open(os.devnull, 'wb'))
subprocess.call(['bin/spatial_covid0d_estrat.exe', 'input/generated-input.txt', '/'.join(['output', output_file]),
                 '3'], stdout=open(os.devnull, 'wb'))
os.chdir("output")
s2, e2, y2, r2, n2, a2, c2, h2, l2, ri2 = read_output(output_file)
os.chdir("..")
os.chdir("scripts")

# Lock, then relax to level 2 on 10 of May, then open on 31st of July
subprocess.call(['python', 'cenario_generator.py', '-i', cenario_folder, '-d', str(day_init), str(day_next_1),
                 str(day_next_2), str(day_next_3), '-m', '3',
                 '-I0', str(ni*g_s), '-R0', str(r_0), '-Rp', str(r_0), '-p', str(g_e0), str(g_e1), str(g_e1), str(g_e0),
                 '-itv', '0', '10', '8', '0'])
os.chdir("..")
subprocess.call(['bin/csv_to_input', cenario_folder], stdout=open(os.devnull, 'wb'))
subprocess.call(['bin/spatial_covid0d_estrat.exe', 'input/generated-input.txt', '/'.join(['output', output_file]),
                 '3'], stdout=open(os.devnull, 'wb'))
os.chdir("output")
s3, e3, y3, r3, n3, a3, c3, h3, l3, ri3 = read_output(output_file)
os.chdir("..")
os.chdir("scripts")

# Lock, then relax to level 1 on 10 of May, relax to level 2 on 30 of June, then open on 31 of July
subprocess.call(['python', 'cenario_generator.py', '-i', cenario_folder, '-d', str(day_init), str(day_next_1),
                 str(day_next_2), str(day_mid), '-m', '3',
                 '-I0', str(ni*g_s), '-R0', str(r_0), '-Rp', str(r_0), '-p', str(g_e0), str(g_e1), str(g_e1), str(g_e1),
                 '-itv', '0', '10', '9', '8', '-ex', '0', str(day_next_3)])
os.chdir("..")
subprocess.call(['bin/csv_to_input', cenario_folder], stdout=open(os.devnull, 'wb'))
subprocess.call(['bin/spatial_covid0d_estrat.exe', 'input/generated-input.txt', '/'.join(['output', output_file]),
                 '3'], stdout=open(os.devnull, 'wb'))
os.chdir("output")
s4, e4, y4, r4, n4, a4, c4, h4, l4, ri4 = read_output(output_file)
os.chdir("..")
os.chdir("scripts")

# Lock, then simply open on 10 of May
subprocess.call(['python', 'cenario_generator.py', '-i', cenario_folder, '-d', str(day_init), str(day_next_1),
                 str(day_next_2), str(day_next_3), '-m', '3',
                 '-I0', str(ni*g_s), '-R0', str(r_0), '-Rp', str(r_0), '-p', str(g_e0), str(g_e1), str(g_e0), str(g_e0),
                 '-itv', '0', '10', '0', '0'])
os.chdir("..")
subprocess.call(['bin/csv_to_input', cenario_folder], stdout=open(os.devnull, 'wb'))
subprocess.call(['bin/spatial_covid0d_estrat.exe', 'input/generated-input.txt', '/'.join(['output', output_file]),
                 '3'], stdout=open(os.devnull, 'wb'))
os.chdir("output")
s5, e5, y5, r5, n5, a5, c5, h5, l5, ri5 = read_output(output_file)
os.chdir("..")
os.chdir("scripts")

# Lock until the 31 of July, then relax to level 3 until end of the year
subprocess.call(['python', 'cenario_generator.py', '-i', cenario_folder, '-d', str(day_init), str(day_next_1),
                 str(day_next_2), str(day_next_3), '-m', '3',
                 '-I0', str(ni*g_s), '-R0', str(r_0), '-Rp', str(r_0), '-p', str(g_e0), str(g_e1), str(g_e1), str(g_e1),
                 '-itv', '0', '10', '10', '7'])
os.chdir("..")
subprocess.call(['bin/csv_to_input', cenario_folder], stdout=open(os.devnull, 'wb'))
subprocess.call(['bin/spatial_covid0d_estrat.exe', 'input/generated-input.txt', '/'.join(['output', output_file]),
                 '3'], stdout=open(os.devnull, 'wb'))
os.chdir("output")
s6, e6, y6, r6, n6, a6, c6, h6, l6, ri6 = read_output(output_file)
os.chdir("..")
os.chdir("scripts")

# Do nothing, wild behaviour scenario
subprocess.call(['python', 'cenario_generator.py', '-i', cenario_folder, '-d', str(day_init), str(day_next_1),
                 str(day_next_2), str(day_next_3), '-m', '3',
                 '-I0', str(ni*g_s), '-R0', str(r_0), '-Rp', str(r_0), '-p', str(g_e0), str(g_e0), str(g_e0), str(g_e0),
                 '-itv', '0', '1', '1', '1'])
os.chdir("..")
subprocess.call(['bin/csv_to_input', cenario_folder], stdout=open(os.devnull, 'wb'))
subprocess.call(['bin/spatial_covid0d_estrat.exe', 'input/generated-input.txt', '/'.join(['output', output_file]),
                 '3'], stdout=open(os.devnull, 'wb'))
os.chdir("output")
s0, e0, y0, r0, n0, a0, c0, h0, l0, ri0 = read_output(output_file)
os.chdir("..")
os.chdir("scripts")

# Somam-se todas as faixas etárias

S1 = np.sum(s1, axis=1)
E1 = np.sum(e1, axis=1)
I1 = np.sum(y1, axis=1)
R1 = np.sum(r1, axis=1)
N1 = np.sum(n1, axis=1)
A1 = np.sum(a1, axis=1)
C1 = np.sum(c1, axis=1)
H1 = np.sum(h1, axis=1)
L1 = np.sum(l1, axis=1)
RI1 = np.sum(ri1, axis=1)

S2 = np.sum(s2, axis=1)
E2 = np.sum(e2, axis=1)
I2 = np.sum(y2, axis=1)
R2 = np.sum(r2, axis=1)
N2 = np.sum(n2, axis=1)
A2 = np.sum(a2, axis=1)
C2 = np.sum(c2, axis=1)
H2 = np.sum(h2, axis=1)
L2 = np.sum(l2, axis=1)
RI2 = np.sum(ri2, axis=1)

S3 = np.sum(s3, axis=1)
E3 = np.sum(e3, axis=1)
I3 = np.sum(y3, axis=1)
R3 = np.sum(r3, axis=1)
N3 = np.sum(n3, axis=1)
A3 = np.sum(a3, axis=1)
C3 = np.sum(c3, axis=1)
H3 = np.sum(h3, axis=1)
L3 = np.sum(l3, axis=1)
RI3 = np.sum(ri3, axis=1)

S4 = np.sum(s4, axis=1)
E4 = np.sum(e4, axis=1)
I4 = np.sum(y4, axis=1)
R4 = np.sum(r4, axis=1)
N4 = np.sum(n4, axis=1)
A4 = np.sum(a4, axis=1)
C4 = np.sum(c4, axis=1)
H4 = np.sum(h4, axis=1)
L4 = np.sum(l4, axis=1)
RI4 = np.sum(ri4, axis=1)

S5 = np.sum(s5, axis=1)
E5 = np.sum(e5, axis=1)
I5 = np.sum(y5, axis=1)
R5 = np.sum(r5, axis=1)
N5 = np.sum(n5, axis=1)
A5 = np.sum(a5, axis=1)
C5 = np.sum(c5, axis=1)
H5 = np.sum(h5, axis=1)
L5 = np.sum(l5, axis=1)
RI5 = np.sum(ri5, axis=1)

S6 = np.sum(s6, axis=1)
E6 = np.sum(e6, axis=1)
I6 = np.sum(y6, axis=1)
R6 = np.sum(r6, axis=1)
N6 = np.sum(n6, axis=1)
A6 = np.sum(a6, axis=1)
C6 = np.sum(c6, axis=1)
H6 = np.sum(h6, axis=1)
L6 = np.sum(l6, axis=1)
RI6 = np.sum(ri6, axis=1)

S0 = np.sum(s0, axis=1)
E0 = np.sum(e0, axis=1)
I0 = np.sum(y0, axis=1)
R0 = np.sum(r0, axis=1)
N0 = np.sum(n0, axis=1)
A0 = np.sum(a0, axis=1)
C0 = np.sum(c0, axis=1)
H0 = np.sum(h0, axis=1)
L0 = np.sum(l0, axis=1)
RI0 = np.sum(ri0, axis=1)

# Plot
t_array = np.linspace(0, t_days - 1, t_days, dtype=int)
tempDate = datetime.datetime(year, month, day)
t_dates = (tempDate + t_array * datetime.timedelta(days=1))
day_1 = (tempDate + day_next_1 * datetime.timedelta(days=1))
day_2 = (tempDate + day_next_2 * datetime.timedelta(days=1))
day_3 = (tempDate + day_next_3 * datetime.timedelta(days=1))

# Infected
_, ax = plt.subplots(figsize=(9, 6))
plt.axvline(day_1, linestyle='--', linewidth=1.0, color='black')
plt.axvline(day_2, linestyle='--', linewidth=1.0, color='black')
plt.axvline(day_3, linestyle='--', linewidth=1.0, color='black')
plt.plot(t_dates, I1, '-', color='blue', label=u'Level 0 until 31st of July')
plt.plot(t_dates, I6, '-', color='orange', label=u'Level 0 until 31st of July, \n then relax to level 3 until '
                                                 u'end of the year')
plt.plot(t_dates, I2, '-', color='m', label=u'Level 0 until 10th of May, level 1 until 31st of July')
plt.plot(t_dates, I4, '-', color='y', label=u'Level 0 until 10th of May, relax to level 1 \n until 30th of June, '
                                            u'relax to level 2 until 31th of July')
plt.plot(t_dates, I3, '-', color='green', label=u'Level 0 until 10th of May, level 2 until 31st of July')
plt.plot(t_dates, I5, '-', color='red', label=u'Simply open on 10th of May')
plt.plot(t_dates, I0, '-', color='black', label=u'Do-nothing scenario')
plt.xlim([datetime.date(year, month, day), datetime.date(2020, 12, 31)])
ax.xaxis.set_major_locator(months)
ax.xaxis.set_minor_locator(weeks)
ax.xaxis.set_major_formatter(month_fmt)
plt.ylabel(u'individuals')
plt.legend(loc='upper right', fontsize='small', framealpha=0.1, fancybox=True)

save_fig('_scenarioA_infected.svg')

# Casualties
_, ax = plt.subplots(figsize=(9, 6))
plt.axvline(day_1, linestyle='--', linewidth=1.0, color='black')
plt.axvline(day_2, linestyle='--', linewidth=1.0, color='black')
plt.axvline(day_3, linestyle='--', linewidth=1.0, color='black')
plt.plot(t_dates, C1, '-', color='blue', label=u'Level 0 until 31st of July')
plt.plot(t_dates, C6, '-', color='orange', label=u'Level 0 until 31st of July, \n then relax to level 3 until the '
                                                 u'end of the year')
plt.plot(t_dates, C2, '-', color='m', label=u'Level 0 until 10th of May, level 1 until 31st of July')
plt.plot(t_dates, C4, '-', color='y', label=u'Level 0 until 10th of May, relax to level 1 \n until 30th of June, '
                                            u'relax to level 2 until 31st of July')
plt.plot(t_dates, C3, '-', color='green', label=u'Level 0 until 10th of May, level 2 until 31st of July')
plt.plot(t_dates, C5, '-', color='red', label=u'Simply open on 10th of May')
plt.plot(t_dates, C0, '-', color='black', label=u'Do-nothing scenario')
plt.xlim([datetime.date(year, month, day), datetime.date(2020, 12, 31)])
ax.xaxis.set_major_locator(months)
ax.xaxis.set_minor_locator(weeks)
ax.xaxis.set_major_formatter(month_fmt)
plt.ylabel(u'individuals')
plt.legend(loc='lower right', fontsize='small', framealpha=0.1, fancybox=True)

save_fig('_scenarioA_casualties.svg')

_, ax = plt.subplots(figsize=(9, 6))
plt.axvline(day_1, linestyle='--', linewidth=1.0, color='black')
plt.axvline(day_2, linestyle='--', linewidth=1.0, color='black')
plt.axvline(day_3, linestyle='--', linewidth=1.0, color='black')
plt.axhline(leitos, linestyle='--', linewidth=1.0, color='red', label=u'Available ICU beds')
plt.plot(t_dates, H1, '-', color='blue', label=u'Level 0 until 31st of July')
plt.plot(t_dates, H6, '-', color='orange', label=u'Level 0 until 31st of July, \n then relax to level 3 until the '
                                                 u'end of the year')
plt.plot(t_dates, H2, '-', color='m', label=u'Level 0 until 10th of May, level 1 until 31st of July')
plt.plot(t_dates, H4, '-', color='y', label=u'Level 0 until 10th of May, relax to level 1 \n until 30th of June, '
                                            u'relax to level 2 until 31st of July')
plt.plot(t_dates, H3, '-', color='green', label=u'Level 0 until 10th of May, level 2 until 31st of July')
plt.plot(t_dates, H5, '-', color='red', label=u'Simply open on 10th of May')
plt.plot(t_dates, H0, '-', color='black', label=u'Do-nothing scenario')
plt.plot(t_dates, L1, '--', color='blue')
plt.plot(t_dates, L6, '--', color='orange')
plt.plot(t_dates, L2, '--', color='m')
plt.plot(t_dates, L4, '--', color='y')
plt.plot(t_dates, L3, '--', color='green')
plt.plot(t_dates, L5, '--', color='red')
plt.plot(t_dates, L0, '--', color='black')
plt.xlim([datetime.date(year, month, day), datetime.date(2020, 12, 31)])
ax.xaxis.set_major_locator(months)
ax.xaxis.set_minor_locator(weeks)
ax.xaxis.set_major_formatter(month_fmt)
plt.ylabel(u'individuals')
plt.legend(loc='upper right', fontsize='small', framealpha=0.1, fancybox=True)

save_fig('_scenarioA_hosp.svg')

_, ax = plt.subplots(figsize=(9, 6))
plt.axvline(day_1, linestyle='--', linewidth=1.0, color='black')
plt.axvline(day_2, linestyle='--', linewidth=1.0, color='black')
plt.axvline(day_3, linestyle='--', linewidth=1.0, color='black')
plt.axhline(N0[-1]*(1-1/r0_post), linestyle='--', linewidth=1.0, color='red', label='Herd immunity level')
plt.plot(t_dates, R1, '-', color='blue', label=u'Level 0 until 31st of July')
plt.plot(t_dates, R6, '-', color='orange', label=u'Level 0 until 31st of July, \n then relax to level 3 until the '
                                                 u'end of the year')
plt.plot(t_dates, R2, '-', color='m', label=u'Level 0 until 10th of May, level 1 until 31st of July')
plt.plot(t_dates, R4, '-', color='y', label=u'Level 0 until 10th of May, relax to level 1 \n until 30th of June, '
                                            u'relax to level 2 until 31st of July')
plt.plot(t_dates, R3, '-', color='green', label=u'Level 0 until 10th of May, level 2 until 31st of July')
plt.plot(t_dates, R5, '-', color='red', label=u'Simply open on 10th of May')
plt.plot(t_dates, R0, '-', color='black', label=u'Do-nothing scenario')
plt.xlim([datetime.date(year, month, day), datetime.date(2020, 12, 31)])
ax.xaxis.set_major_locator(months)
ax.xaxis.set_minor_locator(weeks)
ax.xaxis.set_major_formatter(month_fmt)
plt.ylabel(u'individuals')
plt.legend(loc='lower right', fontsize='small', framealpha=0.1, fancybox=True)

save_fig('_scenarioA_recovered.svg')

plt.draw()
