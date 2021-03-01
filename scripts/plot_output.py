# coding: utf-8
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import numpy as np
import csv
import argparse
import datetime
import os

##### Process command line options
##### Variable parameters, for error estimation within reasonable bounds
parser = argparse.ArgumentParser(description=u'This script plots the optimization results for any scenario.')
parser.add_argument('-d', '--day', type=int, nargs=4, help='Days of measure beginning - four values required',
                    required=True)
parser.add_argument('-s', '--scale_factor', type=float, help='Scale factor accounting for under notification ',
                    required=True)
parser.add_argument('-o', '--output_file', help='Output file with results', required=True)
parser.add_argument('-sc', '--scenario', help='Scenario folder', required=True)
parser.add_argument('-i', '--init_date', type=int, nargs=3, help='Date of the first day of notification [dat], [month],'
                                                                 '[year]', required=True)

args = parser.parse_args()

# Time specifications for plot
months = mdates.MonthLocator()  # every month
weeks = mdates.WeekdayLocator()  # every week
month_fmt = mdates.DateFormatter('%b')

s_0 = float(args.scale_factor)

# Limite entre cenários
day_init = int(args.day[0])
day_next_1 = int(args.day[1])
day_next_2 = int(args.day[2])
day_next_3 = int(args.day[3])
day = int(args.init_date[0])
month = int(args.init_date[1])
year = int(args.init_date[2])

t_days = 400
age_strata = 16
compartments = 11

output_file = args.output_file

# Notification data file
notification_file = 'notification.csv'

# Scenario folder
cenario_folder = args.scenario

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
            cData = Row[1:]
            j = j + 1
        elif j == 4:
            day_one = int(Row[1])
            j = j + 1
        elif j == 5:
            month_one = int(Row[1])
            j = j + 1
        elif j == 6:
            year_one = int(Row[1])
            j = j + 1
        elif j == 7:
            leitos = int(Row[1])
            print(leitos)
            j = j + 1

YData = np.array([xx for xx in YData if xx], dtype=np.float64)
tData = np.array([xx for xx in tData if xx], dtype=np.int64)
CData = np.array([xx for xx in CData if xx], dtype=np.float64)
cData = np.array([xx for xx in cData if xx], dtype=np.int64)

# Return to script directory
os.chdir("..")
os.chdir("..")
os.chdir("..")
os.chdir("scripts")

YData = s_0 * YData

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
t_array = np.array(t_days, dtype=np.float64)

# Leitura dos arquivos de dados fundamentais
curr_dir = os.path.abspath(os.curdir)
print(u'Diretório atual ' + curr_dir)
print(u'Movendo para o diretório de saída (output) ')
os.chdir("..")
os.chdir("output")

curr_dir = os.path.abspath(os.curdir)
print(u'Diretório de saída (output) ' + curr_dir)

with open(output_file, "r") as csvfile:
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

print(u'Voltando para o diretório de script')
os.chdir("..")
os.chdir("scripts")

# Soma todas as faixas etárias

S = np.sum(s, axis=1)
E = np.sum(e, axis=1)
I = np.sum(y, axis=1)
R = np.sum(r, axis=1)
N = np.sum(n, axis=1)
A = np.sum(a, axis=1)
C = np.sum(c, axis=1)
H = np.sum(h, axis=1)
L = np.sum(l, axis=1)
RI = np.sum(ri, axis=1)
Ac = C + RI + I + H

# Criam-se grupos de faixa etária: 1: 0 - 20 / 2: 20 - 55 / 3: 55 ou mais

I1 = np.sum(y[:, 0:3], axis=1)
I2 = np.sum(y[:, 4:9], axis=1)
I3 = np.sum(y[:, 10:16], axis=1)

R1 = np.sum(r[:, 0:3], axis=1)
R2 = np.sum(r[:, 4:9], axis=1)
R3 = np.sum(r[:, 10:16], axis=1)

H1 = np.sum(h[:, 0:3], axis=1)
H2 = np.sum(h[:, 4:9], axis=1)
H3 = np.sum(h[:, 10:16], axis=1)

L1 = np.sum(l[:, 0:3], axis=1)
L2 = np.sum(l[:, 4:9], axis=1)
L3 = np.sum(l[:, 10:16], axis=1)

C1 = np.sum(c[:, 0:3], axis=1)
C2 = np.sum(c[:, 4:9], axis=1)
C3 = np.sum(c[:, 10:16], axis=1)

# Plota os gráficos
t_array = np.linspace(0, t_days - 1, t_days, dtype=int)
tempDate = datetime.datetime(year, month, day)
t_dates = (tempDate + t_array * datetime.timedelta(days=1))
data_tData = (tempDate + (tData - 1) * datetime.timedelta(days=1))
data_cData = (tempDate + (cData - 1) * datetime.timedelta(days=1))
day_1 = (tempDate + day_next_1 * datetime.timedelta(days=1))
day_2 = (tempDate + day_next_2 * datetime.timedelta(days=1))
day_3 = (tempDate + day_next_3 * datetime.timedelta(days=1))
day_last = datetime.date(2020, 12, 31)

plt.figure(1)
ax = plt.subplot(121)
plt.plot(t_dates, S, '-', label='S')
plt.plot(t_dates, I, '-', label='I')
plt.plot(t_dates, R, '-', label='R')
plt.plot(t_dates, C, '-', label='C')
plt.plot(t_dates, H, '-', label='H')
plt.plot(t_dates, L, '-', label='L')
plt.plot(t_dates, Ac, '-', label='Ac')
plt.plot(t_dates, N, '-', label='N')
plt.plot(t_dates, E, '-', label='E')
plt.plot(t_dates, A, '-', label='A')
plt.axvspan(day_1, day_last, alpha=0.1, color='blue')
plt.xlim([datetime.date(year, month, day), datetime.date(2020, 12, 31)])
ax.xaxis.set_major_locator(months)
ax.xaxis.set_minor_locator(weeks)
ax.xaxis.set_major_formatter(month_fmt)
plt.xlabel(u'days')
plt.ylabel(u'individuals')
plt.legend(loc='center left', shadow=True, fontsize='small')

ax = plt.subplot(122)
plt.semilogy(t_dates, S)
plt.semilogy(t_dates, I)
plt.semilogy(t_dates, R)
plt.semilogy(t_dates, C)
plt.semilogy(t_dates, H)
plt.semilogy(t_dates, L)
plt.semilogy(t_dates, Ac)
plt.semilogy(data_tData, YData, 'ok')
plt.semilogy(data_cData, CData, 'or', label=u'Casualties')
plt.axhline(leitos, linestyle='--', linewidth=1.0, color='red')
plt.xlim([datetime.date(year, month, day), datetime.date(2020, 12, 31)])
ax.xaxis.set_major_locator(months)
ax.xaxis.set_minor_locator(weeks)
ax.xaxis.set_major_formatter(month_fmt)
plt.ylim([1, 1.1 * N.max()])
plt.xlabel('time (days)')
plt.axvspan(day_1, day_last, alpha=0.1, color='blue')
max_H = np.max(H)
max_L = np.max(L)
max_I = np.max(I)
max_C = np.max(C)
t_max_I = t_array[np.where(I == max_I)]
t_max_L = t_array[np.where(L == max_L)]
if np.max(L) > leitos:
    t_colap = np.min(t_array[np.where(L > leitos)])
    textstr = '\n'.join([r'$Max(H)=%.2e$' % (max_H,), r'$Max(L)=%.2e$' % (max_L,), r'$Max(I)=%.2e$' % (max_I,),
                         r'$t(max(I))=%.f$ days' % (t_max_I,),
                         r'$t(max(L))=%.f$ days' % (t_max_L,),
                         r'$t(breakdown)=%.f$ days' % (t_colap,),
                         r'Casualties $=%.2e$' % (max_C,)])
else:
    textstr = '\n'.join([r'$Max(H)=%.2e$' % (max_H,), r'$Max(L)=%.2e$' % (max_L,), r'$Max(I)=%.2e$' % (max_I,),
                         r'$t(max(I))=%.f$ days' % (t_max_I,),
                         r'$t(max(L))=%.f$ days' % (t_max_L,),
                         r'$t(breakdown)=$ inf days',
                         r'Casualties $=%.2e$' % (max_C,)])

props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)

# place a text box in upper left in axes coords
plt.text(0.68, 0.15, textstr, transform=plt.gca().transAxes, fontsize='small', verticalalignment='center', bbox=props)
plt.suptitle(cenario_folder[7:])

plt.figure(2)
ax = plt.subplot(121)
plt.plot(t_dates, I1, '-b', label='I: 0 - 20')
plt.plot(t_dates, I2, '-g', label='I: 20 - 55')
plt.plot(t_dates, I3, '-r', label='I: 55 or more')
plt.axvspan(day_1, day_last, alpha=0.1, color='blue')
plt.xlim([datetime.date(year, month, day), datetime.date(2020, 12, 31)])
ax.xaxis.set_major_locator(months)
ax.xaxis.set_minor_locator(weeks)
ax.xaxis.set_major_formatter(month_fmt)
plt.xlabel(u'days')
plt.ylabel(u'individuals')
plt.legend(loc='center right', shadow=True, fontsize='small')

ax = plt.subplot(122)
plt.semilogy(t_dates, I1, '-b', label='I: 0 - 20')
plt.semilogy(t_dates, I2, '-g', label='I: 20 - 55')
plt.semilogy(t_dates, I3, '-r', label='I: 55 or more')
plt.semilogy(t_dates, R1, '--b', label='R: 0 - 20')
plt.semilogy(t_dates, R2, '--g', label='R: 20 - 55')
plt.semilogy(t_dates, R3, '--r', label='R: 55 or more')
plt.axvspan(day_1, day_last, alpha=0.1, color='blue')
plt.xlabel('time (days)')
plt.xlim([datetime.date(year, month, day), datetime.date(2020, 12, 31)])
ax.xaxis.set_major_locator(months)
ax.xaxis.set_minor_locator(weeks)
ax.xaxis.set_major_formatter(month_fmt)
plt.ylim([1, 1.1 * N.max()])
plt.legend(loc='lower right', shadow=True, fontsize='small')

plt.suptitle(cenario_folder[7:])

plt.figure(3)
ax = plt.subplot(121)
plt.plot(t_dates, H1, '-b', label='H: 0 - 20')
plt.plot(t_dates, H2, '-g', label='H: 20 - 55')
plt.plot(t_dates, H3, '-r', label='H: 55 or more')
plt.plot(t_dates, L1, '--b', label='L: 0 - 20')
plt.plot(t_dates, L2, '--g', label='L: 20 - 55')
plt.plot(t_dates, L3, '--r', label='L: 55 or more')
plt.axvspan(day_1, day_last, alpha=0.1, color='blue')
plt.xlim([datetime.date(year, month, day), datetime.date(2020, 12, 31)])
ax.xaxis.set_major_locator(months)
ax.xaxis.set_minor_locator(weeks)
ax.xaxis.set_major_formatter(month_fmt)
plt.xlabel(u'days')
plt.ylabel(u'individuals')
plt.legend(loc='center right', shadow=True, fontsize='small')

ax = plt.subplot(122)
plt.semilogy(t_dates, H1, '-b', label='H: 0 - 20')
plt.semilogy(t_dates, H2, '-g', label='H: 20 - 55')
plt.semilogy(t_dates, H3, '-r', label='H: 55 or more')
plt.semilogy(t_dates, L1, '--b', label='L: 0 - 20')
plt.semilogy(t_dates, L2, '--g', label='L: 20 - 55')
plt.semilogy(t_dates, L3, '--r', label='L: 55 or more')
plt.semilogy(t_dates, C1, '-.b', label='C: 0 - 20')
plt.semilogy(t_dates, C2, '-.g', label='C: 20 - 55')
plt.semilogy(t_dates, C3, '-.r', label='C: 55 or more')
plt.axvspan(day_1, day_last, alpha=0.1, color='blue')
plt.xlabel('time (days)')
plt.xlim([datetime.date(year, month, day), datetime.date(2020, 12, 31)])
ax.xaxis.set_major_locator(months)
ax.xaxis.set_minor_locator(weeks)
ax.xaxis.set_major_formatter(month_fmt)
plt.ylim([1, 1.1 * I.max()])
plt.legend(loc='lower right', shadow=True, fontsize='small')

plt.suptitle(cenario_folder[7:])

_, ax = plt.subplots(figsize=(5, 5))
plt.semilogy(t_dates, S, label='Susceptible')
plt.semilogy(t_dates, I, label='Infected')
plt.semilogy(t_dates, C, 'r', label='Casualties')
plt.semilogy(t_dates, Ac, 'k', label='Accumulated cases')
plt.semilogy(data_tData, YData, 'ok', label='Case data')
plt.semilogy(data_cData, CData, 'or', label='Casualties data')
plt.xlim([datetime.date(year, month, day), datetime.date(2020, 12, 31)])
ax.xaxis.set_major_locator(months)
ax.xaxis.set_minor_locator(weeks)
ax.xaxis.set_major_formatter(month_fmt)
plt.ylim([1, 1.1 * N.max()])
plt.xlabel('time (days)')
plt.axvspan(day_1, day_last, alpha=0.1, color='blue')
plt.legend(loc='lower right', shadow=True, fontsize='small')
plt.suptitle(cenario_folder[7:])

ax.set_aspect('equal')

os.chdir("..")
my_path = os.getcwd()
dirs = os.path.join(my_path, 'figures', cenario_folder)
try:
    os.makedirs(dirs)
except OSError:
    pass
plt.savefig(os.path.join(dirs, cenario_folder[:] + '_fig.svg'), format='svg',
            dpi=300, bbox_inches='tight')
os.chdir("scripts")

plt.draw()
