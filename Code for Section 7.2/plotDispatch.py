# -*- coding: utf-8 -*-

import numpy as np
import pandas as pd
import re
import random
import math
import pickle
import time
import itertools
from scipy import sparse
from scipy.special import comb
from scipy.sparse.linalg import spsolve
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter

data = pd.read_csv('Table-OptDispatchPolicy.csv')
data1 = pd.read_csv('Table-OptLocation-9Units-6.0threshould.csv')
data2 = pd.read_csv('Table-OptLocation-9Units-9.0threshould.csv')
data_9units = pd.merge(data1, data2[['location', 'frac_total_7.5']], on='location', how='outer')

data1 = pd.read_csv('Table-OptLocation-10Units-6.0threshould.csv')
data2 = pd.read_csv('Table-OptLocation-10Units-9.0threshould.csv')
data_10units = pd.merge(data1, data2[['location', 'frac_total_7.5']], on='location', how='outer')

data1 = pd.read_csv('Table-OptLocation-11Units-6.0threshould.csv')
data2 = pd.read_csv('Table-OptLocation-11Units-9.0threshould.csv')
data_11units = pd.merge(data1, data2[['location', 'frac_total_7.5']], on='location', how='outer')

MRT_dispatch = data['MRT'] + 1.5
MRT_9units = data_9units['MRT'] + 1.5
MRT_10units = data_10units['MRT'] + 1.5
MRT_11units = data_11units['MRT'] + 1.5
F6_9units = data_9units['frac_total_4.5']
F6_10units = data_10units['frac_total_4.5']
F6_11units = data_11units['frac_total_4.5']
F9_9units = data_9units['frac_total_7.5']
F9_10units = data_10units['frac_total_7.5']
F9_11units = data_11units['frac_total_7.5']

# Plot Empirical PDF
# plt.figure(figsize=(12, 5))

# plt.subplot(1, 2, 1)
# plt.hist(MRT_dispatch, bins=20, alpha=0.6, color='black', edgecolor='black')
# # plt.title('Empirical PDF')
# plt.xlabel('MRT')
# plt.ylabel('Number of Alternatives')

# Plot Empirical CDF
plt.figure()
sorted_data = np.sort(MRT_dispatch)
cdf = np.arange(1, len(sorted_data)+1) / len(sorted_data)
plt.plot(sorted_data, cdf, marker='.', markersize=8, linestyle='none', label='Policy Performance', color='grey')
plt.step(sorted_data, cdf, where='mid', color='black', label='Empirical CDF')
plt.xlabel('MRT (min)', fontsize=10)
plt.ylabel('Cumulative Distribution', fontsize=10)
bound = np.searchsorted(sorted_data, MRT_dispatch[0], side='right') - 1
xc = sorted_data[bound]
yc = cdf[bound]
# plt.axvline(x=xc, color='grey', linestyle='-.', linewidth=1.0, ymin=0, ymax=yc, clip_on=True)
x_limits = plt.xlim()
x_min, x_max = x_limits
relative_xc = (xc - x_min + 0.05) / (x_max - x_min)
plt.axhline(y=yc, color='dimgrey', linestyle='-.', linewidth=0.7, xmin=0, xmax=relative_xc, clip_on=True)
relative_xc = (sorted_data[0] - x_min + 0.05) / (x_max - x_min)
plt.axhline(y=cdf[0], color='dimgrey', linestyle='-.', linewidth=0.7, xmin=0, xmax=relative_xc, clip_on=True)
current_yticks = plt.yticks()[0]
plt.yticks(list(current_yticks) + [yc])
# plt.plot([xc, xc], [0, yc], color='lightgrey', linestyle='--', linewidth=0.7)
# plt.plot([sorted_data[0] - 0.1, xc], [yc, yc], color='grey', linestyle='-.', linewidth=1.0)
plt.plot(xc, yc, 'bo', markersize=6)  # Add a black dot at the point (xc, yc)
plt.plot(sorted_data[0], cdf[0], 'm*', markersize=8)

# Add an arrow pointing to the point
plt.annotate('Nearest Available Policy, 5.887', xy=(xc, yc), xytext=(xc + 0.1, yc - 0.01),
             arrowprops=dict(arrowstyle='->', facecolor='black', shrinkB=5),
             fontsize=10, color='black')
plt.annotate('Best Policy, 5.885', xy=(sorted_data[0], cdf[0]), xytext=(sorted_data[0] + 0.1, cdf[0] - 0.01),
             arrowprops=dict(arrowstyle='->', facecolor='black', shrinkB=5),
             fontsize=10, color='black')

# # Annotate the point with its coordinates
# plt.text(xc + 0.1, yc - 0.05, f'({xc:.2f}, {yc:.2f})', ha='left', color='black', fontsize=10)
# current_xticks = plt.xticks()[0]
# current_yticks = plt.yticks()[0]
# new_xticks = [tick for tick in current_xticks if tick != 5.9]
# new_yticks = [tick for tick in current_yticks if tick != 0.6]
# new_xticks.append(xc)
# new_yticks.append(yc)
# plt.xticks(new_xticks)
# plt.yticks(new_yticks)
# formatter = FuncFormatter(lambda x, _: f'{x:.2f}')
# plt.gca().xaxis.set_major_formatter(formatter)
# plt.gca().yaxis.set_major_formatter(formatter)
# plt.yticks(list(plt.yticks()[0]) + [yc])
# plt.xticks(list(plt.xticks()[0]) + [xc])
# plt.text(plt.gca().get_xlim()[0] - 0.05, yc, f'{yc:.2f}',
#          ha='right', va='center', color='black', fontsize=10)
# plt.text(xc, plt.gca().get_ylim()[0] - 0.05, f'{xc:.2f}',
#          ha='center', va='top', color='black', fontsize=10)
plt.xlim(sorted_data[0] - 0.05, sorted_data[-1]+0.05)
plt.ylim(-0.05, 1.1)
plt.xticks(fontsize=10)
plt.yticks(fontsize=10)
plt.legend(bbox_to_anchor=(1,0.8))
# plt.xlim(plt.gca().get_xlim()[0], plt.gca().get_xlim()[1] + 0.1)
plt.savefig('dispatch.png', transparent=False, bbox_inches='tight', dpi=400)
plt.tight_layout()
plt.show()

# # Plot 9-11 MRT distribution
# fig, ax = plt.subplots()
# sorted_data = np.sort(MRT_9units)
# cdf = np.arange(1, len(sorted_data)+1) / len(sorted_data)
# ax.plot(sorted_data, cdf, marker='.', linestyle='none', label='9 units')
# ax.step(sorted_data, cdf, where='mid', color='b')
# sorted_data = np.sort(MRT_10units)
# cdf = np.arange(1, len(sorted_data)+1) / len(sorted_data)
# plt.plot(sorted_data, cdf, marker='.', linestyle='none', label='10 Units')
# plt.step(sorted_data, cdf, where='mid', color='orange')
# sorted_data = np.sort(MRT_11units)
# cdf = np.arange(1, len(sorted_data)+1) / len(sorted_data)
# plt.plot(sorted_data, cdf, marker='.', linestyle='none', label='11 Units')
# plt.step(sorted_data, cdf, where='mid', color='green')
# plt.legend()
# plt.xlabel('MRT (min)')
# plt.ylabel('Cumulative Distribution')
# plt.savefig('MRT-CDF9-11.png', transparent=False, bbox_inches='tight', dpi=400)
# plt.tight_layout()
# plt.show()
#
#
# # Plot 9-11 F6 distribution
# fig, ax = plt.subplots()
# sorted_data = np.sort(F6_9units)
# cdf = np.arange(1, len(sorted_data)+1) / len(sorted_data)
# ax.plot(sorted_data, cdf, marker='.', linestyle='none', label='9 units')
# ax.step(sorted_data, cdf, where='mid', color='b')
# sorted_data = np.sort(F6_10units)
# cdf = np.arange(1, len(sorted_data)+1) / len(sorted_data)
# plt.plot(sorted_data, cdf, marker='.', linestyle='none', label='10 Units')
# plt.step(sorted_data, cdf, where='mid', color='orange')
# sorted_data = np.sort(F6_11units)
# cdf = np.arange(1, len(sorted_data)+1) / len(sorted_data)
# plt.plot(sorted_data, cdf, marker='.', linestyle='none', label='11 Units')
# plt.step(sorted_data, cdf, where='mid', color='green')
# plt.legend()
# plt.xlabel('$F(6)$')
# plt.ylabel('Cumulative Distribution')
# plt.savefig('F6-CDF9-11.png', transparent=False, bbox_inches='tight', dpi=400)
# plt.tight_layout()
# plt.show()
# #
# # Plot 9-11 F9 distribution
# fig, ax = plt.subplots()
# sorted_data = np.sort(F9_9units)
# cdf = np.arange(1, len(sorted_data)+1) / len(sorted_data)
# ax.plot(sorted_data, cdf, marker='.', linestyle='none', label='9 units')
# ax.step(sorted_data, cdf, where='mid', color='b')
# sorted_data = np.sort(F9_10units)
# cdf = np.arange(1, len(sorted_data)+1) / len(sorted_data)
# plt.plot(sorted_data, cdf, marker='.', linestyle='none', label='10 Units')
# plt.step(sorted_data, cdf, where='mid', color='orange')
# sorted_data = np.sort(F9_11units)
# cdf = np.arange(1, len(sorted_data)+1) / len(sorted_data)
# plt.plot(sorted_data, cdf, marker='.', linestyle='none', label='11 Units')
# plt.step(sorted_data, cdf, where='mid', color='green')
# plt.legend()
# plt.xlabel('$F(9)$')
# plt.ylabel('Cumulative Proportion')
# plt.savefig('F9-CDF9-11.png', transparent=False, bbox_inches='tight', dpi=400)
# plt.tight_layout()
# plt.show()