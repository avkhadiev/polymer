import numpy as np
import matplotlib.pyplot as plt

sim='20_2'

# NAMES
# R - euclidian separation
# rr - euclidian separation per atom
# sl - SLERP
# p1 - PLERP RATTLE-based
# p2 - PLERP simplified
# l - contour length
# ln - contour length normalized
# c - convergence status (1 = converged)
# f - file
# a - numpy array
# fr - success rate
# avg - mean
# std - standard deviation

# read data
R_f    = sim + "_R.txt"
rr_f   = sim + "_rr.txt"
md_f   = sim + "_md.txt"
sl_c_f = sim + "_sl_c.txt"
sl_l_f = sim + "_sl_l.txt"
p1_c_f = sim + "_p1_c.txt"
p1_l_f = sim + "_p1_l.txt"
p2_c_f = sim + "_p2_c.txt"
p2_l_f = sim + "_p2_l.txt"
R_a    = np.loadtxt(R_f)
rr_a   = np.loadtxt(rr_f)
md_l_a = np.loadtxt(md_f)
sl_c_a = np.loadtxt(sl_c_f)
sl_l_a = np.loadtxt(sl_l_f)
p1_c_a = np.loadtxt(p1_c_f)
p1_l_a = np.loadtxt(p1_l_f)
p2_c_a = np.loadtxt(p2_c_f)
p2_l_a = np.loadtxt(p2_l_f)
c_a = [sl_c_a, p1_c_a, p2_c_a]
l_a = [sl_l_a, p1_l_a, p2_l_a]
print "Files read."

# calculate failure rate
sl_fr = float(np.count_nonzero(sl_c_a == 0)) / np.size(sl_c_a)
p1_fr = float(np.count_nonzero(p1_c_a == 0)) / np.size(p1_c_a)
p2_fr = float(np.count_nonzero(p2_c_a == 0)) / np.size(p2_c_a)
print "%s %5.5f" % ("SLERP  failure rate", sl_fr)
print "%s %5.5f" % ("PLERP1 failure rate", p1_fr)
print "%s %5.5f" % ("PLERP2 failure rate", p2_fr)
# calculate the mean and stdev of euclidian separation
R_avg = np.mean(R_a)
R_std = np.std(R_a)
print "%s %5.5f" % ("Euclidian separation mean ", R_avg)
print "%s %5.5f" % ("Euclidian separation stdev", R_std)
# calculate normalized contour lengths
md_ln_a = md_l_a / R_a
sl_ln_a = sl_l_a / R_a
p1_ln_a = p1_l_a / R_a
p2_ln_a = p2_l_a / R_a
# remove values that did not converge
sl_ln_a[sl_c_a == 0] = np.nan
p1_ln_a[p1_c_a == 0] = np.nan
p2_ln_a[p2_c_a == 0] = np.nan
# calculate the mean and stdev of normalized length of converged paths
md_ln_avg = np.mean(md_ln_a)
sl_ln_avg = np.nanmean(sl_ln_a)
p1_ln_avg = np.nanmean(p1_ln_a)
p2_ln_avg = np.nanmean(p2_ln_a)
md_ln_std = np.std(md_ln_a)
sl_ln_std = np.nanstd(sl_ln_a)
p1_ln_std = np.nanstd(p1_ln_a)
p2_ln_std = np.nanstd(p2_ln_a)
print "%s %5.5f" % ("SLERP  mean  converged path normalized contour length", sl_ln_avg)
print "%s %5.5f" % ("PLERP1 mean  converged path normalized contour length", p1_ln_avg)
print "%s %5.5f" % ("PLERP2 mean  converged path normalized contour length", p2_ln_avg)
print "%s %5.5f" % ("MD     mean  converged path normalized contour length", md_ln_avg)
print "%s %5.5f" % ("SLERP  stdev converged path normalized contour length", sl_ln_std)
print "%s %5.5f" % ("PLERP1 stdev converged path normalized contour length", p1_ln_std)
print "%s %5.5f" % ("PLERP2 stdev converged path normalized contour length", p2_ln_std)
print "%s %5.5f" % ("MD     stdev converged path normalized contour length", md_ln_std)
# make histogram
plt.hist(md_ln_a,                     histtype='step', label = 'MD')
plt.hist(p1_ln_a[~np.isnan(p1_ln_a)], histtype='step', label = 'PLERP 1')
plt.hist(p2_ln_a[~np.isnan(p2_ln_a)], histtype='step', label = 'PLERP 2')
plt.legend()
plt.show()
