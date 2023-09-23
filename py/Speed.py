import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter, LogLocator)

colors = plt.rcParams["axes.prop_cycle"].by_key()["color"]

LOWER_LIM = 1e-13

fname = "data/Speed.txt"
dataf = open(fname, "r")
data = {}
speeds = []
precs = []
for line in dataf.readlines():
	line = line.split(" ")
	speed = float(line[1])
	speed *= 1e6 # rescale to microseconds
	prec = float(line[2])
	data[line[0]] = [speed, prec]
	speeds.append(speed)
	if prec > LOWER_LIM:
		precs.append(prec)
dataf.close()

names = np.array(list(data.keys()))

# group in different categories
# exact measurements have arrows pointing down for precision
exacts = ["ZS", "Diag", "JP"]
exact_mask = np.array([name in exacts for name in names])
exact_speeds = np.array([data[name][0] for name in names[exact_mask]])
exact_precs = np.array([data[name][1] for name in names[exact_mask]])

# these expressions are at least 1% precise
precises = [r"DMP^0", r"DMP^1", "MP", r"AM^2", r"AM^{5/2}", "MF", "AKT", r"DPZ^0", r"DPZ^2", "OS", "Exp"]
precise_mask = np.array([name in precises for name in names])
precise_speeds = np.array([data[name][0] for name in names[precise_mask]])
precise_precs = np.array([data[name][1] for name in names[precise_mask]])

# don't plot these expressions
ignores = ["Zeroth"]
# remaining expressions are fairly imprecise
imprecise_mask = np.array([name not in exacts + precises + ignores for name in names])
imprecises = names[imprecise_mask]
imprecise_speeds = np.array([data[name][0] for name in names[imprecise_mask]])
imprecise_precs = np.array([data[name][1] for name in names[imprecise_mask]])

# actually plot the (im)precise expressions
dpz_mask = np.array(["DPZ" in x for x in names[precise_mask]])
plt.plot(precise_speeds[dpz_mask], precise_precs[dpz_mask], ".", c = colors[1])
plt.plot(precise_speeds[dpz_mask == False], precise_precs[dpz_mask == False], ".", c = colors[0])
#plt.plot(imprecise_speeds, imprecise_precs, "r.")

# log-log
plt.xscale("log")
plt.yscale("log")

# label the axes
plt.xlabel(r"$t{\rm\ [}\mu{\rm s]}$")
plt.ylabel(r"$|\Delta P|/P$")

# set the axes
v = list(plt.axis())
v[0] = 2e-2
v[1] = 1e-0
v[2] = LOWER_LIM
v[3] = 1e-2
plt.axis(v)

# reasonable looking xticks
xticks = [v[0], 1e-1, 1e-0]
plt.xticks(xticks, [r"$%g$" % xtick for xtick in xticks])

# yticks
plt.gca().yaxis.set_minor_locator(LogLocator(base = 10, numticks = 20))
plt.gca().yaxis.set_minor_formatter(FormatStrFormatter(""))

# draw the arrows for the exacts
scale = 2.5
for exact in exact_speeds:
	plt.annotate("", xy = (exact, v[2]), xytext = (exact, scale * v[2]), arrowprops = dict(fc = "black", shrink = 0., width = 1., headwidth = 6., headlength = 8.))

# determine if the label will be to the right or the left, so the text doesn't all overlap
fs = 12
shift = 1.03
va = "center"

rights = ["AJLOS(31)", "AJLOS(48)", "FL", "Vacuum", r"AM^{5/2}", "MF", r"AM^2"]
bottoms = ["MP"]
has = {} # horizontal alignments
vas = {}
shifts = {}
color = {}
for name in names:
	if "DPZ" in name:
		color[name] = colors[1]
	else:
		color[name] = colors[0]
	if name in rights:
		has[name] = "left"
		vas[name] = "center"
		shifts[name] = shift
	elif name in bottoms:
		has[name] = "center"
		vas[name] = "top"
		shifts[name] = 1.
	else:
		has[name] = "right"
		vas[name] = "center"
		shifts[name] = 1. / shift

# label each dot
for name in precises: plt.gca().text(data[name][0] * shifts[name], data[name][1], r"${\rm %s}$" % name, fontsize = fs, va = vas[name], ha = has[name], color = color[name])
#for name in imprecises: plt.gca().text(data[name][0] * shifts[name], data[name][1], r"${\rm %s}$" % name, fontsize = fs, va = vas[name], ha = has[name], color = "r")
for name in exacts: plt.gca().text(data[name][0], shift * scale * v[2], r"${\rm %s}$" % name, fontsize = fs, va = "bottom", ha = "center", color = "k")

# which corner is best and which is worst
theta = 36.
fs = 12
plt.text(0.03, 0.05, r"${\rm BEST}$", ha = "left", va = "bottom", transform = plt.gca().transAxes, rotation = theta, fontsize = fs)
plt.text(0.96, 0.995, r"${\rm WORST}$", ha = "right", va = "top", transform = plt.gca().transAxes, rotation = theta, fontsize = fs)

dx1 = 0.01
dx2 = 0.12
arrowprops = dict(fc = "black", shrink = 0., width = 2., headwidth = 10., headlength = 12.)
plt.annotate(r"", xy = (dx1, dx1), xycoords = "axes fraction", xytext = (dx2, dx2), textcoords = "axes fraction", arrowprops = arrowprops)
plt.annotate(r"", xy = (1 - dx1, 1 - dx1), xycoords = "axes fraction", xytext = (1 - dx2, 1 - dx2), textcoords = "axes fraction", arrowprops = arrowprops)

# version number
plt.text(0.10, 0.01, r"${\rm Nu}$-${\rm Pert}$-${\rm Compare\ v1.2}$", ha = "left", va = "bottom", transform = plt.gca().transAxes, color = "gray", fontsize = 10)
# save the figure
# plt.savefig("fig/Speed.pdf")
plt.show()

