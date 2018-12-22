import matplotlib
import numpy
import matplotlib.pyplot

file1 = "DGV_negative_gain_sizes.csv"
file1label = "negativeDGVgains"

file2 = "pathogenic_gains.csv"
file2label = "positiveDECIPHERgains"


fig, ax = matplotlib.pyplot.subplots()
load1 = numpy.loadtxt(open(file1, "rb"), delimiter=",", skiprows=1)
load2 = load1[load1<50000]
load3 = load2[load2>2000]

load4 = numpy.loadtxt(open(file2, "rb"), delimiter=",", skiprows=1)
load5 = load4[load4<50000]
load6 = load5[load5>2000]

matplotlib.pyplot.hist(load3,log=False,bins=25, alpha=0.5, label=file1label)
matplotlib.pyplot.hist(load6,log=False,bins=25, alpha=0.5, label=file2label)

matplotlib.pyplot.legend(loc='upper right')


ax.set_xlabel("Number of base pairs")
ax.set_ylabel("Frequency")
ax.set_title("CNV Sizes")

matplotlib.pyplot.show()
