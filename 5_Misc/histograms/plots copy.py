import matplotlib
import numpy
import matplotlib.pyplot

matplotlib.pyplot.style.use("ggplot")
file1 = "negatives.csv"
file1label = "non-pathogenic"

file2 = "positives.csv"
file2label = "pathogenic"


fig, ax = matplotlib.pyplot.subplots()
load1 = numpy.loadtxt(open(file1, "rb"), delimiter=",", skiprows=1)
load1 = load1[:5000]
load2 = load1[load1<4000000]
load3 = load2[load2>2000]

load4 = numpy.loadtxt(open(file2, "rb"), delimiter=",", skiprows=1)
load4=load4[:5000]
load5 = load4[load4<4000000]
load6 = load5[load5>2000]
matplotlib.pyplot.hist(load5,log=False,bins=60, alpha=0.8, label=file2label, color='r')
matplotlib.pyplot.hist(load2,log=False,bins=60, alpha=0.5, label=file1label, color='b')


matplotlib.pyplot.legend(loc='upper right')
#matplotlib.pyplot.xlim(1,1000000)

ax.set_xlabel("Number of base pairs")
ax.set_ylabel("Frequency")
ax.set_title("CNV Sizes")

matplotlib.pyplot.show()
