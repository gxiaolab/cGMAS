import numpy as np

############
# This function counts the overlapping substring from a string
# http://stackoverflow.com/questions/2970520/string-count-with-overlapping-occurrences
###########
def occurrences(string, sub):
	count = start = 0
	while True:
		start = string.find(sub, start) + 1
		if start > 0:
			count+=1
		else:
			return count

#################
# This function is for pybedtools interval objects
#	=> convert unicode to string under Python 2; all other values pass through unchanged
# http://pythonhosted.org/pybedtools/intervals.html
################
import sys
def show_value(s):
	if sys.version_info.major == 2:
		if isinstance(s, unicode):
			return str(s)
	return s

#################
# Iterating over every N elements in a list
# http://stackoverflow.com/questions/5389507/iterating-over-every-two-elements-in-a-list
# http://stackoverflow.com/questions/4356329/creating-a-python-dictionary-from-a-line-of-text/4356415#4356415
#################
from itertools import izip

def grouped(iterable, n):
    "s -> (s0,s1,s2,...sn-1), (sn,sn+1,sn+2,...s2n-1), (s2n,s2n+1,s2n+2,...s3n-1), ..."
    return izip(*[iter(iterable)]*n)

#for x, y in grouped(l, 2):
#   print "%d + %d = %d" % (x, y, x + y)



##############################
# Sliding window of a list
# size: window size
# step: sliding step
def sliding_window(iterable, size=2, step=1):
	res = []
	lst=list(window(iterable,size))
	for x in range(0,len(lst),step):
		res.append(lst[x])
	return res

##############################
# group window of a list
# size: window size
# https://coderwall.com/p/zvuvmg/sliding-window-in-python
# e.g.:
# >>> list(window([1,2,3,4,5,6,7,8,9],3))
# [[1, 2, 3], [2, 3, 4], [3, 4, 5], [4, 5, 6], [5, 6, 7], [6, 7, 8], [7, 8, 9]]
def window(iterable, size=2):
	i = iter(iterable)
	win = []
	for e in range(0, size):
		win.append(next(i))
	yield win
	for e in i:
		win = win[1:] + [e]
		yield win


##############################################
# monotonically increasing or decreasing
# https://stackoverflow.com/questions/4983258/python-how-to-check-list-monotonicity
def strictly_increasing(L):
    return all(x<y for x, y in zip(L, L[1:]))

def strictly_decreasing(L):
    return all(x>y for x, y in zip(L, L[1:]))

def non_increasing(L):
    return all(x>=y for x, y in zip(L, L[1:]))

def non_decreasing(L):
    return all(x<=y for x, y in zip(L, L[1:]))
##############################################



#############
# Randomly distribute data points to binned x-axis
#discretize the x-axis to bins: https://stackoverflow.com/questions/31730028/how-can-i-generate-a-random-sample-of-bin-counts-given-a-sequence-of-bin-probabi
# * nbin: number of bins
# * tot: total number of data to be distributed into bins
# * prob: probability for each bin; a list w/ length of nbin+1
# * trials: number of times for random generation to get mean and variance; default 100 times
# * targBINindex: bin of interest to calc mean & var
# return: bg = list of randome bg result for targBINindex bin
#	  bgpeak = mean of bg + var of bg
def distributeBins(nbin,tot,prob,targBINindex,trials=100):
	cdf=np.cumsum(prob)
	bg = []
	for x in range(trials):
		counts = np.histogram(np.random.rand(tot), bins=cdf)[0]
		bg.append(counts[targBINindex])
	bgpeak = np.mean(bg)+np.var(bg)
	return bg, bgpeak

