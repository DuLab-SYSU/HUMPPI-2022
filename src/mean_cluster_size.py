import os, sys
import math
import matplotlib.pyplot as plt


def read_(_out):
	data = []
	with open(_out) as f:
		for line in f:
			c_i = line.strip().split()
			data.append(c_i)
	return data


def cal_mcs(_out):
	clusters = read_(_out)
	c = [len(x) for x in clusters]
	n = sum(c)
	try:
		mcs = sum([c_i * c_i for c_i in c]) / float(n * n)
	except ZeroDivisionError:
		print("Error: %s" % _out)
	print("# of clusters: %s" % len(clusters))
	print("# of nodes: %s" % n)
	print("Mean Cluster Size: %s" % mcs)
	return mcs


def main(file_base):
	results = []
	file_base = file_base
	for i in range(101, 1501):
		print("Inflation: %s" % str(i/100))
		_out = "batch_MCL_out/out." + file_base + ".I" + str(i)
		mcs = cal_mcs(_out)
		results.append((i/100, mcs))

	x_, y_ = zip(*results)
	plt.plot(x_, y_, marker='D')
	plt.show()

if __name__ == "__main__":
	file_base = sys.argv[1]
	main(file_base)
