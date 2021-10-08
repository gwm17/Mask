#!/usr/bin/env python3

import matplotlib.pyplot as plt
import numpy as np
import pickle
import sys

def SetManager(figure):
	dummy = plt.figure()
	manager = dummy.canvas.manager
	manager.canvas.figure = figure
	figure.set_canvas(manager.canvas)

def ViewPyPlots(filename):

	figure = pickle.load(open(filename, "rb"))
	SetManager(figure)

	figure.set_size_inches(12, 12)
	plt.tight_layout()

	plt.show(block=True)


def main():
	if len(sys.argv) == 2:
		ViewPyPlots(sys.argv[1])
	else:
		print("Unable to run ViewPyPlots, incorrect number of commandline arguments -- requires an input pickle file")

if __name__ == '__main__':
	main()