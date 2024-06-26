# Read some data
# We do not need to close all
import pyedflib
import numpy as np
f = pyedflib.EdfReader("data/test_generator.edf")
n = f.signals_in_file
signal_labels = f.getSignalLabels()
sigbufs = np.zeros((n, f.getNSamples()[0]))
for i in np.arange(n):
	sigbufs[i, :] = f.readSignal(i)

