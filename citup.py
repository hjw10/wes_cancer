import sys
import h5py
hf=h5py.File(sys.argv[1],'r')
hf.keys()
opnum=hf["results/optimal/index"][0]
cellfreq=hf["trees/" + str(opnum) + "/clone_freq/block0_values"][:]
tree=hf["trees/" + str(opnum) + "/adjacency_list/block0_values"][:]
print cellfreq
print tree