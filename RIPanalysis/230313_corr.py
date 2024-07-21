import sys
import numpy as np

# Find correlation between two input files

file1 = open(sys.argv[1], "r")
list1 = file1.readlines()
list1_float = [float(x) for x in list1]
file1.close()

file2 = open(sys.argv[2], "r")
list2 = file2.readlines()
list2_float = [float(x) for x in list2]
file2.close()

r_matrix = np.corrcoef(list1_float, list2_float)
r = r_matrix[0,1]

print(r)