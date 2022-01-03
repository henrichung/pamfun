import numpy as np
import pandas as pd
from scipy import sparse

print("reading in data")
A = pd.read_csv("binDTM_nosingletons_.csv")
B = A.to_numpy()
Bt = B.transpose()

print("converting to sparse matrix")
sB = sparse.csr_matrix(B)
sBt = sparse.csr_matrix(Bt)

print("calculating dot product")
C = sB.dot(sBt)
print("saving output")
D = pd.DataFrame(data = sparse.csr_matrix.todense(C))
np.savetxt("bin_cooc_sp.csv", D, delimiter=",")


