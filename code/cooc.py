import numpy as np
import pandas as pd
A = pd.read_csv("d3.csv")
B = A.to_numpy()
Bt = B.transpose()
C = np.dot(B, Bt)
np.savetxt("cooc.csv", C, delimiter=",")

#library(tidyverse)
#A = readr::read_csv("d3.csv") 
#B = A %>% replace(is.na(.), 0)
#readr::write_csv(B, "d3.csv")