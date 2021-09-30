#!/usr/bin/env python3
import jutil as ju
import numpy as np

a=np.array([-1, 0, 0.1, 0.9, 1, 2],dtype=np.float)
e=ju.mix_entropy(a)
print(a)
print(e)
e=ju.mix_entropy(0.1)
print(a)
print(e)
