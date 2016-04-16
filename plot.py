# -*- coding: utf-8 -*-
import scipy.io as sio
import numpy as np
import sys
import matplotlib.pyplot as plt

mat = sio.mmread("bcsstk18.mtx")
plt.spy(mat, precision=1e-1, marker='.', markersize=1)
plt.show()