# -*- coding: utf-8 -*-
import scipy.io as sio
import numpy as np
import sys

argvs = sys.argv
arglen = len(argvs)
if arglen != 2:
	print 'Usage: python %s [MatrixMarketfile]' % argvs[0]
	quit()
target = argvs[1]
(sizex, sizey, nnz, format, field, symmerty) = sio.mminfo(target)
mat = sio.mmread(target)

ftype = mat.getformat()
min_all = mat.min()
max_all = mat.max()
nnz = mat.getnnz()
nnz_row = mat.getnnz(1)
min_nnz = min(nnz_row)
max_nnz = max(nnz_row)
ave_nnz = np.average(nnz_row)
zero = float(sizex*sizex-nnz)/float(sizex*sizex)*100.0

print "対称性 = %s" % (symmerty)
print "格納形式 = %s" % (ftype)
if symmerty == 'symmetric':
	print "次元数 = %d" % (sizex*sizex)
else :
	print "X-サイズ = %d" % (sizex)
	print "Y-サイズ = %d" % (sizey)
print "非ゼロ要素数 = %d" % (nnz)
print "ゼロ要素割合 = %.3f%%" % (zero)
print "行最大非ゼロ要素数 = %d" % (max_nnz)
print "行最小非ゼロ要素数 = %d" % (min_nnz)
print "行平均非ゼロ要素数 = %d" % (ave_nnz)
print "最大要素数 = %.16e" % (max_all)
print "最小要素数 = %.16e" % (min_all)