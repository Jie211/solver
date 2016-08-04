# -*- coding: utf-8 -*-
import scipy.io as sio
import numpy as np
import sys
import os

argvs = sys.argv
arglen = len(argvs)
if arglen != 2:
	print 'Usage: python %s [MatrixMarketfile]' % argvs[0]
	quit()
target = argvs[1]
path = './'+target.replace('.mtx', '')
colvalpath = path+'/ColVal.txt'
ptrpath = path+'/Ptr.txt'
bxpath = path+'/bx.txt'

matcoo = sio.mmread(target)
print '読み込み　完了'
matscr = matcoo.tocsr()
print 'SCR変換　完了'
size = matscr.shape[0]
nnz = matscr.nnz
data = matscr.data
col = matscr.indices
ptr = matscr.indptr
print 'パラメータ　完了'
xvec = []
for i in xrange(size):
	xvec.append(1)
bvec = matscr.dot(xvec)
print '未知ベクトルを1になるように右辺ベクトルを生成　完了'
os.mkdir(path)
f_colval = open(colvalpath, 'w')
f_ptr = open(ptrpath, 'w')
f_bx = open(bxpath, 'w')
print 'ファイルオープン　完了'
head = "%d %d %d\n" % (size, size, nnz)
f_colval.write(head)
f_ptr.write(head)
f_bx.write(head)
for i in xrange(nnz):
	tmp = "%d %.16e\n" % (col[i], data[i])
	f_colval.write(tmp)
for i in ptr:
	tmp = "%d\n" % (i)
	f_ptr.write(tmp)
for i in xrange(size):
	tmp = "%.16e %.16e\n" % (bvec[i], xvec[i])
	f_bx.write(tmp)
print 'ファイル書き込み　完了'
f_colval.close()
f_ptr.close()
f_bx.close()
print 'SCR変換　完了'