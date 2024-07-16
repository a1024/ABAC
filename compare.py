import sys
import os
import timeit
import subprocess
from datetime import datetime
import time

if __name__=='__main__':
	start=time.time()
	print(datetime.now().strftime('%Y-%m-%d %H:%M:%S start'))

	debug=1		# <- CAUTION

	if debug:
		codec1='f26'
		codec2='c03'
		srcpath='C:/Projects/datasets/dataset-kodak-ppm'
		tmppath='C:/Projects/datasets/dataset-kodak-temp'
	else:
		if len(sys.argv)!=5:
			print('Usage:  python %s codec1 codec2 srcpath tmppath'%sys.argv[0])
			exit(0)
		codec1=sys.argv[1]
		codec2=sys.argv[2]
		srcpath=sys.argv[3]
		tmppath=sys.argv[4]

	idx=0
	utotal=0
	ctotal1=0
	ctotal2=0
	etotal1=0
	etotal2=0
	dtotal2=0

	for root, dirs, files in os.walk(srcpath):
		for title in files:
			name, ext=os.path.splitext(title)
			if ext.lower()!='.ppm':
				continue
			srcfn=os.path.join(root, title)
			fn1=os.path.join(tmppath, name+'.1.LSIM')
			fn2=os.path.join(tmppath, name+'.2.LSIM')
			fn3=os.path.join(tmppath, name+'.2.PPM')
			enc1_args=[
				codec1,
				srcfn,
				fn1,
			]
			enc2_args=[
				codec2,
				srcfn,
				fn2,
			]
			dec2_args=[
				codec2,
				fn2,
				fn3,
			]
			elapsed_enc1=timeit.timeit(lambda: subprocess.run(enc1_args), number=1)
			elapsed_enc2=timeit.timeit(lambda: subprocess.run(enc2_args), number=1)
			elapsed_dec2=timeit.timeit(lambda: subprocess.run(dec2_args), number=1)
			usize=os.path.getsize(srcfn)
			csize1=os.path.getsize(fn1)
			csize2=os.path.getsize(fn2)
			print('%4d %10d->%10d->%10d (%+10d)  E1 %16f  E2 %16f D2 %16f  %s'%(
				idx+1, usize, csize1, csize2, csize2-csize1, elapsed_enc1, elapsed_enc2, elapsed_dec2, name
			))
			utotal+=usize
			ctotal1+=csize1
			ctotal2+=csize2
			etotal1+=elapsed_enc1
			etotal2+=elapsed_enc2
			dtotal2+=elapsed_dec2
			idx+=1
	print('\n%4d %10d->%10d->%10d (%+10d)  E1 %16f  E2 %16f D2 %16f  %16f %16f %16f MB/s'%(
		idx,
		utotal, ctotal1, ctotal2, ctotal2-ctotal1,
		etotal1, etotal2, etotal2,
		utotal/(etotal1*1024*1024), utotal/(etotal2*1024*1024), utotal/(dtotal2*1024*1024)
	))
	finish=time.time()
	print('%s finish (%f sec)'%(datetime.now().strftime('%Y-%m-%d %H:%M:%S'), finish-start))
