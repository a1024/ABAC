import sys
import os
import timeit
import subprocess

if __name__=='__main__':

	debug=0		# <- CAUTION

	if debug:
		codec='f26'
		srcpath='C:/Projects/datasets/dataset-kodak-ppm'
		path1='C:/Projects/datasets/dataset-kodak-temp'
		path2='C:/Projects/datasets/dataset-kodak-temp'
	else:
		if len(sys.argv)!=5:
			print('Usage:  python %s codec src path1 path2'%sys.argv[0])
			exit(0)
		codec=sys.argv[1]
		srcpath=sys.argv[2]
		path1=sys.argv[3]
		path2=sys.argv[4]

	idx=0
	utotal=0
	ctotal=0
	etotal=0
	dtotal=0

	for root, dirs, files in os.walk(srcpath):
		for title in files:
			name, ext=os.path.splitext(title)
			if ext.lower()!='.ppm':
				continue
			srcfn=os.path.join(root, title)
			fn1=os.path.join(path1, name+'.LSIM')
			fn2=os.path.join(path2, name+'.PPM')
			enc_args=[
				codec,
				srcfn,
				fn1,
			]
			dec_args=[
				codec,
				fn1,
				fn2,
			]
			elapsed_enc=timeit.timeit(lambda: subprocess.run(enc_args), number=1)
			elapsed_dec=timeit.timeit(lambda: subprocess.run(dec_args), number=1)
			usize=os.path.getsize(srcfn)
			csize=os.path.getsize(fn1)
			print('%4d %10d->%10d %16f %16f %s'%(idx+1, usize, csize, elapsed_enc, elapsed_dec, name))
			utotal+=usize
			ctotal+=csize
			etotal+=elapsed_enc
			dtotal+=elapsed_dec
			idx+=1
	print('\n%4d %10d->%10d %16f %16f    %16f MB/s %16lf MB/s'%(
		idx,
		utotal, ctotal,
		etotal, dtotal,
		utotal/(etotal*1024*1024), utotal/(dtotal*1024*1024)
	))
