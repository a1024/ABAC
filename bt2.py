import sys
import os
import timeit
import subprocess
from datetime import datetime
import time

if __name__=='__main__':
	start=time.time()
	print(datetime.now().strftime('%Y-%m-%d %H:%M:%S'))

	debug=0		# <- CAUTION

	if debug:
		srcpath='C:/Projects/datasets/dataset-kodak-ppm'
		dstpath='C:/Projects/datasets/dataset-kodak-temp'
		srcext='PPM'
		dstext='LSIM'
		args2=['F23', '---s', '---d']
	else:
		if len(sys.argv)<6:
			print('Usage:    python %s srcpath dstpath ext1  ext2  program args...'%sys.argv[0])
			print('')
			print('Example:  python %s srcpath dstpath ppm   jxl   CJXL --quiet --num-threads=-1 -d 0 -e 2 ---s ---d'%sys.argv[0])
			print('Example:  python %s srcpath dstpath png   ppm   FFMPEG -hide_banner -loglevel error -i ---s ---d'%sys.argv[0])
			print('Example:  python %s srcpath dstpath ppm   qlic2 QLIC2 c ---d ---s'%sys.argv[0])
			print('Example:  python %s srcpath dstpath qlic2 ppm   QLIC2 d ---s ---d'%sys.argv[0])
			print('Example:  python %s srcpath dstpath ppm   halic HALIC_ENCODE ---s ---d -mt'%sys.argv[0])
			print('')
			print('args should contain placeholders "---s" for source files, "---d" for destination files')
			exit(0)
		srcpath=sys.argv[1]
		dstpath=sys.argv[2]
		srcext=sys.argv[3]
		dstext=sys.argv[4]
		args2=sys.argv[5:]

	try:
		placeholder1=args2.index('---s')
	except:
		placeholder1=-1
	try:
		placeholder2=args2.index('---d')
	except:
		placeholder2=-1
	if srcext[0]!='.':
		srcext='.'+srcext
	if dstext[0]!='.':
		dstext='.'+dstext

	idx=0
	total1=0
	total2=0
	etotal=0

	for root, dirs, files in os.walk(srcpath):
		for title in files:
			name, ext=os.path.splitext(title)
			if ext.lower()!=srcext.lower():
				continue
			srcfn=os.path.join(root, title)
			dstfn=os.path.join(dstpath, name+dstext)
			enc_args=args2.copy()
			if placeholder1>=0:
				enc_args[placeholder1]=srcfn
			if placeholder2>=0:
				enc_args[placeholder2]=dstfn
			size1=os.path.getsize(srcfn)
			elapsed=timeit.timeit(lambda: subprocess.run(enc_args), number=1)
			size2=os.path.getsize(dstfn)
			if size1<size2:
				usize=size2
			else:
				usize=size1
			print('%4d %11d->%11d (%+11d)  %11.6f sec  %11.6f MB/s  %s'%(
				idx+1, size1, size2, size2-size1, elapsed, usize/(elapsed*1024*1024) if elapsed>0 else float('inf'), name
			))
			total1+=size1
			total2+=size2
			etotal+=elapsed
			idx+=1
			#if idx%50==0:
			#	time.sleep(40)
	if total1<total2:
		utotal=total2
	else:
		utotal=total1
	print('\n%4d %11d->%11d (%+11d)  %11.6f sec  %11.6f MB/s'%(
		idx, total1, total2, total2-total1, etotal, utotal/(etotal*1024*1024) if etotal>0 else float('inf')
	))
	finish=time.time()
	print('%s (%f sec)'%(datetime.now().strftime('%Y-%m-%d %H:%M:%S'), finish-start))
