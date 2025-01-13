import sys
import os

def print_file(idx, file):
	print('%4d %10d %s'%(idx, file[1], file[0]))

if __name__=='__main__':
	debug=0		# <- CAUTION

	if debug:
		path='C:/Projects/datasets/dataset-train/CLIC303'
		count=16
	else:
		if len(sys.argv)!=3:
			print('Usage:  python %s path N'%sys.argv[0])
			print('  Prints N median files by size')
			exit(0)
		path=sys.argv[1]
		count=int(sys.argv[2])
	
	files=[]
	for f in os.listdir(path):
		fn=os.path.join(path, f)
		if os.path.isfile(fn):
			files.append((fn, os.path.getsize(fn)))
	#files = [(f, os.path.getsize(f)) for f in os.listdir(path) if os.path.isfile(f)]

	total=len(files)
	print('Found %d files'%total)
	if total:
		files.sort(key=lambda x: x[1])

		if count>total:
			count=total
		start=(total-count)//2
		end=(total+count)//2

		if start>0:
			print_file(0, files[0])
			if start>1:
				print('...')
		for idx in range(start, end):
			print_file(idx, files[idx])
		if end<total:
			if end<total-1:
				print('...')
			print_file(total-1, files[-1])
