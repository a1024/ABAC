cbm:	cbm.c util.c
	gcc -Wall -O3 -mavx2 -mlzcnt -mbmi util.c cbm.c -o cbm -static -s

cbm2tex:	cbm2tex.c util.c
	gcc -Wall -O3 -mavx2 -mlzcnt -mbmi -static -s util.c cbm2tex.c -o cbm2tex
