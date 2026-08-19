cbm:	cbm.c util.c
	gcc -Wall -O3 -mavx2 -mlzcnt -mbmi -static -s util.c cbm.c -o cbm -lm

cbm2tex:	cbm2tex.c util.c
	gcc -Wall -O3 -mavx2 -mlzcnt -mbmi -static -s util.c cbm2tex.c -o cbm2tex -lm

cb2:	cb2.c
	gcc -Wall -O3 -mavx2 -mlzcnt -mbmi -static -s cb2.c -o cb2 -lm -lole32 -lcomdlg32 -luuid
