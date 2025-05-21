s01:
	gcc -DCODECFUNC=s01_codec -Wall -O3 -s -static main.c s01.c -o s01

s02:
	gcc -DCODECFUNC=s02_codec -Wall -O3 -s -static main.c s02.c -o s02

s02avx2:
	gcc -DCODECFUNC=s02_codec -Wall -O3 -s -static -mavx2 main.c s02.c -o s02

s03:
	gcc -DCODECFUNC=s03_codec -Wall -O3 -s -static -mavx2 main.c s03.c -o s03

s04:
	gcc -DCODECFUNC=s04_codec -Wall -O3 -s -static -mavx2 main.c s04.c -o s04

s05:
	gcc -DCODECFUNC=s05_codec -Wall -O3 -s -static -mavx2 main.c s05.c -o s05

s06:
	gcc -DCODECFUNC=s06_codec -DM_LZCNT -Wall -O3 -s -static -mavx2 -mlzcnt main.c s06.c -o s06
