

FSESRC =\
 fse/debug.c\
 fse/entropy_common.c\
 fse/fse_compress.c\
 fse/fse_decompress.c\
 fse/fseU16.c\
 fse/hist.c\
 fse/huf_compress.c\
 fse/huf_decompress.c


ifeq ($(OS), Windows_NT)
LIBS += -static
else
LIBS += -lm
endif

a01:
	gcc -Wall -O3 -s -mlzcnt -DCODEC_EXT=a01 a01.c sli.c -o a01 $(LIBS)

a02:
	gcc -Wall -O3 -s -mlzcnt -DCODEC_EXT=a02 a02.c sli.c -o a02 $(LIBS)

a03:
	gcc -Wall -O3 -s -mlzcnt -DCODEC_EXT=a03 a03.c sli.c -o a03 $(LIBS)

a04:
	gcc -Wall -O3 -s -mlzcnt -DCODEC_EXT=a04 a04.c fpc/fpc.c sli.c -o a04 $(LIBS)

a04avx2:
	gcc -Wall -O3 -s -mlzcnt -mavx2 -DCODEC_EXT=a04 a04.c fpc/fpc.c sli.c -o a04avx2 $(LIBS)

a05:
	gcc -Wall -O3 -s -mlzcnt -DCODEC_EXT=a05 a05.c fpc/fpc.c $(FSESRC) sli.c -o a05 $(LIBS)

a05avx2:
	gcc -Wall -O3 -s -mlzcnt -mavx2 -DCODEC_EXT=a05 a05.c fpc/fpc.c $(FSESRC) sli.c -o a05avx2 $(LIBS)
