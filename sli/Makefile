
ifeq ($(OS), Windows_NT)
LIBS += -static
else
LIBS += -lm
endif

a01:
	gcc -Wall -O3 -s -mlzcnt -DCODEC_EXT=a01 a01.c sli.c -o a01 $(LIBS)

a02:
	gcc -Wall -O3 -s -mlzcnt -DCODEC_EXT=a02 a02.c sli.c -o a02 $(LIBS)
