CFLAGS = -IC:/Projects/a/OpenCL-Headers
LIBS = C:\Windows\System32\OpenCL.dll
ifeq ($(OS),Windows_NT)
LIBS += -static
else
LIBS += -lm
endif

r01:
	gcc -O3 -Wall -Ihtscodecs -DCODEC_FUNC=r01_codec $(CFLAGS) -mavx2 -mbmi -mlzcnt  r01.c rapid.c util.c  -o r01 $(LIBS) \
	htscodecs/arith_dynamic.c \
	htscodecs/fqzcomp_qual.c \
	htscodecs/htscodecs.c \
	htscodecs/pack.c \
	htscodecs/rANS_static.c \
	htscodecs/rANS_static4x16pr.c \
	htscodecs/rANS_static32x16pr.c \
	htscodecs/rANS_static32x16pr_avx2.c \
	htscodecs/rANS_static32x16pr_sse4.c \
	htscodecs/rle.c \
	htscodecs/tokenise_name3.c \
	htscodecs/utils.c
