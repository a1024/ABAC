

# build:		make codec
# with profiler:	make codec PROF=1


ifeq ($(PROF), 1)
SRC_ZYDIS = -DPROFILER -DDISABLE_MT -DENABLE_PROFILER_DISASSEMBLY -Wno-unused-const-variable -Izydist_amalgamated410 zydist_amalgamated410/Zydis.c
else
SRC_ZYDIS =
endif

ifeq ($(OS), Windows_NT)
LIBS += -static
else
LIBS += -lm
endif

#SRC = c03.c codec.c util.c
#SRC = c07.c codec.c util.c
#SRC = c05.c codec.c util.c
#SRC = c01.c c02.c c03.c c04.c c05.c c06.c c07.c codec.c util.c

#.PHONY: c01 c02 c03 c03-clang c05 c07 c09 c10 c11 c12 c13 c13-debug c13-profile

#build:
#	gcc -O3 -Wall -mavx2 -mbmi -mlzcnt $(SRC_ZYDIS) $(SRC) -o codec -static
#
#build-clang:
#	clang -O3 -Wall -mavx2 -mbmi -mlzcnt $(SRC_ZYDIS) $(SRC) -o codec -static

c01:
	gcc -O3 -Wall -DCODEC_FUNC=c01_codec -mavx2 -mbmi -mlzcnt  $(SRC_ZYDIS) c01.c codec.c util.c  -o c01 $(LIBS)

c02:
	gcc -O3 -Wall -DCODEC_FUNC=c02_codec -mavx2 -mbmi -mlzcnt  $(SRC_ZYDIS) c02.c codec.c util.c  -o c02 $(LIBS)

c03:
	gcc -O3 -Wall -DCODEC_FUNC=c03_codec -mavx2 -mbmi -mlzcnt  $(SRC_ZYDIS) c03.c codec.c util.c  -o c03 $(LIBS)

c03-clang:
	clang -O3 -Wall -DCODEC_FUNC=c03_codec -mavx2 -mbmi -mlzcnt  $(SRC_ZYDIS) c03.c codec.c util.c  -o c03-clang $(LIBS)

c04:
	gcc -O3 -Wall -DCODEC_FUNC=c04_codec -mavx2 -mbmi -mlzcnt  $(SRC_ZYDIS) c04.c codec.c util.c  -o c04 $(LIBS)

c05:
	gcc -O3 -Wall -DCODEC_FUNC=c05_codec -mavx2 -mbmi -mlzcnt  $(SRC_ZYDIS) c05.c codec.c util.c  -o c05 $(LIBS)

c07:
	gcc -O3 -Wall -DCODEC_FUNC=c07_codec -mavx2 -mbmi -mlzcnt  $(SRC_ZYDIS) c07.c codec.c util.c  -o c07 $(LIBS)

c08:
	gcc -O3 -Wall -DCODEC_FUNC=c08_codec -mavx2 -mbmi -mlzcnt  $(SRC_ZYDIS) c08.c codec.c util.c  -o c08 $(LIBS)

c09:
	gcc -O3 -Wall -DCODEC_FUNC=c09_codec -mavx2 -mbmi -mlzcnt  $(SRC_ZYDIS) c09.c codec.c util.c  -o c09 $(LIBS)

c10:
	gcc -O3 -Wall -DCODEC_FUNC=c10_codec -mavx2 -mbmi -mlzcnt  $(SRC_ZYDIS) c10.c codec.c util.c  -o c10 $(LIBS)

c11:
	gcc -O3 -Wall -DCODEC_FUNC=c11_codec -mavx2 -mbmi -mlzcnt  $(SRC_ZYDIS) c11.c codec.c util.c  -o c11 $(LIBS)

c12:
	gcc -O3 -Wall -DCODEC_FUNC=c12_codec -mavx2 -mbmi -mlzcnt  $(SRC_ZYDIS) c12.c codec.c util.c  -o c12 $(LIBS)

c13:
	gcc -O3 -Wall -DCODEC_FUNC=c13_codec -mavx2 -mbmi -mlzcnt  $(SRC_ZYDIS) c13.c codec.c util.c  -o c13 $(LIBS)

c13-debug:
	gcc -g -Wall -DCODEC_FUNC=c13_codec -mavx2 -mbmi -mlzcnt  $(SRC_ZYDIS) c13.c codec.c util.c  -o c13-debug $(LIBS)

c13-profile:
	gcc -O -g -Wall -DCODEC_FUNC=c13_codec -mavx2 -mbmi -mlzcnt  $(SRC_ZYDIS) c13.c codec.c util.c  -o c13-profile $(LIBS)

#c13-gprofile:
#	gcc -g -Wall -DCODEC_FUNC=c13_codec -mavx2 -mbmi -mlzcnt  $(SRC_ZYDIS) c13.c codec.c util.c  -o c13-gprofile $(LIBS) -Wl,--no-as-needed,-lprofiler,--as-needed

c14:
	gcc -O3 -Wall -DCODEC_FUNC=c14_codec -mavx2 -mbmi -mlzcnt  $(SRC_ZYDIS) c14.c codec.c util.c  -o c14 $(LIBS)

c15:
	gcc -O3 -Wall -DCODEC_FUNC=c15_codec -mavx2 -mbmi2 -mlzcnt  $(SRC_ZYDIS) c15.c codec.c util.c  -o c15 $(LIBS)

c16:
	gcc -O3 -Wall -DCODEC_FUNC=c16_codec -mavx2 -mbmi -mlzcnt  $(SRC_ZYDIS) c16.c codec.c util.c  -o c16 $(LIBS)

c17:
	gcc -O3 -Wall -DCODEC_FUNC=c17_codec -mavx2 -mbmi -mlzcnt  $(SRC_ZYDIS) c17.c codec.c util.c  -o c17 $(LIBS)

c18:
	gcc -O3 -Wall -DCODEC_FUNC=c18_codec -mavx2 -mbmi -mlzcnt  $(SRC_ZYDIS) c18.c codec.c util.c  -o c18 $(LIBS)

c18-clang:
	clang -O3 -Wall -DCODEC_FUNC=c18_codec -mavx2 -mbmi -mlzcnt  $(SRC_ZYDIS) c18.c codec.c util.c  -o c18 $(LIBS)

c19:
	gcc -O3 -Wall -DCODEC_FUNC=c19_codec -mavx2 -mbmi -mlzcnt  $(SRC_ZYDIS) c19.c codec.c util.c  -o c19 $(LIBS)

c19-clang:
	clang -O3 -Wall -DCODEC_FUNC=c19_codec -mavx2 -mbmi -mlzcnt  $(SRC_ZYDIS) c19.c codec.c util.c  -o c19 $(LIBS)

c20:
	gcc -O3 -Wall -Ihtscodecs -DCODEC_FUNC=c20_codec -mavx2 -mbmi -mlzcnt  $(SRC_ZYDIS) c20.c codec.c util.c  -o c20 $(LIBS) \
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

c21:
	gcc -O3 -Wall -DCODEC_FUNC=c21_codec -mavx2 -mbmi -mlzcnt  $(SRC_ZYDIS) c21.c codec.c util.c  -o c21 $(LIBS)

c22:
	gcc -O3 -Wall -Ihtscodecs -DCODEC_FUNC=c22_codec -mavx2 -mbmi -mlzcnt  $(SRC_ZYDIS) c22.c codec.c util.c  -o c22 $(LIBS) \
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

c23:
	gcc -O3 -Wall -Ihtscodecs -DCODEC_FUNC=c23_codec -mavx2 -mbmi -mlzcnt  $(SRC_ZYDIS) c23.c codec.c util.c  -o c23 $(LIBS) \
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

c24:
	gcc -O3 -Wall -DCODEC_FUNC=c24_codec -mavx2 -mbmi -mlzcnt  $(SRC_ZYDIS) c24.c codec.c util.c  -o c24 $(LIBS)

c25:
	gcc -O3 -Wall -DCODEC_FUNC=c25_codec -mavx2 -mbmi -mlzcnt  $(SRC_ZYDIS) c25.c codec.c util.c  -o c25 $(LIBS)
