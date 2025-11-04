

# release:	make codec
# profiler:	make codec PROF=1
# debug:	make codec DEBUG=1


ifeq ($(DEBUG), 1)
FLAGS = -Wall -O3 -g -D_DEBUG
else
FLAGS = -Wall -O3 -s
endif
#FLAGS = -Wall -ffunction-sections -fdata-sections -Wl,--gc-sections -O3 -s

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
#	gcc $(FLAGS) -mavx2 -mbmi -mlzcnt $(SRC_ZYDIS) $(SRC) -o codec -static
#
#build-clang:
#	clang $(FLAGS) -mavx2 -mbmi -mlzcnt $(SRC_ZYDIS) $(SRC) -o codec -static

c01:
	gcc $(FLAGS) -DCODEC_EXT=c01 -mavx2 -mbmi -mlzcnt  $(SRC_ZYDIS) c01.c codec.c util.c  -o c01 $(LIBS)

c02:
	gcc $(FLAGS) -DCODEC_EXT=c02 -mavx2 -mbmi -mlzcnt  $(SRC_ZYDIS) c02.c codec.c util.c  -o c02 $(LIBS)

c03:
	gcc $(FLAGS) -DCODEC_EXT=c03 -mavx2 -mbmi -mlzcnt  $(SRC_ZYDIS) c03.c codec.c util.c  -o c03 $(LIBS)

c03-clang:
	clang $(FLAGS) -DCODEC_EXT=c03 -mavx2 -mbmi -mlzcnt  $(SRC_ZYDIS) c03.c codec.c util.c  -o c03-clang $(LIBS)

c04:
	gcc $(FLAGS) -DCODEC_EXT=c04 -mavx2 -mbmi -mlzcnt  $(SRC_ZYDIS) c04.c codec.c util.c  -o c04 $(LIBS)

c05:
	gcc $(FLAGS) -DCODEC_EXT=c05 -mavx2 -mbmi -mlzcnt  $(SRC_ZYDIS) c05.c codec.c util.c  -o c05 $(LIBS)

c07:
	gcc $(FLAGS) -DCODEC_EXT=c07 -mavx2 -mbmi -mlzcnt  $(SRC_ZYDIS) c07.c codec.c util.c  -o c07 $(LIBS)

c08:
	gcc $(FLAGS) -DCODEC_EXT=c08 -mavx2 -mbmi -mlzcnt  $(SRC_ZYDIS) c08.c codec.c util.c  -o c08 $(LIBS)

c09:
	gcc $(FLAGS) -DCODEC_EXT=c09 -mavx2 -mbmi -mlzcnt  $(SRC_ZYDIS) c09.c codec.c util.c  -o c09 $(LIBS)

c10:
	gcc $(FLAGS) -DCODEC_EXT=c10 -mavx2 -mbmi -mbmi2 -mlzcnt  $(SRC_ZYDIS) c10.c codec.c util.c  -o c10 $(LIBS)

c11:
	gcc $(FLAGS) -DCODEC_EXT=c11 -mavx2 -mbmi -mlzcnt  $(SRC_ZYDIS) c11.c codec.c util.c  -o c11 $(LIBS)

c12m:
	gcc $(FLAGS) -DCODEC_EXT=c12 -mavx2 -mbmi -mlzcnt  $(SRC_ZYDIS) c12m.c codec.c util.c  -o c12 $(LIBS)

c12:
	clang $(FLAGS) -DCODEC_EXT=c12 -mavx2 -mbmi -mlzcnt  $(SRC_ZYDIS) c12.c codec.c  -o c12 $(LIBS)

c12debug:
	clang -Wall -g -D_DEBUG -DLOUD -DCODEC_EXT=c12 -mavx2 -mbmi -mlzcnt  $(SRC_ZYDIS) c12.c codec.c  -o c12 $(LIBS)

c12s:
	gcc -Wall -O3 -s -DCODEC_EXT=c12 -mavx2 -mbmi -mlzcnt  c12.c codec.c  -o c12 $(LIBS)

c13:
	gcc $(FLAGS) -DCODEC_EXT=c13 -mavx2 -mbmi -mlzcnt  $(SRC_ZYDIS) c13.c codec.c util.c  -o c13 $(LIBS)

c13-debug:
	gcc -g -Wall -DCODEC_EXT=c13 -mavx2 -mbmi -mlzcnt  $(SRC_ZYDIS) c13.c codec.c util.c  -o c13-debug $(LIBS)

c13-profile:
	gcc -O -g -Wall -DCODEC_EXT=c13 -mavx2 -mbmi -mlzcnt  $(SRC_ZYDIS) c13.c codec.c util.c  -o c13-profile $(LIBS)

#c13-gprofile:
#	gcc -g -Wall -DCODEC_EXT=c13 -mavx2 -mbmi -mlzcnt  $(SRC_ZYDIS) c13.c codec.c util.c  -o c13-gprofile $(LIBS) -Wl,--no-as-needed,-lprofiler,--as-needed

c14:
	gcc $(FLAGS) -DCODEC_EXT=c14 -mavx2 -mbmi -mlzcnt  $(SRC_ZYDIS) c14.c codec.c util.c  -o c14 $(LIBS)

c15:
	gcc $(FLAGS) -DCODEC_EXT=c15 -mavx2 -mbmi2 -mlzcnt  $(SRC_ZYDIS) c15.c codec.c util.c  -o c15 $(LIBS)

c16:
	gcc $(FLAGS) -DCODEC_EXT=c16 -mavx2 -mbmi -mlzcnt  $(SRC_ZYDIS) c16.c codec.c util.c  -o c16 $(LIBS)

c17:
	gcc $(FLAGS) -DCODEC_EXT=c17 -mavx2 -mbmi -mlzcnt  $(SRC_ZYDIS) c17.c codec.c util.c  -o c17 $(LIBS)

c18:
	gcc $(FLAGS) -DCODEC_EXT=c18 -mavx2 -mbmi -mlzcnt  $(SRC_ZYDIS) c18.c codec.c util.c  -o c18 $(LIBS)

c18-clang:
	clang $(FLAGS) -DCODEC_EXT=c18 -mavx2 -mbmi -mlzcnt  $(SRC_ZYDIS) c18.c codec.c util.c  -o c18 $(LIBS)

c19:
	gcc $(FLAGS) -DCODEC_EXT=c19 -mavx2 -mbmi -mlzcnt  $(SRC_ZYDIS) c19.c codec.c util.c  -o c19 $(LIBS)

c19-clang:
	clang $(FLAGS) -DCODEC_EXT=c19 -mavx2 -mbmi -mlzcnt  $(SRC_ZYDIS) c19.c codec.c util.c  -o c19 $(LIBS)

c20:
	gcc $(FLAGS) -Ihtscodecs -DCODEC_EXT=c20 -mavx2 -mbmi -mlzcnt  $(SRC_ZYDIS) c20.c codec.c util.c  -o c20 $(LIBS) \
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
	gcc $(FLAGS) -DCODEC_EXT=c21 -mavx2 -mbmi -mlzcnt  $(SRC_ZYDIS) c21.c codec.c util.c  -o c21 $(LIBS)

c22:
	gcc $(FLAGS) -Ihtscodecs -DCODEC_EXT=c22 -mavx2 -mbmi -mlzcnt  $(SRC_ZYDIS) c22.c codec.c util.c  -o c22 $(LIBS) \
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
	gcc $(FLAGS) -Ihtscodecs -DCODEC_EXT=c23 -mavx2 -mbmi -mlzcnt  $(SRC_ZYDIS) c23.c codec.c util.c  -o c23 $(LIBS) \
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
	clang $(FLAGS) -DCODEC_EXT=c24 -mavx2 -mbmi -mlzcnt  $(SRC_ZYDIS) c24.c codec.c util.c  -o c24 $(LIBS)

c25:
	gcc $(FLAGS) -DCODEC_EXT=c25 -mavx2 -mbmi -mlzcnt  $(SRC_ZYDIS) c25.c codec.c util.c  -o c25 $(LIBS)

c26:
	gcc $(FLAGS) -DCODEC_EXT=c26 -mavx2 -mbmi -mlzcnt  $(SRC_ZYDIS) c26.c codec.c util.c  -o c26 $(LIBS)

c27:
	gcc $(FLAGS) -DCODEC_EXT=c27 -mavx2 -mbmi -mlzcnt  $(SRC_ZYDIS) c27.c codec.c util.c  -o c27 $(LIBS)

c28:
	gcc $(FLAGS) -DCODEC_EXT=c28 -mavx2 -mbmi -mlzcnt  $(SRC_ZYDIS) c28.c codec.c util.c  -o c28 $(LIBS)

c29:
	gcc $(FLAGS) -DCODEC_EXT=c29 -mavx2 -mbmi -mlzcnt  $(SRC_ZYDIS) c29.c codec.c util.c  -o c29 $(LIBS)

c30:
	gcc $(FLAGS) -DCODEC_EXT=c30 -mavx2 -mbmi -mlzcnt  $(SRC_ZYDIS) c30.c codec.c util.c  -o c30 $(LIBS)

c31:
	gcc $(FLAGS) -DCODEC_EXT=c31 -mavx2 -mbmi -mlzcnt  $(SRC_ZYDIS) c31.c codec.c util.c  -o c31 $(LIBS)

c32:
	clang $(FLAGS) -DCODEC_EXT=c32 -mavx2 -mbmi -mlzcnt  $(SRC_ZYDIS) c32.c codec.c util.c  -o c32 $(LIBS)

c33:
	gcc $(FLAGS) -DCODEC_EXT=c33 -mavx2 -mbmi -mlzcnt  $(SRC_ZYDIS) c33.c codec.c util.c  -o c33 $(LIBS)

c34:
	clang $(FLAGS) -DCODEC_EXT=c34 -mavx2 -mbmi -mlzcnt  $(SRC_ZYDIS) c34.c codec.c  -o c34 $(LIBS)

c34-gcc:
	gcc $(FLAGS) -DCODEC_EXT=c34 -mavx2 -mbmi -mlzcnt  $(SRC_ZYDIS) c34.c codec.c  -o c34 $(LIBS)

c36:
	clang $(FLAGS) -DCODEC_EXT=c36 -mavx2 -mbmi -mlzcnt  $(SRC_ZYDIS) c36.c util.c codec.c  -o c36 $(LIBS)

c37:
	clang $(FLAGS) -DCODEC_EXT=c37 -mavx2 -mbmi -mlzcnt  $(SRC_ZYDIS) c37.c codec.c  -o c37 $(LIBS)

c38:
	clang $(FLAGS) -DCODEC_EXT=c38 -mavx2 -mbmi -mlzcnt  $(SRC_ZYDIS) c38.c codec.c  -o c38 $(LIBS)
