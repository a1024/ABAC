
ifeq ($(OS), Windows_NT)
LIBS += -static
else
LIBS += -lm
endif

cb2:	cb2.c
	gcc -Wall -O3 -mavx2 -mlzcnt -mbmi cb2.c -o cb2 $(LIBS) -s -lole32 -lcomdlg32 -luuid

cb2tex:	cb2.c
	gcc -Wall -O3 -mavx2 -mlzcnt -mbmi -DCB2TEX cb2.c -o cb2tex $(LIBS) -s
