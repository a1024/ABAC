PROG=battle.exe

SRC=$(wildcard src/*.c)
CFLAGS=-march=nehalem
build:	CFLAGS += -g -D_DEBUG
build:	$(PROG)
release:	CFLAGS += -O -D_PROFILE
release:	$(PROG)
$(PROG):	$(SRC)
	$(CC) $(CFLAGS) $(SRC) -o $(PROG)
