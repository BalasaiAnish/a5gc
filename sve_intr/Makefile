PREFIX ?= ./toolchain/arm-gnu-toolchain-14.3.rel1-x86_64-aarch64-none-linux-gnu/bin/aarch64-none-linux-gnu
ARM_CC ?= $(PREFIX)-gcc
ARM_CC_FLAGS += "-march=armv9-a"

INC ?= ./inc
SRC ?= ./src/*.c
TEST ?= ./test/default.c

build:
	$(ARM_CC) $(ARM_CC_FLAGS) -I$(INC) $(SRC) $(TEST) -o test.o

configure:
	bear -- $(ARM_CC) $(ARM_CC_FLAGS) -I$(INC) $(SRC) $(TEST) -o test.o

clean:
	rm *.o
