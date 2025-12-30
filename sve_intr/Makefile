RM = /bin/rm

PREFIX ?= ./toolchain/arm-gnu-toolchain-14.3.rel1-x86_64-aarch64-none-linux-gnu/bin/aarch64-none-linux-gnu
ARM_CC ?= $(PREFIX)-gcc
ARM_CC_FLAGS += "-Wall"
ARM_CC_FLAGS += "-march=armv8.2-a+crc+sve+nodotprod+norcpc+noprofile"
ARM_CC_FLAGS += "-mcpu=neoverse-n1+sve"
ARM_CC_FLAGS += "-static"

ARM_QEMU ?= "qemu-aarch64"
ARM_QEMU_FLAGS += "-cpu max,sve-default-vector-length=16"

INCS = ./inc
SRCS += ./src/ntt128.c
CASE ?= butterfly
SRCS += ./test/$(CASE).c

.PHONY: build run all

build:
	$(ARM_CC) $(ARM_CC_FLAGS) -I$(INCS) $(SRCS) -o test.o

run:
	$(ARM_QEMU) -cpu max,sve-default-vector-length=16 ./test.o > test.log

all: build run

asm:
	$(ARM_CC) $(ARM_CC_FLAGS) -I$(INCS) $(SRCS) -S

lsp:
	bear -- $(ARM_CC) $(ARM_CC_FLAGS) -I$(INCS) $(SRCS) -o test.o

clean:
	-$(RM) -f *.o *.s
