SRC_DIR = ./src
BIN_DIR = ./bin
MAKE_BIN_DIR_TARGET = mkdir_bin
INCLUDE_DIR = ./include
ARCH_OPT = -arch=sm_13
FLAG = -D__GOGO_DEBUG__
CFLAGS += $(FLAG)
HW1_EXE = blackScholes
PARAM = params.txt

all: parser.o main.o timer.o black.o bsconfig.o
	nvcc $(ARCH_OPT) $(CFLAGS) $(BIN_DIR)/main.o $(BIN_DIR)/timer.o $(BIN_DIR)/black_scholes.o $(BIN_DIR)/parser.o $(BIN_DIR)/bsconfig.o -o $(BIN_DIR)/blackScholes 

main.o: $(MAKE_BIN_DIR_TARGET) parser.o bsconfig.o
	nvcc $(ARCH_OPT) $(CFLAGS) -I${INCLUDE_DIR} -c $(SRC_DIR)/main.cpp -o $(BIN_DIR)/main.o

timer.o: $(MAKE_BIN_DIR_TARGET)
	g++ -I${INCLUDE_DIR} $(CFLAGS) -c $(SRC_DIR)/timer.cpp -o  $(BIN_DIR)/timer.o
	
black.o: $(MAKE_BIN_DIR_TARGET) bsconfig.o
	nvcc $(ARCH_OPT) $(CFLAGS) -I${INCLUDE_DIR} -c $(SRC_DIR)/black_scholes.cu -o  $(BIN_DIR)/black_scholes.o
	
parser.o: $(MAKE_BIN_DIR_TARGET)
	nvcc $(ARCH_OPT) $(CFLAGS) -I${INCLUDE_DIR} -c $(SRC_DIR)/parser.cpp -o  $(BIN_DIR)/parser.o
	
bsconfig.o: $(MAKE_BIN_DIR_TARGET)
	g++ -I${INCLUDE_DIR} $(CFLAGS) -c $(SRC_DIR)/bsconfig.cpp -o  $(BIN_DIR)/bsconfig.o
	
$(MAKE_BIN_DIR_TARGET) :
	mkdir -p $(BIN_DIR)

clean:
	\rm -f $(BIN_DIR)/*

cuda_4096_run:
	number=1 ; while [[ $$number -le 10 ]] ; \
	do \
		$(BIN_DIR)/$(HW1_EXE) $(PARAM) 4096 0; \
		((number = number + 1)) ; \
	done
cuda_65536_run:
	number=1 ; while [[ $$number -le 10 ]] ; \
	do \
		$(BIN_DIR)/$(HW1_EXE) $(PARAM) 65536 0; \
		((number = number + 1)) ; \
	done
cuda_17367040_run:
	number=1 ; while [[ $$number -le 10 ]] ; \
	do \
		$(BIN_DIR)/$(HW1_EXE) $(PARAM) 17367040 0; \
		((number = number + 1)) ; \
	done

corr_fixed:
	for number in 256 4096 65536 17367040; do \
		$(BIN_DIR)/$(HW1_EXE) $(PARAM) $$number 2; \
	done

corr_only1:
	for number in 256 4096 65536 17367040; do \
		$(BIN_DIR)/$(HW1_EXE) $(PARAM) $$number 1; \
	done
