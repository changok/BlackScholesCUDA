SRC_DIR = ./src
BIN_DIR = ./bin
MAKE_BIN_DIR_TARGET = mkdir_bin
INCLUDE_DIR = ./include
ARCH_OPT = -arch=sm_13

all: parser.o main.o timer.o black.o
	nvcc $(ARCH_OPT) $(BIN_DIR)/main.o $(BIN_DIR)/timer.o $(BIN_DIR)/black_scholes.o $(BIN_DIR)/parser.o -o $(BIN_DIR)/blackScholes 

main.o: $(MAKE_BIN_DIR_TARGET) parser.o
	g++ -I${INCLUDE_DIR} -c $(SRC_DIR)/main.cpp -o $(BIN_DIR)/main.o

timer.o: $(MAKE_BIN_DIR_TARGET)
	g++ -I${INCLUDE_DIR} -c $(SRC_DIR)/timer.cpp -o  $(BIN_DIR)/timer.o
	
black.o: $(MAKE_BIN_DIR_TARGET)
	nvcc $(ARCH_OPT) -I${INCLUDE_DIR} -c $(SRC_DIR)/black_scholes.cu -o  $(BIN_DIR)/black_scholes.o
	
parser.o: $(MAKE_BIN_DIR_TARGET)
	g++ -I${INCLUDE_DIR} -c $(SRC_DIR)/parser.cpp -o  $(BIN_DIR)/parser.o
	
$(MAKE_BIN_DIR_TARGET) :
	mkdir -p $(BIN_DIR)

clean :
	\rm -f $(BIN_DIR)/*
