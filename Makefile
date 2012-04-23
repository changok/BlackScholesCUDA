SRC_DIR = ./src
BIN_DIR = ./bin
MAKE_BIN_DIR_TARGET = mkdir_bin
INCLUDE_DIR = ./include

all: main.o timer.o black.o

main.o: $(MAKE_BIN_DIR_TARGET)
	g++ -I${INCLUDE_DIR} -c $(SRC_DIR)/main.cpp -o $(BIN_DIR)/main.o

timer.o: $(MAKE_BIN_DIR_TARGET)
	g++ -I${INCLUDE_DIR} -c $(SRC_DIR)/timer.cpp -o  $(BIN_DIR)/timer.o
	
black.o: $(MAKE_BIN_DIR_TARGET)
	nvcc -I${INCLUDE_DIR} -c $(SRC_DIR)/black_scholes.cu -o  $(BIN_DIR)/black_scholes.o
	
$(MAKE_BIN_DIR_TARGET) :
	mkdir -p $(BIN_DIR)

clean :
	\rm -f $(BIN_DIR)/*
