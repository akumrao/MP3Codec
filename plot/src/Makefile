CXX        = g++
FLAG       = -Wall -I$(INCLUDE) -g -std=c++11
SRC        = src
OBJ        = tmp
BIN        = ./
INCLUDE    = ./

PLOT_OBJ = \
	$(OBJ)/linespec.o \
	$(OBJ)/plot.o \
	$(OBJ)/main.o \

all: plot

plot: $(PLOT_OBJ)
	$(CXX) -o $(BIN)/$@ $^ 

$(OBJ)/%.o: $(SRC)/%.cpp
	$(CXX) $(FLAG) -c $< -o $@

.phony: clean
clean:
	rm -rf $(OBJ)/*  
