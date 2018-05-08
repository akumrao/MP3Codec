CXX        = g++
FLAG       = -Wall -I$(INCLUDE) -O2 -std=c++14
SRC        = src
OBJ        = tmp
BIN        = ./
INCLUDE    = include

EGGPLOT_OBJ = \
	$(OBJ)/linespec.o \
	$(OBJ)/plotPoint.o \
	$(OBJ)/MatLabPlot.o \
	$(OBJ)/main.o \

all: plot

plot: $(EGGPLOT_OBJ)
	$(CXX) -o $(BIN)/$@ $^ 

$(OBJ)/%.o: $(SRC)/%.cpp
	$(CXX) $(FLAG) -c $< -o $@

.phony: clean
clean:
	rm -rf $(OBJ)/*  
