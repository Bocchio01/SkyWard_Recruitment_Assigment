.PHONY: compile run clean

CC = gcc
FLAGS = -Wall -std=c99 -O2
SRCS = *.c
TARGET = OnBoard_Software

all: run

compile: $(SRCS)
	$(CC) $(FLAGS) $(SRCS) -o $(TARGET)
	@echo Compilation complete

run: compile
	@$(TARGET).exe
	@echo Done runnig complete

clean:
	del /Q $(TARGET).exe