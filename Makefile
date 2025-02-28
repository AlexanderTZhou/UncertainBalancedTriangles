CC = g++
CFLAGS = -g -Wall -std=c++11
TARGET = baltri

all: $(TARGET)

$(TARGET) : $(TARGET).cpp
	$(CC) $(CFLAGS) -o $(TARGET) $(TARGET).cpp
