
TARGET = agc
OBJ = main.o agc.o audio-process.o fixed-point-fft.o vad.o

all : $(TARGET)
$(TARGET) : $(OBJ)
	g++ $^ -o $@

%.o : %.cc
	g++ -O3 -c $< -o $@

clean :
	rm -r *.o
.PHONY : clean

