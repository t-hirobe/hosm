.SUFFIXES: .cpp .o
#CFLAGS = -O3
CFLAGS = -O3 -pg -DNDEBUG -Wall
INCLUDES = -I$(HOME)\include\fftw -I$(HOME)\src\boost_1_51_0 -I$(HOME)\lab\windwave\hosm\hosm
LDFLAGS = -L$(HOME)\include\fftw -L$(HOME)\src\boost_1_51_0\libs
CXX = g++
OBJS = WCreater.o HosmCalculator.o GradientCalculator.o SpectrumProperty.o Logger.o SpectrumCreater.o util.o Property.o WaveBinData.o StatisticLogger.o FFT2D.o Exception.o
export

all :
	@echo 'make <prog name>'

hosm: hosm.o $(OBJS)
	$(CXX) -pg -o $@ $(OBJS) $< $(INCLUDE) $(LDFLAGS) -lfftw3-3

subsystem :
	$(MAKE) -C exception


.cpp.o:
	$(CXX) -c $< $(CFLAGS) $(INCLUDES)

clean:
	rm *.o
