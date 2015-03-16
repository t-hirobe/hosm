# ./exception/Makefile

INCLUDES = -I$(HOME)\include\fftw -I$(HOME)\src\boost_1_51_0 -I$(HOME)\lab\windwave\hosm\hosm
LDFLAGS = -L$(HOME)\include\fftw -L$(HOME)\src\boost_1_51_0\libs

TARGETS = All.o

all : $(TARGETS)

%.o : %.cpp
	$(CXX) -c  $*.cpp $(CFLAGS) $(INCLUDES)

.cpp.o:
	$(CXX) -c $< $(CFLAGS) $(INCLUDES)

clean :
	rm $(TARGETS)
