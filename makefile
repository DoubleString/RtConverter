CC= g++
### Directories
DEBUG = ./Debug/
Json = ./src/Json/
Rnxobs = ./src/Rnxobs/
RtConverter = ./src/RtConverter/
Rtklib =./src/Rtklib/
VrsObs =./src/VrsObs/
Rnxbrd = ./src/Rnxbrd/
Encoder = ./src/Encoder/

### Patterns 
C_SRC= $(wildcard *.cpp $(Json)/*.cpp  $(Encoder)/*.cpp $(Rnxobs)/*.cpp $(RtConverter)/*.cpp $(Rtklib)/*.cpp $(VrsObs)/*.cpp $(Rnxbrd)/*.cpp)
C_OBJ= $(patsubst %.cpp,%.o,$(C_SRC))
### TMP OBJECT
TMP_OBJ= $(patsubst %.cpp,%.o,$(notdir $(C_SRC)))
F_OBJ= $(patsubst %.o,$(DEBUG)%.o,$(TMP_OBJ))

### Target
TARGET=RtConverter
.PHONY: all clean
all:$(TARGET)
##sudo apt-get install  libcurl4-openssl-dev
### Compiler Flags
LCFLAGS= -g -c -std=c++11 
#######################################################################################
$(TARGET):$(C_OBJ)
	$(CC) -o $(TARGET) $(TMP_OBJ)  -lpthread -lm
	mv $(TMP_OBJ) $(DEBUG) 
%.o:%.cpp
	$(CC) $(LCFLAGS) $<
clean:
	rm $(F_OBJ) $(TARGET) $(TMP_OBJ)
