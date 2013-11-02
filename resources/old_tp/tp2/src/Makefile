OBJECTOS   = main.o
EJECUTABLE = main
CXXFLAGS   = -g -Wall

all: $(OBJECTOS)
	g++ $(CXXFLAGS) -o $(EJECUTABLE) $(OBJECTOS)

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c -o $@ $<

clean:
	@rm -Rf $(EJECUTABLE).exe $(EJECUTABLE) $(OBJECTOS)

new: clean all
