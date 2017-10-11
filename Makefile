CXXFLAGS := -O3 -Wall 
INCLUDES := -Iboost -I/usr/local/include/openbabel-2.0
LIBS=-L/usr/local/lib -lm -lopenbabel

all: fragment analyze

fragment: fragment.C
	g++ $(CXXFLAGS) -o fragment fragment.C utilities.C $(INCLUDES) $(LIBS)

analyze: analyze.C
	g++ $(CXXFLAGS) -o analyze analyze.C utilities.C $(INCLUDES) $(LIBS)

clean:
	rm fragment analyze
