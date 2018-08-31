build: molopt-sf.cpp
	g++ -O3 -o molopt.exe molopt-sf.cpp

PHONY:.clean
clean:
	rm -r -f coords alphas *.txt *.csv molopt.exe
