out: main.cpp
	g++ -o out main.cpp -pthread -lbcm2835 -lfftw3

clean:
	rm *.o out -rf
