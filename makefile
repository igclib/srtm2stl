debug:
	g++ -std=gnu++20 -g -Iinclude/ src/*.cpp -o bin/srtm2stl

release:
	g++ -std=gnu++20 -O2 -Iinclude/ src/*.cpp -o bin/srtm2stl	

run:
	bin/srtm2stl

