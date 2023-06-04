compile_dll:
	g++ -c main.cpp solver.cpp -fPIC -O
	g++ -shared -static -o fokker_plank_solver.dll main.o solver.o

compile_so:
	g++ -c main.cpp solver.cpp -fPIC -O
	g++ -shared -static-libstdc++ -o fokker_plank_solver.so main.o solver.o
