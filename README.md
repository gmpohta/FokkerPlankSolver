# FokkerPlankSolver

See Readme.pdf

## Install project
First you must to install g++ compiler and Python.
If you are working on windows:
```
$ git clone https://github.com/gmpohta/FokkerPlankSolver.git
$ cd FokkerPlankSolver
$ cd fokker_plank_solver
$ g++ -c main.cpp solver.cpp -fPIC -O
$ g++ -shared -static -o fokker_plank_solver.dll main.o solver.o
$ cd ..
$ cd Python
$ python -m venv ./venv
$ venv/Scripts/activate
$ pip install -r requirements.txt
```

If you are working on linux:
```
$ git clone https://github.com/gmpohta/FokkerPlankSolver.git
$ cd FokkerPlankSolver
$ cd fokker_plank_solver
$ g++ -c main.cpp solver.cpp -fPIC -O
$ g++ -shared -static-libstdc++ -o fokker_plank_solver.so main.o solver.o
$ cd ..
$ cd Python
$ python3 -m venv ./venv
$ source venv/bin/activate
$ pip install -r requirements.txt
```

Run tests:
```
$ cd tests
$ python convergence.py
```
