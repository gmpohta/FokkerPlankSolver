# FokkerPlankSolver

See Readme.pdf

## Install project
First you must to install g++ compiler and python. 
Compile dynamic library:
```
$ git clone https://github.com/gmpohta/FokkerPlankSolver.git
$ cd FokkerPlankSolver
$ cd fokker_plank_solver
```
if you are working on windows:
```
$ g++ -c main.cpp solver.cpp -fPIC -O
$	g++ -shared -static -o fokker_plank_solver.dll main.o solver.o
```
if you are working on linux:
```
$ make compile_so
```
Next create python project:
```
$ cd ..
$ cd Python
```
if you are working on windows:
```
$ python -m venv ./venv
$ venv/Scripts/activate
```
if you are working on linux:
```
$ python3 -m venv ./venv
$ source venv/bin/activate
```
Install python libs:
```
$ pip install -r requirements.txt
```
Run tests:
```
$ cd tests
$ python convergence.py
```
