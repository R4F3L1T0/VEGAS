#!/bin/bash

g++ -std=c++11 -O3 vegas.cpp
./a.out
python plotter.py