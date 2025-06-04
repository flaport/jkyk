#!/bin/bash
gcc -lm -O3 -ffast-math  "./$1.c" -o "./$1".out && time "./$1".out
