# Oblivous Rider

## Introduction

The application simulates the situation where an anonymous rider needs to find the closest driver around to pick her up, hinted by the [sample code](https://github.com/ldsec/lattigo/blob/master/examples/bfv/main.go) of Lattigo.

Running this program will give you a report of the execution where you can see the randomly chose coordinates of the rider and some drivers and the final result. It will also report the closest rider, and the corresponding distance.

## How to Use

1. First, you need to get Microsoft SEAL prepared. I recommend to build the library from source and install it [globally]('https://github.com/microsoft/SEAL#installing-microsoft-seal) in your system so that the `CMakeList.txt` file doesn't need to be modified and can be used directly.

2. `cmake -S . -B build` to build out of source.

3. `cd build && make`

4. `./main`

5. You may also tweak the parameters chosen in Step 1 in `main`function to see what's happening.
