cmake -DDEBUG=OFF -DCMAKE_CXX_COMPILER:STRING="C:/MinGW/bin/g++.exe" -G "MinGW Makefiles" -S . -B build
cmake --build ./build
./build/src/MO_Labs.exe