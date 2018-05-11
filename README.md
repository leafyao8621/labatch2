compile w/
gcc -o test src/backend/*.c -lm

When run, pass in data file name, parameter file name and output file name

e.g.
./test data/lab2in.dat lab2par.dat out/lab2outc.dat

compile gui version w/

gcc -o lab2g src/backend/Lab2data.c src/frontend/*.c `pkg-config --cflags gtk+-3.0` `pkg-config --libs gtk+-3.0` -lm

run w/

./lab2g

