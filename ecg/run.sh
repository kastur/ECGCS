#!/bin/sh
sh script/convert_ecg.sh p04 > temp/data.txt
make main 
./main temp/data.txt > temp/out.txt
cd temp
echo 'plot "out.txt" using 2 with lines' > plot_script
gnuplot plot_script

