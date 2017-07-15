#~/bin/bash
export LD_LIBRARY_PATH=$PWD/build
SYS="lv"
US="_"
EXP="rgpe"
FILE=$SYS$US$EXP$US
echo "Running reg_guess_plots_exp LV Script:"
for i in {2,3,4,5,6,7}; do
    ./bin/reg_guess_plots_exp \
    -s lotka_volterra \
    -t 0:.001:1 \
    -i 0:.001:1 \
    -u 1,2,1 \
    -g 5,6,5 \
    -o 10,50,120 \
    -y 1,1 \
    -d $i > results/$EXP/$FILE$i.txt
    (cd results/$EXP/ && python3 graph.py $i $FILE)
done
echo "Finished"
