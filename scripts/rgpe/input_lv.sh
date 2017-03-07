#~/bin/bash
export LD_LIBRARY_PATH=$PWD/build
SYS="lv"
US="_"
EXP="rgpe"
FILE=$SYS$US$EXP$US
echo "Running reg_guess_plots_exp LV Script:"
for i in {2,4,10,20,35,50}; do
    ./bin/reg_guess_plots_exp -s lotka_volterra -t 0:.001:1 -i 0:.001:1 -u .48,.026,.93 -g .5,1.0,.5 -o 1,1,1 -y 35,4 -d $i > results/$EXP/$FILE$i.txt
    (cd results/$EXP/ && python3 graph.py $i $FILE)
done
echo "Finished"
