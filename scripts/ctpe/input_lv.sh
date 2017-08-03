#~/bin/bash
export LD_LIBRARY_PATH=$PWD/build
SYS="lv"
US="_"
EXP="ctpe"
FILE=$SYS$US$EXP$US
echo "Running reg_guess_plots_exp LV Script:"
for i in {2,3,4,5}; do
    ./bin/curvature_test_plots_exp \
    -s lotka_volterra \
    -t 0:.001:1 \
    -i 0:.001:1 \
    -u 1,2,1 \
    -o 5,6,5 \
    -y 1,1 \
    -d $i \
    > results/$EXP/$FILE$i.txt
    (cd results/$EXP/ && python3 graph.py $i $FILE)
done
echo "Finished"
