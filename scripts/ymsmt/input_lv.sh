#~/bin/bash
export LD_LIBRARY_PATH=$PWD/build
SYS="lotka_volterra"
US="_"
EXP="ymsmt"
FILE=$SYS$US$EXP$US
echo "Running reg_guess_plots_exp LV Script:"
for i in {2,3,4,5,7}; do
    ./bin/reg_ode_soln_msmt_plots_exp \
    -s $SYS \
    -t 0:.1:1 \
    -i 0:.1:1 \
    -u 1,2,1 \
    -o 5,6,5 \
    -y 1,1 \
    -d $i \
    > results/$EXP/$FILE$i.txt
    (cd results/$EXP/ && python3 graph.py $i $FILE)
done
echo "Finished"
