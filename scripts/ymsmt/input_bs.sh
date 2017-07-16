#~/bin/bash
export LD_LIBRARY_PATH=$PWD/build
SYS="bistable_switch"
US="_"
EXP="ymsmt"
FILE=$SYS$US$EXP$US
echo "Running reg_guess_plots_exp BS Script:"
for i in {2,10,20,25}; do
    ./bin/reg_ode_soln_msmt_plots_exp \
    -s $SYS \
    -t 0:.01:3 \
    -i 0:.01:3 \
    -u 150,3.2,2 \
    -o 200,4,3 \
    -y 25,25 \
    -d $i \
    > results/$EXP/$FILE$i.txt
    (cd results/$EXP/ && python3 graph.py $i $FILE)
done
echo "Finished"
