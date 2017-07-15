#~/bin/bash
export LD_LIBRARY_PATH=$PWD/build
SYS="toggle_switch_config1"
US="_"
EXP="ymsmt"
FILE=$SYS$US$EXP$US
echo "Running reg_guess_plots_exp TSC1 Script:"
for i in {20,40,82,124}; do
    ./bin/reg_ode_soln_msmt_plots_exp \
    -s $SYS \
    -t 0:.005:6 \
	-i 0:.005:6 \
	-u 5,3,2,1,5,6,4,2,1,5 \
	-o 2,1,1,1,1,2,1,1,1,2 \
	-y 20,20,20,20 \
    -n 0.005 \
    -d $i \
    > results/$EXP/$FILE$i.txt
    (cd results/$EXP/ && python3 graph.py $i $FILE)
done
echo "Finished"
