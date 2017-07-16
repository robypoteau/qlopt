#~/bin/bash
export LD_LIBRARY_PATH=$PWD/build
SYS="toggle_switch_config1"
US="_"
EXP="g_findplot"
FILE=$SYS$US$EXP
echo "Running g_reg_find_plot_exp TSC1 Script:"
./bin/g_reg_find_plot_exp \
-s $SYS \
-t 0:.005:6 \
-i 0:.005:6 \
-u 5,3,2,1,5,6,4,2,1,5 \
-o 2,1,1,1,1,2,1,1,1,2 \
-y 20,20,20,20 \
> results/$EXP/$FILE.txt
echo "Finished"
