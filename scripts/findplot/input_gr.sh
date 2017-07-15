#~/bin/bash
export LD_LIBRARY_PATH=$PWD/build
SYS="general_repressilator"
US="_"
EXP="findplot"
FILE=$SYS$US$EXP
echo "Running reg_find_plot_exp GR Script:"
./bin/reg_find_plot_exp \
-s $SYS \
-t 0:.01:10 \
-i 0:.01:10 \
-u 216,0.216,2,0.2,0.2,216,0.216,2,0.2,0.2,216,0.216,2,0.2,0.2 \
-o 175,.5,1,.5,.5,175,.5,1,.5,.5,175,.5,1,.5,.5 \
-y 5,1,10,1,15,1 \
> results/$EXP/$FILE.txt
echo "Finished"
