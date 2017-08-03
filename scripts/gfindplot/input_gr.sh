#~/bin/bash
export LD_LIBRARY_PATH=$PWD/build
SYS="general_repressilator"
US="_"
EXP="g_findplot"
FILE=$SYS$US$EXP
echo "Running g_reg_find_plot_exp GR Script:"
./bin/g_reg_find_plot_exp \
-s $SYS \
-t 0:.01:10 \
-i 0:.01:10 \
-u 216,0.216,2,0.2,0.2,216,0.216,2,0.2,0.2,216,0.216,2,0.2,0.2 \
-o 175,.5,1,.5,.5,175,.5,1,.5,.5,,175,.5,1,.5,.5 \
-g 250,.4,3,.4,.4,250,.4,3,.4,.4,250,.4,3,.4,.4 \
-y 5,1,10,1,15,1 \
> results/$EXP/$FILE.txt
echo "Finished"
