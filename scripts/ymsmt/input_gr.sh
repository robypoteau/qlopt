#~/bin/bash
export LD_LIBRARY_PATH=$PWD/build
SYS="general_repressilator"
US="_"
EXP="ymsmt"
FILE=$SYS$US$EXP$US
echo "Running reg_guess_plots_exp GR Script:"
for i in {3,6,9}; do
    ./bin/reg_ode_soln_msmt_plots_exp \
    -s $SYS \
    -t 0:.01:10 \
	-i 0:.01:10 \
	-u 216,0.216,2,0.2,0.2,216,0.216,2,0.2,0.2,216,0.216,2,0.2,0.2 \
	-o 150,.5,.5,.5,.5,150,.5,.5,.5,.5,,150,.5,.5,.5,.5 \
	-y 5,1,10,1,15,1 \
    -d $i \
    > results/$EXP/$FILE$i.txt
    (cd results/$EXP/ && python3 graph.py $i $FILE)
done
echo "Finished"
