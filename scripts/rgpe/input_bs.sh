#~/bin/bash
export LD_LIBRARY_PATH=$PWD/build
SYS="bistable_switch"
US="_"
EXP="rgpe"
FILE=$SYS$US$EXP$US
echo "Running reg_guess_plots_exp BS Script:"
for i in {2,10,20,25}; do
    ./bin/reg_guess_plots_exp \
    -s $SYS \
    -t 0:.01:3 \
    -i 0:.01:3 \
    -u 150,3.2,2 \
    -g 135,2.0,1.5 \
    -o 200,4,3 \
    -y 25,25 \
    -n 0.01 \
    -d $i > results/$EXP/$FILE$i.txt
    (cd results/$EXP/ && python3 graph.py $i $FILE)
done
echo "Finished"
