#~/bin/bash
export LD_LIBRARY_PATH=$PWD/build
echo "Running RPE1 Eight Script:"
for i in {100,6,11,16,22,26,33,26,44}; do
    ./bin/reg_param_exp1 -s toggle_switch -t 0:.005:6 -i 0:.005:6 -u 5,3,2,1,5,6,4,2,1,5 -o 3,2,1,2,2,3,2,1,2,2 -y 20,20,20,20 -d $i > results/rpe1/ts$i.txt
    #(cd results/rpe1/ && python3 graph.py $i "ts")
done
echo "Finished"
