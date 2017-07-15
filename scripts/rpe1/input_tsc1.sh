#~/bin/bash
export LD_LIBRARY_PATH=$PWD/build
echo "Running RPE1 TSC1 Script:"
for i in {20,40,82,124}; do
    ./bin/reg_param_exp1 \
    -s toggle_switch_config1 \
    -t 0:.005:6 \
    -i 0:.005:6 \
    -u 5,3,2,1,5,6,4,2,1,5 \
    -o 2,1,1,1,1,2,1,1,1,2 \
    -y 20,20,20,20 \
    -n 0.005 \
    -d $i \
    > results/rpe1/tsc1$i.txt
    (cd results/rpe1/ && python3 graph.py $i "tsc1")
done
echo "Finished"
