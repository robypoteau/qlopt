#~/bin/bash
echo "Running Toggle Switch Script:"
./bin/reg_exp \
-s toggle_switch \
-t 0:.05:12 -i 0:.05:12 \
-u 1,1,10,2,1,1,1,1,10,2,1,1 \
-g .5,.5,5,.5,.5,.5,.5,.5,5,.5,.5,.5 \
-o .5,.5,5,.5,.5,.5,.5,.5,5,.5,.5,.5 \
-y 5,2,5,2 \
-n 0.00001 \
> results/reg1/ts_rng.txt
echo "Finished Running Toggle Switch Script"
