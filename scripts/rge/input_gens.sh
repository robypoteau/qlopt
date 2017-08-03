#~/bin/bash
echo "Running Toggle Switch Script:"
./bin/reg_guess_exp -s toggle_switch \
-t 0:.1:12 -i 0:.1:12 \
-u 1,1,10,2,1,1,1,1,10,2,1,1 \
-g 2,2,15,4,2,2,2,2,15,4,2,2 \
-o .5,.5,5,.5,.5,.5,.5,.5,5,.5,.5,.5 \
-y 5,2,5,2 \
-n 0.00005 \
> results/rge/ts_rng.txt
echo "Finished Running Toggle Switch Script"
#-u 1,1,10,2,1,1,1,1,10,2,1,1 \
#-g 2,2,15,4,2,2,2,2,15,4,2,2 \
#-o .5,.5,5,.5,.5,.5,.5,.5,5,.5,.5,.5 \
