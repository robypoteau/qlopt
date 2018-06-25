#~/bin/bash
echo "Running Toggle Switch Script:"
./bin/main \
-s toggle_switch_config4 \
-t 0:.1:6 \
-i 0:.1:6 \
-u 1,2,10,1,1,1,2,10,1,1 \
-o .5,1,5,.5,.5,.5,1,5,.5,.5 \
-y 5,2,5,2 \
-k 5 \
-p 9 \
-r 0 \
-n 1 \
-l 0.00001
echo "Finished Running Toggle Switch Script"
