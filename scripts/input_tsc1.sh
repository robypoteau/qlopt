#~/bin/bash
echo "Running Toggle Switch N Script:"
./bin/main \
-s toggle_switch_config1 \
-t 0:.005:6 \
-i 0:.005:6 \
-u 5,3,2,1,5,6,4,2,1,5 \
-o 2,2,2,2,2,2,2,2,2,2 \
-y 20,20,20,20 \
-k 1 \
-p 10 \
-r 1 \
-l 0.005
echo "Finished Running Toggle Switch N Script"
