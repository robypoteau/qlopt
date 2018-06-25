#~/bin/bash
echo "Running Toggle Switch Script:"
./bin/main \
-s toggle_switch \
-m \
-t 0:.1:24 \
-i 0:.1:24 \
-u 1,1,10,2,1,1,1,1,10,2,1,1 \
-o .5,.5,5,1,.5,.5,.5,.5,5,1,.5,.5 \
-y 5,2,5,2 \
-k 1 \
-p 4 \
-r 1 \
-l 0.0001
echo "Finished Running Toggle Switch Script"
