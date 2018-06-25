#~/bin/bash
echo "Running Bistable Script:"
./bin/main \
-s bistable_switch \
-m \
-t 0:.01:3 \
-i 0:.1:3 \
-u 150,3.2,2.0 \
-o 125,2.2,1.0 \
-y 14,13 \
-k 3 \
-d 1 \
-p 5 \
-n 1 \
-r 0
echo "Finished Running BS Script"
