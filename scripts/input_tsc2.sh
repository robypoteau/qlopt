#~/bin/bash
echo "Running Toggle Switch Second Configuration Script:"
	./bin/main \
	-s toggle_switch_config2 \
	-t 0:.005:6  \
	-i 0:.005:6  \
	-u 5,3,4,1,5,6,4,4,1,5 \
	-o 2,2,2,2,2,2,2,2,2,2 \
	-y 20,20,20,20 \
	-k 1 \
	-p 10 \
	-r 0 \
	-l 0.0005 
echo "Finished Running Toggle Switch Script"