#~/bin/bash
echo "Running Toggle Switch N Script:"
	./bin/reg_exp \
	-s toggle_switch_config3 \
	-t 0:.005:6 \
	-i 0:.005:6 \
	-u 5,3,4,2,1,5,6,4,4,2,1,5 \
	-o 4,2,3,1,1,4,5,3,3,1,1,4 \
	-y 20,20,20,20 \
	-n .00025 \
	> results/reg1/tsc3.txt
echo "Finished Running Toggle Switch N Script"
