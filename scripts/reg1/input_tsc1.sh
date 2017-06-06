#~/bin/bash
echo "Running Toggle Switch N Script:"
	./bin/reg_exp \
	-s toggle_switch_config1 \
	-t 0:.005:6 \
	-i 0:.005:6 \
	-u 5,3,2,1,5,6,4,2,1,5 \
	-o 2,1,1,1,1,2,1,1,1,2 \
	-y 20,20,20,20 \
	-n .005 \
	> results/reg1/tsc1.txt
echo "Finished Running Toggle Switch N Script"
