#~/bin/bash
echo "Running Toggle Switch N Script:"
	./bin/reg_exp \
	-s toggle_switch_config2 \
	-t 0:.005:6 \
	-i 0:.005:6 \
	-u 5,3,4,1,5,6,4,4,1,5 \
	-o 2,1,1,1,1,2,1,1,1,2 \
	-y 20,20,20,20 \
	-n 0.0005 \
	> results/reg1/tsc2.txt
echo "Finished Running Toggle Switch N Script"
