#~/bin/bash
echo "Running Toggle Switch N Script:"
	./bin/reg_guess_exp \
	-s toggle_switch_noise \
	-t 0:.005:6 \
	-i 0:.005:6 \
	-u 5,3,2,1,5,6,4,2,1,5 \
	-g 4,2,3,.5,4,5,3,4,.5,4 \
	-o 1,2,3,2,2,1,2,3,2,2 \
	-y 20,20,20,20 \
	-n 0.00005 \
	> results/rge/tsten_range.txt
echo "Finished Running Toggle Switch N Script"
