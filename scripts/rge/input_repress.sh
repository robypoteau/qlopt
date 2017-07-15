#~/bin/bash
echo "Running Repressilator Script:"
./bin/reg_guess_exp \
    -s repressilator \
	-t 0:.01:10 \
	-i 0:.01:10 \
	-u 216,0.216,2,0.2 \
	-g 250,.4,3,.4 \
	-o 100,.75,.75,.75 \
	-y 5,1,10,1,15,1 \
	-n 0.0009 \
	> results/rge/repress_rng.txt
echo "Finished Running Repressilator Script"
