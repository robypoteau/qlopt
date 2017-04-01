#~/bin/bash
echo "Running Lorenz Script:"
./bin/reg_guess_exp \
	-s lorenz \
	-t 0:.005:2 \
	-i 0:.005:2  \
	-u 10,28,2.666 \
	-g 8,24,2 \
	-o 28,49,7 \
	-y 1,1,1 \
	-n 0.00005 \
	> results/rge/lor_rng.txt
echo "Finished Running Lorenz Script"
