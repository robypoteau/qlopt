#~/bin/bash
echo "Running General Repressilator Script:"
./bin/reg_guess_exp \
	-s general_repressilator \
	-t 0:.01:10 \
	-i 0:.01:10 \
	-u 216,0.216,2,0.2,0.2,216,0.216,2,0.2,0.2,216,0.216,2,0.2,0.2 \
	-o 175,.5,1,.5,.5,175,.5,1,.5,.5,175,.5,1,.5,.5 \
	-g 250,.4,3,.4,.4,250,.4,3,.4,.4,250,.4,3,.4,.4 \
	-y 5,1,10,1,15,1 \
	-n 0.00001 \
	> results/rge/gr_rng.txt
echo "Finished Running General Repressilator Script"
