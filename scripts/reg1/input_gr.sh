#~/bin/bash
echo "Running General Repressilator Script:"
./bin/reg_exp \
	-s general_repressilator \
	-t 0:.01:10 \
	-i 0:.01:10 \
	-u 216,0.216,2,0.2,0.2,216,0.216,2,0.2,0.2,216,0.216,2,0.2,0.2 \
	-o 175,.5,1,.5,.5,1575,.5,1,.5,.5,,175,.5,1,.5,.5 \
	-y 5,1,10,1,15,1 \
	-n 0.0000001 \
	 > results/reg1/gr_rng.txt
echo "Finished Running General Repressilator Script"
