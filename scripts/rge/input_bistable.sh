#~/bin/bash
echo "Running Bistable Script:"
./bin/reg_guess_exp \
    -s bistable_switch \
    -t 0:.01:3 \
    -i 0:.01:3 \
    -u 150,3.2,2 \
    -g 135,2.0,1.5 \
    -o 200,1,1 \
    -y 25,25 \
    -n 0.00005 \
    > results/rge/bs_rng.txt
echo "Finished Running BS Script"
# ~best without regs
#./bin/reg_guess_exp -s bistable_switch -t 0:.1:3 -i 0:.1:3 -u 150,3.2,2 -g 135,2.2,1.5 -o 115,1.5,1 -y 25,25 -n 0.0001 > results/rge/bs_rng.txt
#-o 200,4,3 \
