#~/bin/bash
echo "Running Lotka Volterra Script:"
./bin/reg_guess_exp2 -s lotka_volterra -t 0:.01:1 -i 0:.01:1 -u 1,2,1 -g 6,5,6 -o 12,15,12 -y 1,1 > results/lv_results/lv_rng2.txt
echo "Finished Running LV Script"
