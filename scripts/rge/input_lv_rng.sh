#~/bin/bash
echo "Running Lotka Volterra Script:"
./bin/reg_guess_exp -s lotka_volterra -t 0:.01:1 -i 0:.01:1 -u 1,2,1 -g 5,6,5 -o 10,50,120 -y 1,1 -n .000009 > results/lv_results/lv_rng.txt
echo "Finished Running LV Script"
