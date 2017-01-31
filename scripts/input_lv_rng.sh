#~/bin/bash
echo "Running Lotka Volterra Script:"
./bin/prog -s lotka_volterra -t 0:.1:1 -i 0:.1:1 -u .48,.026,.93 -o 1,1,1 -y 35,5 -k 1 -p 7 -d 1 -r > results/lv_results/lv_rng.txt
echo "Finished Running LV Script"
