#~/bin/bash
echo "Running Lotka Volterra Script:"
for i in {1..15}; do
	echo $i
	./bin/prog -s lotka_volterra -t 0:.1:1 -i 0:.1:1 -u .48,.026,.93 -o .05,.05,.05 -y 35,4 -k 1 -p 7 -d 1 > results/lv_results/lv$i.txt
done
echo "Finished Running LV Script"
