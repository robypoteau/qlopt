#~/bin/bash
echo "Running Lotka Volterra Script:"
for i in {9..10}; do
	echo $i
	./bin/prog -s lotka_volterra -t 0:.005:2 -i 0:.005:2 -u .48,.026,.93 -o .05,.05,.05 -y 35,4 -k 100 -n $i -p 3 -b > results/lv_results/lv$i.txt
done
echo "Finished Running LV Script"