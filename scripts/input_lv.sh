#~/bin/bash
echo "Running Lotka Volterra Script:"
for i in {12..13}; do
	echo $i
	./bin/prog -s lotka_volterra -t 0:.005:4 -i 0:.005:4 -u .48,.026,.93 -o .05,.05,.05 -y 35,4 -k 10 -n $i -p 5 -b -d 2 > results/lv_results/lv$i.txt
done
echo "Finished Running LV Script"