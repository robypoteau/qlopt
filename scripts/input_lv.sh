#~/bin/bash
echo "Running Lotka Volterra Script:"
for i in {1..2}; do
	echo $i
	./bin/main \
	-s lotka_volterra \
	-t 0:.001:1 \
	-i 0:.001:1 \
	-u 1,2,1 \
	-o 10,9,10 \
	-y 1,1 \
	-k 1 \
	-p 7 \
	-d 1 \
	> results/lv_results/lv$i.txt
done
echo "Finished Running LV Script"
#-u .48,.026,.93 -o .05,.05,.05
