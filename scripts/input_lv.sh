#~/bin/bash
echo "Running Lotka Volterra Script:"
for i in {0..0}; do
	echo $i
	./bin/main \
	-s lotka_volterra \
	-i 0:1:10 \
	-u .48,.026,.93 \
	-o .05,.05,.05 \
	-y 35,4 \
	-k 1 \
	-p 7 \
	-d 1 \
	-r 0 \
	-l 0.0001 \
	> results/lv_results/lv$i.txt
done
echo "Finished Running LV Script"
#	-u 1,2,1 \ -o 10,9,10 \
