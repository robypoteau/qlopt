#~/bin/bash
echo "Running Toggle Switch N Script:"
for i in {1..1}; do
	echo $i
	./bin/main \
	-s toggle_switch_config2 \
	-t 0:.005:6  \
	-i 0:.005:6  \
	-u 5,3,4,1,5,6,4,4,1,5 \
	-o 2,2,2,2,2,2,2,2,2,2 \
	-y 20,20,20,20 \
	-k 1 \
	-p 10 \
	-n $i \
	-r \
	> results/ts_results/tsc2_$i.txt
done
echo "Finished Running Toggle Switch N Script"
