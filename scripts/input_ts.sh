#~/bin/bash
echo "Running Toggle Switch Script:"
for i in {1..1}; do
	echo $i
	./bin/main \
	-s toggle_switch_config1 \
	-i 0:.06:6 \
	-u 1,2,1,1,1,2,1,1,10,10 \
	-o .5,1,.5,.5,.5,1,.5,.5,5,5 \
	-g 1,2,1,1,1,2,1,1,10,10 \
	-y 5,2,5,2 \
	-k 1 \
	-p 17 \
	-r 0 \
	-d 1 \
	-l 0.0 \
	> results/ts_results/ts$i.txt
done
echo "Finished Running Toggle Switch Script"	
