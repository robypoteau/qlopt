#~/bin/bash
echo "Running Repressilator Script:"
for i in {1..1}; do
	echo $i
	./bin/main \
	-s general_repressilator \
	-i 0:.02:10 \
	-u 216,0.216,2,0.2,0.2,216,0.216,2,0.2,0.2,216,0.216,2,0.2,0.2 \
	-o 175,.5,1,.5,.5,175,.5,1,.5,.5,175,.5,1,.5,.5 \
	-g 250,.4,3,.4,.4,250,.4,3,.4,.4,250,.4,3,.4,.4 \
	-y 5,1,10,1,15,1 \
	-k 1 \
	-p 22 \
	-d 1 \
	-r 1 \
	-l 0.000001 \
	> results/repress_results/repress$i.txt
done
echo "Finished Running Repressilator Script"

#with lsq sinusoidal divs 15
