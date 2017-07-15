#~/bin/bash
echo "Running Lorenz No Chaos Script:"
for i in {1..2}; do
	echo $i
	./bin/main \
	-s lorenz \
	-t 0:.005:2\
	-i 0:.005:2 \
	-u 10,20,2.666 \
	-o 5,15,5  \
	-y 1,1,1 \
	-k 100 \
	-p 120 \
	-d 1 \
	-n $i \
	-r \
	> results/lor_results/lornc$i.txt
done
echo "Finished Running Lorenz No Chaos Script"
