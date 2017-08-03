#~/bin/bash
echo "Running Lorenz Script:"
for i in {1..1}; do
	echo $i
	./bin/main \
	-s lorenz \
	-t 0:.005:2 \
	-i 0:.005:2 \
	-u 10,28,2.666 \
	-o 5,15,5  \
	-y 1,1,1 \
	-k 100 \
	-p 38 \
	-d 1 \
	-n $i \
	> results/lor_results/lor$i.txt
done
echo "Finished Running Lorenz Script"
