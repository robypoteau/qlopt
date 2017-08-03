#~/bin/bash
echo "Running Lorenz Script:"
for i in {2..2}; do
	echo $i
	./bin/prog -s lorenz -t 0:.005:2 -i 0:.005:2 -u 10,28,2.666 -o 5,15,5 -y 1,1,1 -k 100 -p 26 d 1 -n $i -b > results/lor_results/lor$i.txt
done
echo "Finished Running Lorenz Script"
