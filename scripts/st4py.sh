#~/bin/bash
export OPTFLAGS=$(/usr/bin/python2.7-config --cflags)
echo $OPTFLAGS
export OPTLIBS=$(/usr/bin/python2.7-config --ldflags)
echo $OPTLIBS
make
