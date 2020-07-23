# Qlopt

[![Build Status](https://travis-ci.org/robypoteau/qlopt.svg?branch=master)](https://travis-ci.org/robypoteau/qlopt)
<!---[![Build status](https://ci.appveyor.com/api/projects/status/7oe6e10sjbqt0g6w?svg=true)](https://ci.appveyor.com/project/jgoldfar/qlopt)-->


This is a software package developed to pursue numerical calculations reported 
in the paper “Identification of parameters in systems biology” by 
Ugur G. Abdulla & Roby Poteau, submitted to Mathematical Biosciences.

Professor Ugur G. Abdulla,Ph.D.,Dr.Sci.  
Department of Mathematical Sciences  
Florida Institute of Technology  
Melbourne, FL  
Email:abdulla@fit.edu

Roby Poteau  
Ph.D. Student in Operations Research  
Department of Mathematical Sciences  
Florida Institute of Technology  
Email:rpoteau2010@my.fit.edu


### Prerequisites

I use [GNU Make](https://www.gnu.org/software/make/) which is bundled in build-essential package

```
sudo apt-get install build-essential
```

[Eigen 3.3.7](http://eigen.tuxfamily.org/index.php?title=Main_Page)

```
git clone https://github.com/eigenteam/eigen-git-mirror.git eigen_source
mkdir eigen_build
cd eigen_build
cmake ../eigen_source
make install            # may need sudo to run this
```

[Boost 1.69.0](https://dl.bintray.com/boostorg/release/1.69.0/source/boost_1_69_0.tar.gz)

```
wget https://dl.bintray.com/boostorg/release/1.69.0/source/boost_1_69_0.tar.gz
tar xf boost_1_69_0.tar.gz
cd boost_1_69_0
sudo ./bootstrap.sh
sudo ./b2 install
```

[SUNDIALS](https://computing.llnl.gov/projects/sundials/download/sundials-4.1.0.tar.gz)

```
wget https://computing.llnl.gov/projects/sundials/download/sundials-4.1.0.tar.gz
tar zxf sundails-4.1.0.tar.gz
cd sundails-4.1.0
mkdir build/
cd build
cmake ..
make
make install
```

### Installation
Make sure Eigen is installed. Download or clone the repository into a 
folder and run `make` command in the project's root directory 
(you can also use `sudo make` if the command doesn't work).

### Usage
When you run `make` it will create executables in the bin folder from source files 
in the examples folder. You can use them as a guide to apply this code to your
own examples.
When you create your own examples, just drop it into the examples
folder and run the make file in the root directory this will create an
executable which you can run from the bin folder. 

To run the executable, you need libparamid in your system's library
```
# At QLOPT root directory
cd build
sudo cp libparamid.* /usr/lib/
sudo ldconfig
```

## Questions
If you have any questions feel free to e-mail us at:
abdulla@fit.edu; rpoteau2010@my.fit.edu

## License
This project is licensed under the GNU General Public License v3.0 License - see the [LICENSE.md](LICENSE.md) file for details

## Notes
I developed in a debian Linux environment.
