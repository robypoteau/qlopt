# Qlopt

[![Build Status](https://travis-ci.org/jgoldfar/qlopt.svg?branch=master)](https://travis-ci.org/jgoldfar/qlopt)
[![Build status](https://ci.appveyor.com/api/projects/status/7oe6e10sjbqt0g6w?svg=true)](https://ci.appveyor.com/project/jgoldfar/qlopt)


This is a software package developed to pursue numerical calculations reported 
in the paper “Identification of parameters in systems biology” by 
Ugur G. Abdulla & Roby Poteau, submitted to Mathematical Biosciences.

> Professor Ugur G. Abdulla,Ph.D.,Dr.Sci.  
> Department of Mathematical Sciences  
> Florida Institute of Technology  
> Melbourne, FL  
> Email: [abdulla@fit.edu](abdulla@fit.edu)

> Roby Poteau  
> Ph.D. Student in Operations Research  
> Department of Mathematical Sciences  
> Florida Institute of Technology  
> Email: [rpoteau2010@my.fit.edu](rpoteau2010@my.fit.edu)


## Prerequisites

* [GNU Make](https://www.gnu.org/software/make/)

* [Eigen 3](http://eigen.tuxfamily.org/index.php?title=Main_Page)

* [Boost 1.65+](https://www.boost.org)

### OSX Build Instructions

The version of Eigen and Boost installed by [Homebrew](https://brew.sh/) is sufficient to build qlopt on an OSX platform:

```shell
brew install eigen boost
```

You'll need the command-line build tools available as well, though those are a Homebrew build requirement.

### Ubuntu/Debian
If you're running a new enough Debian/Ubuntu system, you can install those libraries through your package manager.
The versions available on Ubuntu Bionic (18.04 LTS) or Debian Buster are known to work, and can be installed by

```shell
sudo apt-get install libeigen3-dev libboost-math-dev build-essential
```

On Ubuntu Xenial or Debian Stretch, Boost is too old, so you'll need to build it from source (see below) and on distributions older than that, both Boost and Eigen must be built from source.
You'll need `build-essential` no matter which version of Debian/Ubuntu you're using.

### Build Dependencies from Source

To build [Eigen 3.3.7](https://github.com/eigenteam/eigen-git-mirror) from source local to your project, run

```shell
git clone https://github.com/eigenteam/eigen-git-mirror.git
mkdir -p eigen-git-mirror/build && cd eigen-git-mirror/build
git checkout tags/3.3.7 -b eigen-3.3.7
cmake -DCMAKE_INSTALL_PREFIX=/path/to/project-root/usr ..
make install
```

*Note*: Modify the prefix above to point to your project root, and follow the modified build instructions below.

To install [Boost 1.69.0](https://dl.bintray.com/boostorg/release/1.69.0/source/boost_1_69_0.tar.gz) local to your project, run

```shell
wget https://dl.bintray.com/boostorg/release/1.69.0/source/boost_1_69_0.tar.gz
tar xf boost_1_69_0.tar.gz
cd boost_1_69_0
./bootstrap.sh --prefix=/path/to/project-root/usr
./b2 -d0 -q --with-math --with-test --with-system install
```

*Note*: Modify the prefix above to point to your project root.

If you install Boost or Eigen this way, compile Qlopt with

```shell
make all CXXFLAGS_EXTRA="-I`pwd`/usr/include"
```

### Installation and Usage

Download or clone the repository into a folder and run `make` command in the project's root directory.

When you run `make` it will create executables in the `bin` folder from source files 
in the `examples` folder. You can use them as a guide to apply this code to your
own examples.

When you create your own example, just drop it into the examples
folder and run the `make` in the root directory; this will create an
executable which you can run from the `bin` folder. 

## Questions
If you have any questions feel free to e-mail us at:
abdulla@fit.edu; rpoteau2010@my.fit.edu

## License
This project is licensed under the GNU General Public License v3.0 License - see the [LICENSE.md](LICENSE.md) file for details
