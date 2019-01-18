# Qlopt

[![Build Status](https://travis-ci.org/jgoldfar/qlopt.svg?branch=master)](https://travis-ci.org/jgoldfar/qlopt)
[![Build status](https://ci.appveyor.com/api/projects/status/7oe6e10sjbqt0g6w?svg=true)](https://ci.appveyor.com/project/jgoldfar/qlopt)


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

[Eigen 3.3.4](http://bitbucket.org/eigen/eigen/get/3.3.4.tar.bz2) - the 
installation instruction can be found on there website 
[eigen.tuxfamily.org/](http://eigen.tuxfamily.org/).

I also use [GNU Make](https://www.gnu.org/software/make/).

### Installation
Make sure Eigen is installed. Download or clone the repository into a 
folder and run make command in the project's root directory.

### Usage
When you run `make` it will create executables in the bin folder from source files 
in the examples folder. You can use them as a guide to apply this code to your
own examples.
When you create your own examples, just drop it into the examples
folder and run the make file in the root directory this will create an
executable which you can run from the bin folder. 

## Questions
If you have any questions feel free to e-mail us at:
abdulla@fit.edu; rpoteau2010@my.fit.edu

## License
This project is licensed under the GNU General Public License v3.0 License - see the [LICENSE.md](LICENSE.md) file for details

## Notes
I developed in a debian Linux environment.
