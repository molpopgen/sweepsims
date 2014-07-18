	sweepsims - A collection of programs to simulate selective sweeps under different scenarios



  Copyright (C) 2006-2014 Kevin Thornton

  sweepsims is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

Comments are welcome.

	- Kevin Thornton <krthornt@uci.edu>

#Dependencies

libsequence (http://github.com/molpopgen/libsequence)

#Installation

./configure

make

sudo make install

##If dependent libraries are in non-standard locations.  For example, "/opt":

CXXFLAGS=-I/opt/include LDFLAGS="$LDFLAGS -L/opt/lib" ./configure

make 

sudo make install

##Installing somewhere other than /usr/local/bin

./configure --prefix=/path/to/where/you/want/it

For example,

./configure --prefix=$HOME

make 

make install

will result in msstats being in ~/bin.

#Usage

Over the years, I've done lots of modeling of selective sweeps, and this repo has been my "think-tank" for developing the necessary simulations. It is here on git now. To figure out what a program does, try to run it w/no arguments and see what it asks you to do. If a program segfaults when run w/no arguments, the reason is that the program is not checking for no arguments. You'd need to read the code to figure out the arguments. (Such programs have not been used for research, and may be in various states of disarray.)

Some of the programs have been used in various publications over the years, and it is good form to cite the appropriate literature if and when you use them. The programs and the corresponding citations are:

rsweep_stochCG and rsweep_stochCGdist were developed for, and used in:

```
@article{Jensen:2008eo, 
	author = {Jensen, Jeffrey D and Thornton, Kevin R and Andolfatto, Peter}, 
	title = {{An approximate bayesian estimator suggests strong, recurrent selective sweeps in Drosophila.}}, 
	journal = {PLoS Genetics},
	year = {2008}, 
	volume = {4}, 
	number = {9}, 
	pages = {e1000198} }
```
sweep_der and sweep_der_randomX were developed for, and used in:

```
@article{Thornton:2007du,
         author = {Thornton, Kevin R and Jensen, Jeffrey D}, 
	 title = {{Controlling the false-positive rate in multilocus genome scans for selection.}}, 
	 journal = {Genetics}, 
	 year = {2007}, 
	 volume = {175}, 
	 number = {2}, 
	 pages = {737--750}, 
	 month = feb }
```
The programs listed above are (reasonably) well-tested. As with all simulations of selection in the coalescent, there are often strong assumptions in the models that are usually not enforced in the code. This means that the user can inadvertently do "dumb" things by inputting parameter values that are not sensible given the assumptions, and the program will happily do non-rigorous things. Examples include the programs based on deterministic trajectories of beneficial mutations. These models require that 2Ns >= 1000 at the very least. Modeling weaker selection is not a good idea, and one should use the programs based on stochastic trajectories. I strongly encourage the users to go back and refresh their memories with the original literature whenever possible. Caveat emptor, and all that.