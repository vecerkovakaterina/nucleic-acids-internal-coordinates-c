# Efficient computation of rigid base coordinates in C

Python module for computation of internal coordinates in DNA and RNA. Takes data from MD simulations in the form of a PDB file and creates output files containing base  coordiantes *shear, stretch, stagger, buckle, propeller, opening* and base step coordinates *shift, slide, rise, roll, tilt, twist*. The coordinates can be computed in either of the three definitions - 3DNA, Curves+ or cgDNA.

This project was created as a a part of my bachelor thesis.


# Compilation

To be able to use the program, compile the main file with any C compiler

    gcc main.c -o coordinates

# Usage
The usage of the program is displayed whenever required arguments are not provided

    USAGE: ./coordinates -t [3dna|curves|cgdna] -n number_of_snapshots -i /path/input/filename -o /path/output/filename

where input file name is a common name for a PDB files numbered 1 to n. And the output file name is a common name for the output files. For example

    $ ./coordinates -t 3dna -n 10 -i /my_dir/simulation -o /my_other_dir/coordinates

would compute base and base step coordinates according to the 3DNA definition for PDB files named simulation.pdb.1 through simulation.pdb.10.

# Output files
The output files are created for base and base step coordinates separately. For base coordinates, each base pair in the simulated oligomer will have its own output file, where tab separated columns represent *shear, stretch, stagger, buckle, propeller and opening* and rows represent its coordinates through the snapshots. The same goes for base step coordinates, each pair of base pairs will have its own file with columns *shift, slide, rise, roll, tilt, twist* and rows with coordinates for every snapshot. In the example above for an oligomer of length 5 with 10 snapshots the output files created are coordiantes_base_prm_1.out through coordiantes_base_prm_5.out each with 10 rows and coordinates_base_step_1.out through coordinates_base_step_4.out also with 10 rows.

# Licence
A [source file](http://lh3lh3.users.sourceforge.net/download/eigeng.c) that impelments functions to obtain eigenvectors from matrices distrubuted under the MIT lincence is used in the header file `lh3.h`.

The MIT License

Copyright (c) 1996 Frank Uhlig et al. 2009 Genome Research Ltd (GRL).

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.



