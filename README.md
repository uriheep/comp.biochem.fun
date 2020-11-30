# funDFT

1. Have you tried to use statistical distributions for the initial
values of coefficients in the development of molecular orbitals on the basis set functions?
Would you like to match the obtained results with the known data and to study convergence/divergence
of iterations?
2. Have you tried to generate massive amounts of simulation results and feed it to a machine learning
engine to infer a new heuristic / law / molecular force-field?
3. Have you dreamed of cracking and hacking molecular modeling without fear to 'break something' or 'being looking ignorant' ?

If you are reading this line, Welcome to **funDFT** ! :)

Dream BIG and have FUN!

**funDFT** is an implementation of the SCF method
with the focus on a reasonable balance between clear code and performance considerations.
Intended to serve demonstrative / didactic / research purposes.

## Benchmarking

to be reported

## Dependences:

1. C++ [ Armadillo ]( http://arma.sourceforge.net/docs.html ) library for linear algebra
2. C++ [ *Rapid*JSON ]( https://rapidjson.org/ ) library for parsing JSON files with the basis set functions from the [ Basis Set Exchange ]( https://www.basissetexchange.org/ ) 
3. C++ [ readcif ]( https://www.rbvi.ucsf.edu/chimerax/docs/devel/bundles/mmcif/mmcif_cpp/readcif_cpp/docs/api.html ) library for parsing CIF files with molecular structures
4. C++ [ Boost ]( https://www.boost.org/ ) library for unit tests ( only `unit_test_framework` is needed for unit tests; compiled statically )
5. C++ STL
6. [ scons ]( https://scons.org/ ) for building

## Build:

1. `scons` to build all the code in this project

## Example:

to be reported


