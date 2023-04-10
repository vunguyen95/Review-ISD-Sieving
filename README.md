# Anonymous-ISD-Sieving
Anonymous github folder for double-blind peer review for the paper "New sieving-style ISD"

# isd-sieving
Many versions of the algorithms were tried: (main.cc, new_main.cc, impl1.cc, impl2.cc (not presented in this folder)).

The `main.cc' refers to the basic (but not trivial implementation).

Currently using **new_main.cc** for the implementation. To compile: 

``` g++ -o new new_main.cc vector.cc misc.cc ```

``` ./new ```

Parameters to change in new_main.cc include $length = k+ l$; $nbrParities = \ell$, $p, 2p$, and $nbrSamples = M = \frac{2}{q} \cdot \delta$.

The folder Submit_to_Pub contains estimate and success_probability script (in Python). Results obtained in the tables are presented in the txt files. Numbers 2 and 5 refer to the c_label values used in the estimation. 

