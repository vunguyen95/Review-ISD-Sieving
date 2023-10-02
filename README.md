# ISD-Sieving
Our github repo for peer review process for the paper "New sieving-style ISD". A few words about this repository: this is a basic implementation for the Algorithm Merge_Set (Algorithm 2, p.10). Our goal is to verify if we can find a `target' random vector using this routine. We generate a random vector and parity check matrix for each run and observe if this vector can be produced using our described algorithm. More importantly, the memory requirement should be consistent with the theoretical analysis. We stated in the paper that the prospect of implementing a full-scale experiment (medium-size instance, including the outer iteration) are what we are aiming in the near future.

# isd-sieving
Many versions of the algorithms were tried: (main.cc, new_main.cc, impl1.cc, impl2.cc (not presented in this folder)).

The `main.cc' refers to the basic (but not trivial implementation).

Currently using **new_main.cc** for the implementation. To compile: 

``` g++ -o new new_main.cc vector.cc misc.cc ```

``` ./new ```

Parameters to change in new_main.cc include $length = k+ l$; $nbrParities = \ell$, $p, 2p$, and $nbrSamples = M = \frac{2}{q} \cdot \delta$.

The folder Submit_to_Pub contains estimate and success_probability script (in Python). Results obtained in the tables are presented in the txt files. Numbers 2 and 5 refer to the c_label values used in the estimation. 

