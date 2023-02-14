# Anonymous-ISD-Sieving
Anonymous github folder for double-blind peer review for the paper "New sieving-style ISD"

# isd-sieving
Many versions of the algorithms were tried: (main.cc, new_main.cc, impl1.cc, impl2.cc).

The `main.cc' refers to the basic (but not trivial implementation). The impl1 and impl2 refers to other versions not mentioned in the paper.

Currently using **new_main.cc** for the implementation. To compile: 

``` g++ -o new new_main.cc vector.cc misc.cc ```

``` ./new ```

Parameters to change in new_main.cc include $length = k+ l$; $nbrParities = \ell$, $p, 2p$, and $nbrSamples = M = \frac{2}{q} \cdot \delta$.

