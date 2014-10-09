FACILITY LOCATION PROBLEM

The code is a dirty one that approximates a solution to the Capacitated Facility Location Problem. No dataset is provided but a small example is incorporated in the script (see ZIMPL examples directory, from https://code.google.com/p/python-zibopt/source/browse/trunk/examples/facility-location.py?r=170).

It has been tested against benchmark dataset examples found in a coursera.org training. The current algorithm deviates from an acceptable minimal approximation when the file gets bigger than 100-200 customers.
Despite of that, the approximation to a minimal acceptable benchmark was still acceptable for the cases subjected to test (about 5%-10% worst)* when not better. It could be 0%-30% above the optimum.

(*)OBS: measured as benchmark/observed

It consists in a greedy comparison using a simple heuristic that uses LOCAL (neighbour) information only. The three best heuristics are included as commented lines in the code.
Depending of the problem, a heuristic could probe better than the rest.

The code is some steps away from the correct solution to solve all the example datasets.

For information about this code please contact me at varikvi(alt)yahoo(dot)com