## Description
This repository computes generators and relations of the rational cohomology of Hilbert schemes from *link*. 

## Run computations
1. To compute all monomials in terms of basis elements run 
> ``python recursion_formula.py``

It is possible to pass a dimension and save path by adding ``--dimension d --save_path ./path-to/file.pkl``

2. To compute minimal relations run
> ``python compute_relations.py``

It is possible to pass a dimension, load path and save path by adding ``--dimension d --load_path ./path-to/monomials.pkl --save_path ./path-to/relations.pkl`` By default the relations are saved in a txt file with TeX syntax.