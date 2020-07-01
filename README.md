# LaWeCSE
Computes BBP-MCISs; allows to skip vertices as in Largest Weight Common Subtree Embedding

1) download and build OGDF; at least an out-of-source build into a 'Release'-folder, standard settings should suffice
(https://ogdf.uos.de/download/ ; Catalpa 2020-02-09) , see below
2) set correct OGDF_Path in CMakeLists.txt of 'LaWeCSE', this should be the OGDF-directory containing the README.md
3) if you did not create a Release-Folder as ouf-of-source build, set that path accordingly in the 1., 2., and 4. row from below
4) run cmake, then make

For a list of parameters type ./LaWeCSE --help

Example:
./LaWeCSE -c example.fog 1 2

An exemplary weight function is 'weight_function', to be added with '-l weight_function'.

The labels are _string_ labels, the same is true for the graph file

An exemplary graph file 'example.fog' is also attached. The file format is 3 lines per graph:
1) #name, used in output of LaWeCSE; number of nodes N; number of edges M
2) label node 1; ... ; label node N 
3) M triples of: start node; end node; edge label

It also accepts .gml files.

------
With this software you may compute maximum common induced subgraphs between two input graphs, as described in 'Finding Largest Common Substructures of Molecules in Quadratic Time', http://link.springer.com/chapter/10.1007/978-3-319-51963-0_24

Additionally this software supports skipped vertices, cf. 'Largest Weight Common Subtree Embeddings with Distance Penalties', MFCS 2018, http://drops.dagstuhl.de/opus/frontdoor.php?source_opus=9636

In case the input is not outerplanar, we compute a mapping between blocks of which at least one is not outerplanar based on a clique approach (see 'An algorithm for reporting maximal c-cliques' from Frédéric Cazals, Chinmay Karande)

---
Exemplary OGDF Setup

To install OGDF, first unpack it (link above), then within the OGDF directory perform the following steps
1) mkdir Release
2) cd Release
3) ccmake ../
4) type c, then g (make sure 'release' type build is selected)
5) cmake .
6) _instead of_ performing steps 3 to 5, you may just type cmake ../ for a standard build.
7) make

Then, this Release directory is the one to specify in the LaWeCSE cmake as above.

