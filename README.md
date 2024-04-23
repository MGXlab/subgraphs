# Cuttlefish - subgraph extraction
This Rust application is designed to extract subgraphs from a larger graph based on two user-defined parameters: `n` and `c`. Here, `n` represents the fraction of nodes that must overlap with the target region, and `c` represents the number of nodes to include in the context.

This tool is particularly useful for visualizing subgraphs of LMEMs (Longest Matching Extension Matches) identified by [Graphite](https://github.com/rickbeeloo/GraphiteV2). In complex graphs, it can be challenging to distinguish the region of interest due to clutter. By applying this simple filtering process, the application helps reduce complexity and generates a new GFA (Graphical Fragment Assembly) file that provides a clear and easy-to-visualize representation of the desired subgraph.

For example:
![Example](test_data/subgraph_over_range.png)

## Installation
1. Install Rust, first check if rust is installed using `rustc --verion`, if not, install it ([read here]((https://www.rust-lang.org/tools/install))).
2. Then `git clone https://github.com/rickbeeloo/subgraphs`
3. Then move into the dir, `cd subgraphs` and run `cargo build --release`
4. The binary is now at `./target/release/subgraphs`

## Quickstart
First built a graph using [Cuttlefish](https://github.com/COMBINE-lab/cuttlefish). Then run it using something like:
`./target/release/subgraphs -g test_data/graph/test_graph.cf_seq -s test_data/graph/test_graph.cf_seg -i Reference:1_Sequence:seq1 -f 4 -t 8 -n 1.0 -c 4 -o test.gfa  --coords -p coords.tsv -x -a node_colors.csv  -k 7 -m test_data/`
- `-g` The `cf_seq` file 
- `-s` the `cf_seg` file
- `-i` the reference id of which we extract the subpath
- `-f` **from** at which index (0-based) we start in the path (included)
- `-t`**to** till which index the subpath extends (included)
- `-n` the fraction of nodes a path should traverse to end up in the final GFA, so `1.0` means all nodes in the subpath, but `0.8` means `80%` of the subpath nodes 
- `-c` how many nodes left and right of the subpath to include in the final graph. Usually something like `10-100` makes sense to get a clue of the context
- `-o` where to write the .GFA of the subgraph
- `--coords` + `-p`, to write the matching coordinates in the other sequences, this can be ommitted but might be nice to extract all subsequences for further analysis
- `-x` + `-a` We write a `.csv` that we can load in bandage that, for each node, has a list of the genomes we found it in. We can use that to color or just look at.


## Method

The subgraph extraction is quite straightforward, we first search the graph file for the target identifier (`-i`). From this path we extract the region from `-f` to `-t` and create a hashset of the nodes. We then read the graph file again, extract the path, and create a bool vector specifying whether the node is present in the hashset or not. We scan the bool vector for regions of true of at least size `n` (`n` = fraction of nodes query path). If the region satisfies `n` we now have a region that has nodes present in our query and is long enough to satisfy our criteria, this however may contain duplicates. We deduplicate the region and intersect it with the query node set to validate if `n` is actually satisfied. If so, we extract the context up to `c` to the left and right and write the matching path to a [GFA](http://gfa-spec.github.io/GFA-spec/GFA1.html). 

Ps. first ever Rust script, so optimizations are welcome :)
## Authors

- [@rickbeeloo](https://www.github.com/rickbeeloo)

