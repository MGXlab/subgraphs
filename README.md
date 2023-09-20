
# Graphite - subgraph extraction

This Rust application is designed to extract subgraphs from a larger graph based on two user-defined parameters: `n` and `c`. Here, `n` represents the fraction of nodes that must overlap with the target region, and `c` represents the number of nodes to include in the context.

This tool is particularly useful for visualizing subgraphs of LMEMs (Longest Matching Extension Matches) identified by [Graphite](https://github.com/rickbeeloo/GraphiteV2). In complex graphs, it can be challenging to distinguish the region of interest due to clutter. By applying this simple filtering process, the application helps reduce complexity and generates a new GFA (Graphical Fragment Assembly) file that provides a clear and easy-to-visualize representation of the desired subgraph.

For example:
![Example](test_data/subgraph_over_range.png)


## Documentation


```
USAGE:
    subgraphs [FLAGS] --context <context> --to <end> --graph <graph-file> --frac <node-fraq> --output <output> --seq <seq-file> --from <start> --target <target>

FLAGS:
    -d, --debug      
    -h, --help       Prints help information
    -V, --version    Prints version information

OPTIONS:
    -c, --context <context>     
    -t, --to <end>              
    -g, --graph <graph-file>    
    -n, --frac <node-fraq>      
    -o, --output <output>       
    -s, --seq <seq-file>        
    -f, --from <start>          
    -i, --target <target>
```
## Method

The subgraph extraction is quite straightforward, we first search the graph file for the target identifier (`-i`). From this path we extract the region from `-f` to `-t` and create a hashset of the nodes. We then read the graph file again, extract the path, and create a bool vector specifying whether the node is present in the hashset or not. We scan the bool vector for regions of true of at least size `n` (`n` = fraction of nodes query path). If the region satisfies `n` we now have a region that has nodes present in our query and is long enough to satisfy our criteria, this however may contain duplicates. We deduplicate the region and intersect it with the query node set to validate if `n` is actually satisfied. If so, we extract the context up to `c` to the left and right and write the matching path to a [GFA](http://gfa-spec.github.io/GFA-spec/GFA1.html). 

Ps. first ever Rust script, so optimizations are welcome :)
## Authors

- [@rickbeeloo](https://www.github.com/rickbeeloo)

