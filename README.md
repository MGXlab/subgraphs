
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
## Authors

- [@rickbeeloo](https://www.github.com/rickbeeloo)

