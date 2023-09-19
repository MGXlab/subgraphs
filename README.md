This Rust application is designed to extract subgraphs from a larger graph based on two user-defined parameters: n and c. Here, n represents the fraction of nodes that must overlap with the target region, and c represents the number of nodes to include in the context.

This tool is particularly useful for visualizing subgraphs of LMEMs (Longest Matching Extension Matches) identified by Graphite. In complex graphs, it can be challenging to distinguish the region of interest due to clutter. By applying this simple filtering process, the application helps reduce complexity and generates a new GFA (Graphical Fragment Assembly) file that provides a clear and easy-to-visualize representation of the desired subgraph.

TODO: add CLI