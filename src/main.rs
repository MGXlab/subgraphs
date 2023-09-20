use std::io::{BufWriter, Write};
use std::collections::HashSet;
use std::collections::HashMap;
use std::cmp;
use std::fs::File;
use structopt::StructOpt;

struct GFAWriter {
    writer: BufWriter<File>,
    sequences: HashMap<usize, String>,
}

impl GFAWriter {
    fn new(out_path: &str, sequences: HashMap<usize, String>) -> Result<Self, std::io::Error> {
        let f = File::create(out_path)?;
        let mut writer = BufWriter::new(f);

        // Write the header line
        writer.write_all(b"H\tVN:Z:2.0\n")?;

        Ok(GFAWriter {
            writer,
            sequences
        })
    }

    fn add_to_gfa(&mut self, node_path: &[usize]) -> Result<(), std::io::Error> {
        for &node_id in node_path {
            if let Some(sequence) = self.sequences.get(&node_id) {
                let segment_line = format!("S\t{}\t*\t{}\n", node_id, sequence.len());
                self.writer.write_all(segment_line.as_bytes())?;
            } else {
                eprintln!("Error: Sequence not found for node_id {}", node_id);
            }
        }

        let mut prev_node_id: Option<usize> = None;

        for &node_id in node_path {
            if let Some(prev_id) = prev_node_id {
                if let Some(overlap) = self.sequences.get(&prev_id).map(|seq| seq.len()) {
                    let link_line = format!("L\t{}\t+\t{}\t+\t{}M\n", prev_id, node_id, overlap);
                    self.writer.write_all(link_line.as_bytes())?;
                }
            }

            prev_node_id = Some(node_id);
        }

        Ok(())
    }

}

// No alloc reader from https://stackoverflow.com/questions/45882329/read-large-files-line-by-line-in-rust
mod my_reader {
    use std::{
        fs::File,
        io::{self, prelude::*},
    };

    pub struct BufReader {
        reader: io::BufReader<File>,
    }

    impl BufReader {
        pub fn open(path: impl AsRef<std::path::Path>) -> io::Result<Self> {
            let file = File::open(path)?;
            let reader = io::BufReader::new(file);

            Ok(Self { reader })
        }

        pub fn read_line<'buf>(
            &mut self,
            buffer: &'buf mut String,
        ) -> Option<io::Result<&'buf mut String>> {
            buffer.clear();

            self.reader
                .read_line(buffer)
                .map(|u| if u == 0 { None } else { Some(buffer) })
                .transpose()
        }
    }
}

fn find_tab(l: &str) -> usize {
    l.chars().position(|c| c == '\t').unwrap()
}

fn load_nodes(node_file: &str) -> HashMap<usize, String> {
    let mut seq_map: HashMap<usize, String> = HashMap::new();
    let mut reader = my_reader::BufReader::open(node_file).unwrap();
    let mut buffer = String::new();
    while let Some(line_result) = reader.read_line(&mut buffer) {
        let line = line_result.unwrap().trim();
        let tab_index: usize = find_tab(&line);
        let identifier: usize = line[..tab_index].parse::<usize>().unwrap();
        let seq: String = line[tab_index + 1..].to_string();
        seq_map.insert(identifier, seq);
    }
    seq_map
}


fn find_target<'a>(graph_path: &'a str, target_id:&str, start:usize, stop:usize) ->  Result<Vec<usize>, ()> {
    let mut reader = my_reader::BufReader::open(graph_path).unwrap();
    let mut buffer = String::new();

    while let Some(line_result) = reader.read_line(&mut buffer) {
        let line = line_result.unwrap().trim(); 
        let tab_index: usize = find_tab(&line);
        if &line[..tab_index] == target_id {
            let v: Vec<&str> = line[tab_index+1..].split(' ').collect();

            let lmem_target: Result<Vec<usize>, _> = v[start-1..stop]
            .iter()
            .map(|&s| s[..s.len() - 1].parse())
            .collect();
            return Ok(lmem_target.unwrap())
        }
    }
    Err(())
}


fn find_overlaps(graph_path: &str, target_path: &Vec<usize>, c: usize, n: f64, gfa_writer: &mut GFAWriter) {

    // Input
    let mut reader = my_reader::BufReader::open(graph_path).unwrap();
    let mut buffer = String::new();

    // Get hash set of target region
    let target_set: HashSet<usize> = target_path.iter().cloned().collect();

    let n_abs: usize = (target_set.len() as f64 * n) as usize;
    while let Some(line_result) = reader.read_line(&mut buffer) {
        let line = line_result.unwrap().trim();
        let tab_index = find_tab(&line);
        let identifier: &str = &line[..tab_index];

        let node_path: Vec<usize> = line[tab_index + 1..]
            .split(' ')
            .map(|s| s[..s.len() - 1].parse::<usize>())
            .collect::<Result<Vec<usize>, _>>()
            .unwrap();

        let match_mask: Vec<bool> = node_path.iter()
            .map(|x| target_set.contains(x))
            .collect();

        // We can find the chains satisfying n_abs and then extending bidirectional
        let mut i: usize = 0;
        while i < match_mask.len() {
    
            if match_mask[i] == true {
                
                let r_boundary = cmp::min(i + target_path.len() + c, node_path.len());
               
                let set_bools = match_mask[i..r_boundary].iter().filter(|&&x| x).count();

                if set_bools >= n_abs {
         
                    let candidate_set: HashSet<usize> = node_path[i..r_boundary].iter().cloned().collect();
                    let shared: usize = (&candidate_set & &target_set).len();

                    if shared >= n_abs {
                        gfa_writer
                        .add_to_gfa(&node_path[i..r_boundary])
                        .expect("Error adding node path to GFA");
                        println!("Found for {} at {} with {} / {}; set bools {}", identifier, i, shared, target_set.len(), set_bools);
                        i = r_boundary;
                    }
                }
            }
            i += 1;
        }

        // Clear the buffer for the next line
        buffer.clear();
    }
}


#[derive(Debug, StructOpt)]
#[structopt(name = "Graphtie - subgraphs", about = "Extracting an LMEM subgraph")]
struct Opt {
    #[structopt(short, long)]
    debug: bool,

    #[structopt(short = "g", long = "graph", about = "Path to graph file, cf.seq")]
    graph_file: String,

    #[structopt(short = "s", long = "seq", about = "Path to node sequence file, cf.seg")]
    seq_file: String,

    #[structopt(short="n", long = "frac", about = "Fraction of nodes in LMEM to be present")]
    node_fraq: f64,

    #[structopt(short="i", long = "target", about = "Identifier of target path")]
    target: String,

    #[structopt(short="f", long = "from", about = "Start position on path")]
    start: usize,

    #[structopt(short="t", long = "to", about = "End position on path")]
    end: usize,

    #[structopt(short="c", long = "context", about = "Context to include, in nodes")]
    context: usize,

    #[structopt(short="o", long = "output", about = "Output path to write GFA to")]
    output: String,

}

fn main() {
    let opt = Opt::from_args();
    let target_path = find_target(&opt.graph_file, &opt.target, opt.start, opt.end).unwrap();    
    let node_map = load_nodes(&opt.seq_file);
    let mut gfa_writer = GFAWriter::new(&opt.output, node_map.clone()).expect("Error creating GFA file");

    find_overlaps(&opt.graph_file, &target_path, opt.context, opt.node_fraq, &mut gfa_writer);

}
