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

    fn add_to_gfa(&mut self, node_str_path: &[&str]) -> Result<(), std::io::Error> {
        for node_id_str in node_str_path {
            let node_id:usize = node_id_str[..node_id_str.len() - 1].parse::<usize>().unwrap();
            if let Some(sequence) = self.sequences.get(&node_id) {
                let segment_line = format!("S\t{}\t*\t{}\n", node_id, sequence.len());
                self.writer.write_all(segment_line.as_bytes())?;
            } else {
                eprintln!("Error: Sequence not found for node_id {}", node_id);
            }
        }

        let mut prev_node_id: Option<&str> = None;

        for node_id_str in node_str_path {
            if let Some(prev_id) = prev_node_id {
                let node_id:usize = node_id_str[..node_id_str.len() - 1].parse::<usize>().unwrap();
                if let Some(overlap) = self.sequences.get(&node_id).map(|seq| seq.len()) {

                    // Extract the id and traversal dir
                    let prev_id: &str = &prev_node_id.unwrap()[..prev_node_id.unwrap().len()-1];
                    let prev_dir: &str = &prev_node_id.unwrap().chars().last().unwrap().to_string();
                    let node_dir: &str = &node_id_str.chars().last().unwrap().to_string();

                    let link_line = format!("L\t{}\t{}\t{}\t{}\t{}M\n", prev_id, prev_dir, node_id, node_dir, overlap);
                    self.writer.write_all(link_line.as_bytes())?;
                }
            }

            prev_node_id = Some(&node_id_str);
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
        let identifier = &line[..tab_index];
        //println!("{} - {}", identifier, target_id);
        if identifier == target_id {
            println!("Target found");
            let v: Vec<&str> = line[tab_index+1..].split(' ').collect();

            let lmem_target: Result<Vec<usize>, _> = v[start..stop]
            .iter()
            .map(|&s| s[..s.len() - 1].parse())
            .collect();
            return Ok(lmem_target.unwrap())
        }
    }
    Err(())
}

fn calc_location(node_path: &Vec<usize>, from_idx: usize, to_idx: usize, sequences: &HashMap<usize, String>, k: usize) ->  Result<(usize, usize), ()> {
    println!("From: {}, to: {}", from_idx, to_idx);
    let mut genomic_position = 1;
    let mut start = 0;
    for (i, node_id) in node_path.iter().enumerate() {
        let seq_len = sequences.get(&node_id).unwrap().len();
        if i == from_idx {
            start = genomic_position;
        } else if i == to_idx {
            let stop = genomic_position + seq_len - 1;
            return Ok((start, stop));
        }
        genomic_position += seq_len - k + 1;
    }
    println!("Complete fail");
    return Err(());
}

struct Config {
    graph_file: String, 
    context_size: usize, 
    node_fraction: f64, 
    write_coords: bool,
    coord_path: String,
    write_colors: bool, 
    color_path: String,
    k: usize,
}

fn find_overlaps(c: &Config, target_path: &Vec<usize>, gfa_writer: &mut GFAWriter,  sequences: &HashMap<usize, String>) {
    // Input
    let mut reader = my_reader::BufReader::open(&c.graph_file).unwrap();
    let mut buffer = String::new();

    let mut color_handle: Option<BufWriter<File>> = if c.write_colors {
        Some(BufWriter::new(File::create(&c.color_path).unwrap()))
    } else {
        None
    };

    let mut coord_handle: Option<BufWriter<File>> = if c.write_coords {
        Some(BufWriter::new(File::create(&c.coord_path).unwrap()))
    } else {
        None
    };


    // We might to dump the node colors based on our selection so we can color them on different metadata
    // For this lets for now just hashmap the node ids to the identifiers
    let mut color_map: HashMap<usize, Vec<String>> = HashMap::new();

    // Get hash set of target region
    let target_set: HashSet<usize> = target_path.iter().cloned().collect();

    let n_abs: usize = (target_set.len() as f64 * c.node_fraction) as usize;
    while let Some(line_result) = reader.read_line(&mut buffer) {
        let line = line_result.unwrap().trim();
        let tab_index = find_tab(&line);
        let identifier: &str = &line[..tab_index];

        let node_str_path: Vec<&str> =  line[tab_index + 1..].split(' ').collect(); // We can keep the original GFA formatting

        let node_path: Vec<usize> = node_str_path.iter().map(|s| s[..s.len() - 1].parse::<usize>())
            .collect::<Result<Vec<usize>, _>>()
            .unwrap();

        let match_mask: Vec<bool> = node_path.iter()
            .map(|x| target_set.contains(x))
            .collect();

        // We can find the chains satisfying n_abs and then extending bidirectional
        let mut i: usize = 0;
        while i < match_mask.len() {
    
            if match_mask[i] == true {
                
                let r_boundary = cmp::min(i + target_path.len(), node_path.len()-1);
                                                
                let set_bools = match_mask[i..r_boundary].iter().filter(|&&x| x).count();

                if set_bools >= n_abs {
         
                    let candidate_set: HashSet<usize> = node_path[i..r_boundary].iter().cloned().collect();
                    let shared: usize = (&candidate_set & &target_set).len();

                    if shared >= n_abs {

                        // Now extract with the context, c
                        let r_boundary = cmp::min(i + target_path.len() + c.context_size, node_path.len()-1);
                        let l_boundary = cmp::max(i as i32 - c.context_size as i32, 0) as usize;
                        println!("Max left: {}",l_boundary);

                        gfa_writer
                        .add_to_gfa(&node_str_path[l_boundary..r_boundary])
                        .expect("Error adding node path to GFA");

       
                        // Check if we also should calculate the actual genomic positions 
                        if c.write_coords {
                            let (start, end) = calc_location(&node_path, l_boundary, r_boundary, sequences, c.k).unwrap();
                            let s = format!("{}\t{}\t{}\t{}\t{}\n", identifier, shared, target_set.len(), start, end);

                            if let Some(writer) = &mut coord_handle {
                                writer.write_all(s.as_bytes()).unwrap();
                            }
                        }
                        
                       // println!("{}: {}", identifier, node_path[i..r_boundary].len());
                        // we only have to keep track of the node colors when we write them in the end
                        if c.write_colors {
                            for node_id in node_path[i..r_boundary].iter() {
                                color_map.entry(*node_id) 
                                    .or_insert(Vec::new())
                                    .push(identifier.to_string());
                            }
                        }
                        
                        i = r_boundary;
                    }
                }
            }
            i += 1;
        }

        if c.write_colors {
            if let Some(writer) = &mut color_handle {
                for (node_id, identifiers) in &color_map {
                    //println!("node id: {} with identifiers: {}", node_id, identifiers.len());
                    for identifier in identifiers.iter() {
                        //println!("{}", identifier);
                        let s = format!("{}\t{}\n", node_id, identifier);
                        writer.write_all(s.as_bytes()).unwrap();
                    }
                }                
            }            
        }

        // Clear the buffer for the next line
        buffer.clear();
    }
}


#[derive(Debug, StructOpt)]
#[structopt(name = "Graphtie - subgraphs", about = "Extracting an LMEM subgraph")]
struct Opt {

    // flags
    #[structopt(short, long)]
    debug: bool,

    // other args
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

    // Sub path positions

    #[structopt(short = "y", long = "coords")]
    coords: bool,

    #[structopt(short = "k", long = "kmer")]
    kmer: Option<usize>,

    #[structopt(short = "p", long = "coord_path")]
    coord_path: Option<String>,

    // Colors
    
    #[structopt(short = "x", long = "colors")]
    colors: bool,

    #[structopt(short = "a", long = "color_path")]
    color_path: Option<String>,


}

fn main() {
    let opt = Opt::from_args();
    let target_path = find_target(&opt.graph_file, &opt.target, opt.start, opt.end).unwrap();    
    println!("Target path: {:?}", target_path);
    let node_map = load_nodes(&opt.seq_file);
    let mut gfa_writer = GFAWriter::new(&opt.output, node_map.clone()).expect("Error creating GFA file");
    

    let conf: Config = Config{
        graph_file: opt.graph_file,
        node_fraction: opt.node_fraq,
        context_size: opt.context, 
        write_coords: opt.coords, 
        k : opt.kmer.unwrap_or(0),
        coord_path: opt.coord_path.unwrap_or("".to_string()), 
        write_colors: opt.colors,
        color_path: opt.color_path.unwrap_or("".to_string()),
    };

    find_overlaps(&conf, &target_path, &mut gfa_writer, &node_map);

}


