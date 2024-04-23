use std::collections::HashSet;
use std::collections::HashMap;
use std::cmp;
use structopt::StructOpt;
use std::path::Path;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::PathBuf;
use std::fs::{File, remove_file};

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
        let seq: String = line[tab_index..].to_string();
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

struct GFA {
    link_writer: BufWriter<File>,
    sequences: HashMap<usize, String>,
    observed_segments: HashSet<usize>,
    output: PathBuf,
    k: usize,
    tmp: String,
}

impl GFA {

    fn new(output: &str, sequences: HashMap<usize, String>, tmp: &str, k: usize) -> Result<Self, std::io::Error>  {
        // First create the tmp output files for the GFA 
        let tmp_link_path = Path::new(tmp).join("tmp_links.txt");
        println!("Tmp link path: {:?}", tmp_link_path);
        let link_writer = BufWriter::new(File::create(tmp_link_path).expect("Failed to create tmp link file"));
        let observed_segments: HashSet<usize> = HashSet::new();
        Ok(GFA {
            link_writer,
            observed_segments,
            output: Path::new(tmp).join(output.to_string()),
            sequences,
            k,
            tmp: tmp.to_string()
        })
    }

    fn split_node(node_str: &str) -> (usize, char) {
        let node_id:usize = node_str[..node_str.len() - 1].parse::<usize>().unwrap();
        let prev_dir: char = *&node_str.chars().last().unwrap();
        (node_id, prev_dir)
    }

    fn add_path(&mut self, node_path: &[&str]) {
        assert!(node_path.len() > 1, "Single node path not supported yet");
        for (i, node_str) in node_path.iter().enumerate() {
            if i == 0 {continue;};
            let (node_id, node_sign) = GFA::split_node(node_str);
            self.observed_segments.insert(node_id);
            let (prev_id, prev_sign ) = GFA::split_node(node_path[i-1]);
            let link_line = format!("L\t{}\t{}\t{}\t{}\t{}M\n", prev_id, prev_sign, node_id, node_sign, self.k-1); // k-1 overlap I suppose for all?
            self.link_writer.write_all(link_line.as_bytes()).expect("Couldn't write link to file");
        }
    }

    fn finalize(&mut self) {
        // We first open the output file, dump the segments then read the link file and copy them over
        let mut writer = BufWriter::new(File::create(&self.output).expect("Couldn't create output file"));
        println!("Writing segments...");
        for seen_segment in &self.observed_segments {
            let seq: &str = self.sequences.get(&seen_segment).expect("Node not present");
            let segment_line = format!("S\t{}\t{}\n", seen_segment, seq);
            writer.write_all(segment_line.as_bytes()).expect("Couldnt write link to file");
        }
        println!("Copying links to output file...");
        let tmp_link_path = Path::new(&self.tmp).join("tmp_links.txt");
        println!("Tmp link: {:?}", tmp_link_path);
        self.link_writer.flush().unwrap();
        let source_file = File::open(&tmp_link_path).expect("Could not open tmp link file");
        let source_reader = BufReader::new(source_file);
        for line in source_reader.lines() {
            let line = line.unwrap();
            writer.write_all(line.as_bytes()).unwrap();
            writer.write_all(b"\n").unwrap(); 
        }
        println!("Cleaning tmp");
        remove_file(&tmp_link_path).unwrap();
        println!("Done writing GFA!😊");
    }

}



fn find_overlaps(c: &Config, target_path: &Vec<usize>,  sequences: &HashMap<usize, String>, output: &str) {
    // Input
    let mut reader = my_reader::BufReader::open(&c.graph_file).unwrap();
    let mut buffer = String::new();

    // GFA output 
    let mut gfa = GFA::new(output, sequences.clone(), &c.tmp, c.k ).unwrap();


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
                        

                        gfa.add_path(&node_str_path[l_boundary..r_boundary]);

                        // Check if we also should calculate the actual genomic positions 
                        if c.write_coords {
                            println!("ID: {}, from: {}, to: {}", identifier, l_boundary, r_boundary);
                            let (start, end) = calc_location(&node_path, l_boundary, r_boundary, sequences, c.k).unwrap();
                            let s = format!("{}\t{}\t{}\t{}\t{}\n", identifier, shared, target_set.len(), start, end);

                            if let Some(writer) = &mut coord_handle {
                                writer.write_all(s.as_bytes()).unwrap();
                            }
                        }
                        
                        // we only have to keep track of the node colors when we write them in the end
                        // helps to visualize in bandage
                        if c.write_colors {
                            for node_id in node_path[l_boundary..r_boundary].iter() {
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
                    for identifier in identifiers.iter() {
                        let s = format!("{}\t{}\n", node_id, identifier);
                        writer.write_all(s.as_bytes()).unwrap();
                    }
                }                
            }            
        }

        // Clear the buffer for the next line
        buffer.clear();
    }
    println!("Finalizing!");
    gfa.finalize();
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

    #[structopt(short="m", long = "tmp", about = "Tepm dir")]
    tmp: String,

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

struct Config {
    graph_file: String, 
    context_size: usize, 
    node_fraction: f64, 
    write_coords: bool,
    coord_path: String,
    write_colors: bool, 
    tmp: String,
    color_path: String,
    k: usize,
}


fn main() {
    let opt = Opt::from_args();
    let target_path = find_target(&opt.graph_file, &opt.target, opt.start, opt.end).unwrap();    
    println!("Target path: {:?}", target_path);
    let node_map = load_nodes(&opt.seq_file);

    let conf: Config = Config{
        graph_file: opt.graph_file,
        node_fraction: opt.node_fraq,
        context_size: opt.context, 
        write_coords: opt.coords, 
        k : opt.kmer.unwrap_or(0),
        coord_path: opt.coord_path.unwrap_or("".to_string()), 
        write_colors: opt.colors,
        color_path: opt.color_path.unwrap_or("".to_string()),
        tmp: opt.tmp,
    };
    

    find_overlaps(&conf, &target_path,  &node_map, &opt.output);

}


