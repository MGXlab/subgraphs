use std::io::{BufWriter, Write};
use std::collections::HashSet;
use std::collections::HashMap;
use std::cmp;
use std::fs::File;

struct GFAWriter {
    writer: BufWriter<File>,
    sequences: HashMap<usize, String>,
}

impl GFAWriter {
    fn new(out_path: &str, sequences: HashMap<usize, String>) -> Result<Self, std::io::Error> {
        let f = File::create(out_path)?;
        let mut writer = BufWriter::new(f);

        // Write the header line (optional but recommended)
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

fn main() {
    let f = "/media/codegodz/LaCie/Paper1/graph_timed.cf_seq".to_string();
    let node_file = "/media/codegodz/LaCie/Paper1/graph_timed.cf_seg".to_string();
    let target_id = "Reference:277_Sequence:CP047082".to_string();
    let target_start = 9153;
    let target_end = 9257;
    let c = 100;
    let n:f64 = 0.9;
    let output_path = "output.gfa".to_string();

    let target_path = find_target(&f, &target_id, target_start, target_end).unwrap();
    
    let node_map = load_nodes(&node_file);

    let mut gfa_writer = GFAWriter::new(&output_path, node_map.clone()).expect("Error creating GFA file");


    find_overlaps(&f, &target_path, c, n, &mut gfa_writer);

}
