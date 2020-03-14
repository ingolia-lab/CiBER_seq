use std::collections::HashMap;
use std::io::BufRead;
use std::io::BufReader;
use std::io::Read;
use std::path::{Path,PathBuf};

use neighborhood::*;

pub struct CLI {
    pub input: String,
    pub output_base: String,
}

impl CLI {
    pub fn output_filename(&self, name: &str) -> PathBuf {
        let base_ref: &Path = self.output_base.as_ref();
        let mut namebase = base_ref
            .file_name()
            .map_or(std::ffi::OsString::new(), std::ffi::OsStr::to_os_string);
        namebase.push(name);
        base_ref.with_file_name(namebase)
    }

    pub fn run(&self) -> Result<(), failure::Error> {
        let reader: Box<dyn Read> = if self.input == "-" {
            Box::new(std::io::stdin())
        } else {
            Box::new(std::fs::File::open(&self.input)?)
        };
        let mut barcode_reader = BufReader::new(reader);

        let barcode_counts = CLI::count_barcodes(&mut barcode_reader)?;

        let nbhds_raw = Neighborhood::gather_neighborhoods(barcode_counts);
        let nbhds: Vec<_> = nbhds_raw.into_iter().map(|n| n.into_sorted()).collect();

        SortedNeighborhood::write_tables(&self.output_base, nbhds.iter())?;

        Ok(())
    }

    pub fn count_barcodes<R: BufRead>(barcode_reader: &mut R) -> Result<HashMap<Vec<u8>, usize>, failure::Error> {
        let mut barcode_counts = HashMap::new();
        
        for line_res in barcode_reader.lines() {
            let line = line_res?;
            
            let barcode = line.into_bytes();
            let barcode_count = barcode_counts.entry(barcode).or_insert(0);
            *barcode_count += 1;
        }
        
        Ok(barcode_counts)
    }
}

