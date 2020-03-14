extern crate barcode_assign;
#[macro_use]
extern crate clap;
extern crate failure;

use std::path::PathBuf;

use clap::{App, Arg};

use barcode_assign::bc_grna::*;

fn main() {
    let matches = App::new("bc-grna")
        .version("1.0")
        .author("Nick Ingolia <ingolia@berkeley.edu>")
        .about("Barcode / guide RNA assignment")
        .arg(
            Arg::with_name("bambyname")
                .short("b")
                .long("bam-by-name")
                .value_name("SORTED-BY-NAME-BAM")
                .help("BAM format alignment sorted by name")
                .takes_value(true)
                .required(true),
        )
        .arg(
            Arg::with_name("outbase")
                .short("o")
                .long("outbase")
                .value_name("OUTBASE")
                .help("Output base filename")
                .takes_value(true)
                .required(true),
        )
        .arg(
            Arg::with_name("alignstart")
                .long("align-start")
                .value_name("START")
                .help("Alignment starting position (0-based)")
                .takes_value(true)
                .required(true),
        )
        .arg(
            Arg::with_name("minreads")
                .long("min-reads")
                .value_name("MINREADS")
                .help("Minimum number of reads for good assignment")
                .takes_value(true)
                .default_value("3"),
        )
        .arg(
            Arg::with_name("minqual")
                .long("min-qual")
                .value_name("MINQUAL")
                .help("Minimum (in any position) quality score to include read")
                .takes_value(true)
                .default_value("30"),
        )
        .arg(
            Arg::with_name("mintargetpurity")
                .long("min-target-purity")
                .value_name("PURITY")
                .help("Minimum fraction of reads that must align to the same guide")
                .takes_value(true)
                .default_value("0.901"),
        )
        .arg(
            Arg::with_name("reverse")
                .long("reverse")
                .help("Alignment reverse strand"),
        )
        .get_matches();

    let config = Config {
        bowtie_bam: PathBuf::from(matches.value_of("bambyname").unwrap()),
        align_start: value_t!(matches.value_of("alignstart"), usize).unwrap_or_else(|e| e.exit()),

        is_reverse: matches.is_present("reverse"),

        out_base: PathBuf::from(matches.value_of("outbase").unwrap()),
        
        min_reads: value_t!(matches.value_of("minreads"), usize).unwrap_or_else(|e| e.exit()),
        min_qual: value_t!(matches.value_of("minqual"), u8).unwrap_or_else(|e| e.exit()),
        min_purity: value_t!(matches.value_of("mintargetpurity"), f64).unwrap_or_else(|e| e.exit()),
    };

    if let Err(ref e) = barcode_to_grna(&config) {
        println!("error: {}", e);

        ::std::process::exit(1);
    }
}

