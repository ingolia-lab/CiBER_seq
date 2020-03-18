extern crate clap;
extern crate failure;

extern crate barcode_assign;

use std::io;
use std::io::Write;
use std::process;

use clap::{App, Arg};

use barcode_assign::bc_frag::*;

fn main() {
    match wrapper() {
        Err(e) => {
            io::stderr().write(format!("{}\n", e).as_bytes()).unwrap();
            process::exit(1);
        }
        _ => (),
    };
}

fn wrapper() -> Result<(), failure::Error> {
    let cli = get_cli()?;
    let config = Config::new(&cli)?;
    bc_frag(config)
}

fn get_cli() -> Result<CLI, failure::Error> {
    let matches = App::new("bc-frag")
        .version("0.0.0")
        .author("Nick Ingolia <ingolia@berkeley.edu>")
        .about("Compute barcode-to-fragment assignment from genomic alignment")
        .long_about("ZZZ")
        .arg(
            Arg::with_name("bambyname")
                .short("b")
                .long("bam-by-name")
                .value_name("SORTED_BY_NAME.BAM")
                .help("BAM format input file, sorted by read name")
                .takes_value(true)
                .required(true),
        )
        .arg(
            Arg::with_name("outbase")
                .short("o")
                .long("output")
                .value_name("OUTBASE")
                .help("Filename base for output")
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
        .get_matches();

    Ok(CLI {
        bambyname: matches.value_of("bambyname").unwrap().to_string(),
        outbase: matches.value_of("outbase").unwrap().to_string(),
        minreads: matches.value_of("minreads").unwrap().to_string(),
        minqual: matches.value_of("minqual").unwrap().to_string(),
        minpurity: matches.value_of("mintargetpurity").unwrap().to_string(),
        fwd_strand: None,
        position: None,
    })
}
