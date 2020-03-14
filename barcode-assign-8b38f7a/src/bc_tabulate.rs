extern crate barcode_assign;
extern crate clap;

use std::io::Write;

use clap::{App, Arg};

use barcode_assign::bc_tabulate::*;

fn main() {
    let matches = App::new("bc-tabulate")
        .version("1.0")
        .author("Nick Ingolia <ingolia@berkeley.edu>")
        .about("Tabulate barcode counts from many different samples")
        .arg(
            Arg::with_name("inputs")
                .index(1)
                .multiple(true)
                .takes_value(true)
                .required(true),
        )
        .arg(
            Arg::with_name("output")
                .short("o")
                .long("output")
                .value_name("OUTPUT.TXT")
                .help("Output filename for tab-delimited count table")
                .takes_value(true)
                .required(true),
        )
        .arg(
            Arg::with_name("mintotal")
                .short("m")
                .long("mintotal")
                .value_name("#READS")
                .help("Minimum total read count across all samples to report barcode")
                .takes_value(true),
        )
        .arg(
            Arg::with_name("minsamples")
                .short("n")
                .long("minsamples")
                .value_name("#SAMPLES")
                .help("Minimum number of samples with barcode to report")
                .takes_value(true),
        )
        .arg(
            Arg::with_name("mininsample")
                .long("mininsample")
                .value_name("#READS")
                .help("Minimum read count in at least one sample to report barcode")
                .takes_value(true),
        )
        .arg(
            Arg::with_name("omitfile")
                .long("omitfile")
                .value_name("OMITTED.TXT")
                .help("Output filename for omitted barcodes")
                .takes_value(true),
        )
        .get_matches();

    fn parse_int(value: &str) -> usize {
        value
            .parse::<usize>()
            .unwrap_or_else(|e| panic!("Bad whole number argument {:?}: {}", value, e))
    }

    let cli = CLI {
        inputs: matches
            .values_of("inputs")
            .unwrap()
            .map(String::from)
            .collect(),
        output: matches.value_of("output").unwrap().to_string(),
        mintotal: matches.value_of("mintotal").map(parse_int),
        minsamples: matches.value_of("minsamples").map(parse_int),
        mininsample: matches.value_of("mininsample").map(parse_int),
        omitfile: matches.value_of("omitfile").map(String::from),
    };

    match cli.run() {
        Ok(_) => (),
        Err(e) => {
            std::io::stderr()
                .write(format!("{}\n", e).as_bytes())
                .unwrap();
            std::process::exit(1);
        }
    }
}
