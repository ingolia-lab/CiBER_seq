extern crate barcode_assign;
extern crate clap;

use clap::{App, Arg};

use barcode_assign::bc_seqs::*;

fn main() {
    let matches = App::new("bc-seqs")
        .version("1.0")
        .author("Nick Ingolia <ingolia@berkeley.edu>")
        .about("Match R1 barcode sequences with R2 insert sequences")
        .arg(
            Arg::with_name("barcodes")
                .short("b")
                .long("barcodes")
                .value_name("BARCODE-FQ")
                .help("FastQ file of barcode sequences")
                .takes_value(true)
                .required(true),
        )
        .arg(
            Arg::with_name("sequences")
                .short("s")
                .long("sequences")
                .value_name("SEQUENCE-FQ")
                .help("FastQ file of insert sequences")
                .takes_value(true)
                .required(true),
        )
        .arg(
            Arg::with_name("outbase")
                .short("o")
                .long("outbase")
                .value_name("OUTBASE")
                .help("Output filename base (barcoded sequences in OUTBASE_barcoded.fq)")
                .takes_value(true)
                .required(true),
        )
        .arg(
            Arg::with_name("no-bclist")
                .long("no-barcode-list")
                .help("Do not write list of barcodes")
                .conflicts_with("bclist"),
        )
        .arg(
            Arg::with_name("bclist")
                .short("l")
                .long("barcode-list")
                .value_name("BARCODES-TXT")
                .help("Filename for list of barcodes (default OUTBASE-barcodes.txt)")
                .takes_value(true)
                .conflicts_with("no-bclist"),
        )
        .arg(
            Arg::with_name("no-bcfreq")
                .long("no-barcode-freq")
                .help("Do not write table of barcode frequencies")
                .conflicts_with("bcfreq"),
        )
        .arg(
            Arg::with_name("bcfreq")
                .short("f")
                .long("barcode_freq")
                .value_name("BARCODE-FREQS-TXT")
                .help("Filename for barcode frequency table (default OUTBASE-barcode-freqs.txt)")
                .takes_value(true)
                .conflicts_with("no-bcfreq"),
        )
        .arg(
            Arg::with_name("neighborhood")
                .short("n")
                .long("neighborhood")
                .value_name("NBHD_BASE")
                .help("Group barcodes into single-mismatch neighborhoods")
                .takes_value(true),
        )
        .get_matches();

    let outbase = matches.value_of("outbase").unwrap();
    let barcodes = matches
        .value_of("bclist")
        .map_or(outbase.to_string() + "-barcodes.txt", |s| s.to_string());
    let barcode_freqs = matches
        .value_of("bcfreq")
        .map_or(outbase.to_string() + "-barcode-freqs.txt", |s| {
            s.to_string()
        });

    let config = Config {
        barcode_fastq: matches.value_of("barcodes").unwrap().to_string(),
        sequ_fastq: matches.value_of("sequences").unwrap().to_string(),
        out_fastq: outbase.to_string() + "_barcoded.fq",
        out_barcodes: if matches.is_present("no-bclist") {
            None
        } else {
            Some(barcodes)
        },
        out_barcode_freqs: if matches.is_present("no-bcfreq") {
            None
        } else {
            Some(barcode_freqs)
        },
        neighborhood: matches.value_of("neighborhood").map(|s| String::from(s)),
    };

    match bc_seqs(config) {
        Ok(_) => (),
        Err(e) => panic!(e),
    }
}
