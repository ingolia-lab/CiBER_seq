extern crate barcode_assign;
extern crate bio;
extern crate clap;
extern crate failure;

use clap::{App, Arg};

use barcode_assign::bc_collapse::CLI;

fn main() {
    let matches = App::new("bc-collapse")
        .version("1.0")
        .author("Nick Ingolia <ingolia@berkeley.edu>")
        .about("Collapse barcode sequences")
        .arg(
            Arg::with_name("input")
                .short("i")
                .long("input")
                .value_name("BARCODES.TXT")
                .help("Text file of barcode sequences")
                .takes_value(true)
                .required(true),
        )
        .arg(
            Arg::with_name("output_base")
                .short("o")
                .long("outbase")
                .value_name("OUTPUT_BASE")
                .help("Base name for output files")
                .takes_value(true)
                .required(true),
        )
        .get_matches();

    let cli = CLI {
        input: matches.value_of("input").unwrap().to_string(),
        output_base: matches.value_of("output_base").unwrap().to_string(),
    };

    match cli.run() {
        Ok(_) => (),
        Err(e) => panic!(e),
    }
}

