use std::collections::HashMap;
use std::fs;
use std::io::Write;
use std::path::{Path, PathBuf};
use std::str::FromStr;

use failure;

use rust_htslib::bam;
use rust_htslib::bam::Read as BamRead;

use assign::*;
use barcode_group::*;
use depth::*;
use frag_purity::*;
use purity::*;

pub struct CLI {
    pub bambyname: String,
    pub outbase: String,

    pub minreads: String,
    pub minqual: String,
    pub minpurity: String,

    pub fwd_strand: Option<bool>,
    pub position: Option<String>,
}

pub struct Config {
    bambyname: String,
    outbase: PathBuf,
    min_reads: usize,
    min_qual: u8,
    min_purity: f64,
    #[allow(dead_code)]
    fwd_strand: Option<bool>,
    #[allow(dead_code)]
    position: Option<i32>,
}

impl Config {
    pub fn new(cli: &CLI) -> Result<Self, failure::Error> {
        Ok(Config {
            bambyname: cli.bambyname.to_string(),
            outbase: PathBuf::from(&cli.outbase),

            min_reads: usize::from_str(&cli.minreads)?,
            min_qual: u8::from_str(&cli.minqual)?,
            min_purity: f64::from_str(&cli.minpurity)?,

            fwd_strand: cli.fwd_strand,
            position: cli
                .position
                .as_ref()
                .map(|pos| i32::from_str(&pos))
                .transpose()?,
        })
    }

    pub fn outfile<T>(&self, file_suffix: T) -> PathBuf
    where
        T: std::convert::AsRef<std::ffi::OsStr>,
    {
        let mut filename = self
            .outbase
            .file_name()
            .map_or_else(|| std::ffi::OsString::new(), |str| str.to_os_string());
        filename.push(file_suffix);
        self.outbase.with_file_name(filename)
    }
}

pub fn bc_frag(config: Config) -> Result<(), failure::Error> {
    let mut input = if config.bambyname == "-" {
        bam::Reader::from_stdin()?
    } else {
        bam::Reader::from_path(Path::new(&config.bambyname))?
    };

    let mut fates: HashMap<Fate, usize> = HashMap::new();

    let mut good_assign = fs::File::create(config.outfile("barcode-good-assign.txt"))?;
    write!(good_assign, "barcode\ttarget\n")?;

    let mut good_assign_bed = fs::File::create(config.outfile("barcode-good-assign.bed"))?;

    let mut any_assign = fs::File::create(config.outfile("barcode-any-assign.txt"))?;
    write!(any_assign, "barcode\n")?;

    let mut bad_assign = fs::File::create(config.outfile("barcode-bad-assign.txt"))?;
    write!(bad_assign, "barcode\tDefect\tDetails\n")?;

    let mut depth_out = fs::File::create(config.outfile("barcode-depth.txt"))?;
    write!(depth_out, "{}\n", Depth::header())?;

    let mut purity_out = fs::File::create(config.outfile("barcode-purity.txt"))?;
    write!(purity_out, "{}\n", GrnaPurity::header())?;

    let mut fidelity_out = fs::File::create(config.outfile("barcode-fidelity.txt"))?;
    write!(fidelity_out, "{}\n", AssignMatch::header())?;

    let header = bam::Header::from_template(input.header());
    let header_view = bam::HeaderView::from_header(&header);

    let targets_result: std::result::Result<Vec<&str>, std::str::Utf8Error> = header_view
        .target_names()
        .iter()
        .map(|name_u8| ::std::str::from_utf8(name_u8))
        .collect();
    let targets = targets_result?;

    let barcode_groups = BarcodeGroups::new_with_read_names(&mut input)?;

    for barcode_group in barcode_groups {
        let (bc, qall) = barcode_group?;
        let bc_str = ::std::str::from_utf8(bc.as_slice()).unwrap_or("???");

        let (qvec, depth) = filter_by_quality(qall, config.min_qual);
        write!(depth_out, "{}\n", depth.line(bc_str))?;

        let fate = if depth.n_good() < config.min_reads {
            Fate::NoDepth
        } else {
            let purity = FragPurity::new(qvec.iter());

            write!(purity_out, "{}\t{}\n", bc_str, purity.line(&targets))?;

            if let Some(ap) = purity.primary_pos() {
                if purity.align_purity() < config.min_purity {
                    write!(
                        bad_assign,
                        "{}\tPurity\t{:.02}\n",
                        bc_str,
                        purity.align_purity()
                    )?;
                    Fate::NoPurity
                } else {
                    // Get full alignments!
                    write!(
                        good_assign,
                        "{}\t{}:{}{}\n",
                        bc_str,
                        ap.target(&targets),
                        ap.pos(),
                        ap.strand_chr()
                    )?;
                    write!(
                        good_assign_bed,
                        "{}\t{}\t{}\t{}\t{}\t{}\n",
                        ap.target(&targets),
                        ap.pos(),
                        ap.pos() as u32 + 1,
                        bc_str,
                        purity.n_primary(),
                        ap.strand_chr()
                    )?;
                    Fate::Good
                }
            } else {
                write!(
                    bad_assign,
                    "{}\tUnaligned\t{:.02}\n",
                    bc_str,
                    purity.read_purity()
                )?;
                Fate::NoMatch
            }
        };

        *fates.entry(fate).or_insert(0) += 1;
    }

    let mut fates_out = fs::File::create(config.outfile("barcode-fates.txt"))?;
    write!(
        fates_out,
        "Good\t{}\n",
        fates.get(&Fate::Good).unwrap_or(&0)
    )?;
    write!(
        fates_out,
        "No Match\t{}\n",
        fates.get(&Fate::NoMatch).unwrap_or(&0)
    )?;
    write!(
        fates_out,
        "Mixed\t{}\n",
        fates.get(&Fate::NoPurity).unwrap_or(&0)
    )?;
    write!(
        fates_out,
        "Too Few\t{}\n",
        fates.get(&Fate::NoDepth).unwrap_or(&0)
    )?;

    Ok(())
}

#[derive(Debug, Clone, Hash, PartialEq, Eq)]
enum Fate {
    Good,
    #[allow(dead_code)]
    NoAlignFidelity,
    #[allow(dead_code)]
    NoPosFidelity,
    NoMatch,
    NoPurity,
    NoDepth,
}
