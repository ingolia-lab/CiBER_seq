use std::fs;
use std::collections::HashMap;

use std::io::Write;
use std::path::PathBuf;
use std::str;

use rust_htslib::bam;
use rust_htslib::bam::Read;

use assign::*;
use barcode_group::*;
use depth::*;
use purity::*;

#[derive(Debug)]
pub struct Config {
    pub bowtie_bam: PathBuf,
    pub align_start: usize,
    pub is_reverse: bool,
    pub out_base: PathBuf,
    pub min_reads: usize,
    pub min_qual: u8,
    pub min_purity: f64,
}

impl Config {
    pub fn outfile<T>(&self, file_suffix: T) -> PathBuf
    where
        T: std::convert::AsRef<std::ffi::OsStr>,
    {
        let mut filename = self
            .out_base
            .file_name()
            .map_or_else(|| std::ffi::OsString::new(), |str| str.to_os_string());
        filename.push(file_suffix);
        self.out_base.with_file_name(filename)
    }
}

pub fn barcode_to_grna(config: &Config) -> Result<(), failure::Error> {
    let mut fates: HashMap<Fate, usize> = HashMap::new();

    let mut good_assign = fs::File::create(config.outfile("barcode-grna-good.txt"))?;
    write!(good_assign, "barcode\tguide\n")?;
    let mut bad_assign = fs::File::create(config.outfile("barcode-bad-assign.txt"))?;
    write!(bad_assign, "barcode\tDefect\tDetails\n")?;

    let mut depth_out = fs::File::create(config.outfile("barcode-depth.txt"))?;
    write!(depth_out, "{}\n", Depth::header())?;

    let mut purity_out = fs::File::create(config.outfile("barcode-purity.txt"))?;
    write!(purity_out, "{}\n", GrnaPurity::header())?;

    let mut fidelity_out = fs::File::create(config.outfile("barcode-fidelity.txt"))?;
    write!(fidelity_out, "{}\n", AssignMatch::header())?;

    let mut bam_reader = bam::Reader::from_path(&config.bowtie_bam)?;
    let header = bam::Header::from_template(bam_reader.header());
    let header_view = bam::HeaderView::from_header(&header);

    let targets_result: std::result::Result<Vec<&str>, std::str::Utf8Error> = header_view
        .target_names()
        .iter()
        .map(|name_u8| str::from_utf8(name_u8))
        .collect();
    let targets = targets_result?;

    let barcode_groups = BarcodeGroups::new_with_read_names(&mut bam_reader)?;

    for barcode_group in barcode_groups {
        let (bc, qall) = barcode_group?;
        let bc_str = str::from_utf8(bc.as_slice()).unwrap_or("???");

        let (qvec, depth) = filter_by_quality(qall, config.min_qual);
        write!(depth_out, "{}\n", depth.line(bc_str))?;

        let fate = if depth.n_good() < config.min_reads {
            Fate::NoDepth
        } else {
            let purity = GrnaPurity::new(qvec.iter())?;
            write!(purity_out, "{}\n", purity.line(bc_str))?;

            if purity.target_purity() < config.min_purity {
                write!(
                    bad_assign,
                    "{}\tPurity\t{:.02}\n",
                    bc_str,
                    purity.target_purity()
                )?;
                Fate::NoPurity
            } else {
                if let ReadAssign::Match(assign) = purity.primary_assign() {
                    write!(fidelity_out, "{}\n", assign.line(bc_str, &targets)?)?;

                    if (assign.pos() == config.align_start as i64)
                        && (assign.is_reverse() == config.is_reverse)
                        && assign.is_cigar_perfect()
                        && assign.is_md_perfect()
                    {
                        write!(good_assign, "{}\t{}\n", bc_str, assign.target(&targets))?;
                        Fate::Good
                    } else {
                        write!(bad_assign, "{}\tFidelity\t{}\n", bc_str, assign.field())?;
                        Fate::NoFidelity
                    }
                } else {
                    write!(bad_assign, "{}\tUnaligned\tN/A\n", bc_str)?;
                    Fate::NoMatch
                }
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
        "Bad Match\t{}\n",
        fates.get(&Fate::NoFidelity).unwrap_or(&0)
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
    NoFidelity,
    NoMatch,
    NoPurity,
    NoDepth,
}
