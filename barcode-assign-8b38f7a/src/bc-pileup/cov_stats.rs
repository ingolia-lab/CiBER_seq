use std::cmp::min;
use std::collections::HashSet;
use std::fmt::{self, Display, Formatter};
use std::iter::repeat;

use rust_htslib::bam;

use aln_pos::{AlnPos, AlnPosCons};
use offset_vector::OffsetVector;

#[derive(Debug, Clone, Hash, PartialEq, Eq, PartialOrd, Ord, Copy)]
pub enum CoverClass {
    Good,
    Gapped,
    Heterog,
}

impl Display for CoverClass {
    fn fmt(&self, f: &mut Formatter) -> fmt::Result {
        write!(
            f,
            "{}",
            match *self {
                CoverClass::Good => "Good",
                CoverClass::Gapped => "Gapped",
                CoverClass::Heterog => "Heterog",
            }
        )
    }
}

#[derive(Debug, Clone)]
pub struct Cover {
    wildtype: usize,
    mutant: usize,
    heterog: usize,
    none: usize,
}

impl Cover {
    pub fn new<'a, I>(aln_cons: I, _covstart: usize, _covend: usize) -> Self
    where
        I: Iterator<Item = (usize, &'a AlnPosCons)>,
    {
        let mut cover = Cover {
            wildtype: 0,
            heterog: 0,
            none: 0,
            mutant: 0,
        };
        for (_pos, apc) in aln_cons {
            if apc.is_wildtype() {
                cover.wildtype += 1;
            } else if apc.is_mutant() {
                cover.mutant += 1;
            } else if apc.is_uncovered() {
                cover.none += 1;
            } else {
                panic!("Unclassifiable AlnPosCons {:?}\n", apc);
            }

            if apc.is_heterog() {
                cover.heterog += 1;
            }
        }
        cover
    }

    pub fn wildtype(&self) -> usize {
        self.wildtype
    }
    pub fn heterog(&self) -> usize {
        self.heterog
    }
    pub fn none(&self) -> usize {
        self.none
    }
    pub fn mutant(&self) -> usize {
        self.mutant
    }

    pub fn classify(&self, max_none: usize, max_heterog: usize) -> CoverClass {
        if self.none > max_none {
            CoverClass::Gapped
        } else if self.heterog > max_heterog {
            CoverClass::Heterog
        } else {
            CoverClass::Good
        }
    }
}

/// Statistics on how many (distinct, non-singleton) barcodes have a
/// read starting at a certain position. When multiple reads for a
/// single barcode have the same read start position, this is counted
/// just once.
pub struct ReadStartStats {
    starts: Vec<u32>,
    nom_cover: Vec<u32>,
}

impl ReadStartStats {
    /// Creates a new read start tabulator for a sequence of a
    /// specified length.
    pub fn new(len: usize) -> Self {
        ReadStartStats {
            starts: vec![0; len],
            nom_cover: vec![0; len],
        }
    }

    /// Increment read start counts for all distinct starts in the
    /// collection of `rust_htslib::bam::Record` objects.
    pub fn add_read_group<'a, I>(&mut self, record_iter: I) -> (u32, u32)
    where
        I: Iterator<Item = &'a bam::Record>,
    {
        let mut seen = HashSet::new();
        let mut nomcov = vec![0; self.nom_cover.len()];

        for r in record_iter {
            let pos = r.pos() as usize;
            if !seen.contains(&pos) && (pos < self.starts.len()) {
                self.starts[pos] += 1;
                for i in pos..min(pos + r.seq().len() - 1, nomcov.len() - 1) {
                    nomcov[i] = 1;
                }
                seen.insert(pos);
            }
        }

        let mut nnomcov = 0;
        for i in 0..(nomcov.len() - 1) {
            nnomcov += nomcov[i];
            self.nom_cover[i] += nomcov[i];
        }

        (seen.len() as u32, nnomcov)
    }

    /// Returns a tab-delimited table of read starts counted at each
    /// position.
    pub fn table(&self) -> String {
        let mut buf = String::new();

        for (pos, nstart) in self.starts.iter().enumerate() {
            buf.push_str(&pos.to_string());
            buf.push('\t');
            buf.push_str(&nstart.to_string());
            buf.push('\t');
            buf.push_str(&self.nom_cover[pos].to_string());
            buf.push('\n');
        }

        buf
    }
}

pub struct CoverStats(OffsetVector<CoverAt>);

impl CoverStats {
    pub fn new(start: usize, len: usize, max: usize) -> Self {
        CoverStats(OffsetVector::new(
            start,
            repeat(CoverAt::new(max)).take(len).collect(),
        ))
    }

    pub fn add_coverage<'a, I>(&mut self, pos_iter: I, het_fract: f64)
    where
        I: Iterator<Item = (usize, &'a AlnPos)>,
    {
        for (pos, ap) in pos_iter {
            if let Some(at) = self.0.get_mut(pos) {
                at.add_coverage(ap, het_fract);
            }
        }
    }

    pub fn table(&self) -> String {
        let mut buf = String::new();

        for (pos, at) in self.0.pos_iter() {
            buf.push_str(&pos.to_string());

            for cvg in at.cover() {
                buf.push('\t');
                buf.push_str(&cvg.to_string());
            }

            buf.push('\n');
        }

        buf
    }
}

#[derive(Debug, Clone)]
struct CoverAt {
    cover: Vec<u32>,
}

impl CoverAt {
    fn new(max: usize) -> Self {
        CoverAt {
            cover: vec![0; max + 1],
        }
    }

    fn _max(&self) -> usize {
        self.cover.len()
    }

    fn cover(&self) -> &[u32] {
        &self.cover
    }

    fn add_coverage(&mut self, ap: &AlnPos, het_fract: f64) {
        let cvg = if AlnPosCons::new(ap, het_fract).is_good() {
            ap.nttl()
        } else {
            0
        };
        let idx = if cvg >= self.cover.len() {
            self.cover.len() - 1
        } else {
            cvg
        };
        self.cover[idx] += 1;
    }
}
