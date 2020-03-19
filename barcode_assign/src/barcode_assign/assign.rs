use std::fmt;
use std::str;

use rust_htslib::bam;
use rust_htslib::bam::record::{Cigar, CigarString};

#[derive(Debug, Clone, Hash, PartialEq, Eq)]
pub struct AssignPos {
    tid: u32,
    pos: i64,
    is_reverse: bool,
}

impl AssignPos {
    pub fn new(r: &bam::Record) -> Option<Self> {
        if r.tid() < 0 {
            None
        } else {
            Some(AssignPos {
                tid: r.tid() as u32,
                pos: r.pos(),
                is_reverse: r.is_reverse(),
            })
        }
    }

    #[allow(dead_code)]
    pub fn tid(&self) -> u32 {
        self.tid
    }

    pub fn target(&self, targets: &[&str]) -> String {
        targets.get(self.tid as usize).unwrap_or(&"???").to_string()
    }

    pub fn pos(&self) -> i64 {
        self.pos
    }
    pub fn is_reverse(&self) -> bool {
        self.is_reverse
    }
    pub fn strand_chr(&self) -> char {
        if self.is_reverse {
            '-'
        } else {
            '+'
        }
    }
}

impl fmt::Display for AssignPos {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(
            f,
            "{}:{}{}",
            self.tid,
            self.pos,
            if self.is_reverse { "-" } else { "+" }
        )
    }
}

#[derive(Debug, Clone, Hash, PartialEq, Eq)]
pub enum ReadAssign {
    NoMatch,
    Match(AssignMatch),
}

impl ReadAssign {
    pub fn new(r: &bam::Record) -> Result<Self, failure::Error> {
        if let Some(assign_pos) = AssignPos::new(r) {
            let md_aux = r.aux(b"MD").ok_or_else(|| failure::err_msg("No MD tag"))?;
            let md: Vec<u8> = match md_aux {
                bam::record::Aux::String(md) => Ok(md.to_vec()),
                _ => Err(failure::err_msg("MD tag not a string")),
            }?;

            let cigar = (*r.cigar()).clone();
            let len = ReadAssign::cigar_string_len(&cigar);

            let assign_match = AssignMatch {
                assign_pos: assign_pos,
                len: len,
                cigar: cigar,
                md: md,
            };
            Ok(ReadAssign::Match(assign_match))
        } else {
            Ok(ReadAssign::NoMatch)
        }
    }

    fn cigar_string_len(cigar_str: &CigarString) -> u32 {
        cigar_str.iter().map(ReadAssign::cigar_len).sum()
    }

    fn cigar_len(cigar: &Cigar) -> u32 {
        match cigar {
            Cigar::Match(ref l) => *l,
            Cigar::Ins(_) => 0,
            Cigar::Del(ref l) => *l,
            Cigar::RefSkip(ref l) => *l,
            Cigar::SoftClip(_) => 0,
            Cigar::HardClip(_) => 0,
            Cigar::Pad(_) => 0,
            Cigar::Equal(ref l) => *l,
            Cigar::Diff(ref l) => *l,
        }
    }

    pub fn assign_match(&self) -> Option<&AssignMatch> {
        match self {
            ReadAssign::NoMatch => None,
            ReadAssign::Match(m) => Some(m),
        }
    }

    pub fn is_no_match(&self) -> bool {
        match self {
            ReadAssign::NoMatch => true,
            _ => false,
        }
    }

    pub fn tid(&self) -> Option<u32> {
        match self {
            ReadAssign::Match(am) => Some(am.tid()),
            ReadAssign::NoMatch => None,
        }
    }

    pub fn assign_pos(&self) -> Option<&AssignPos> {
        self.assign_match().map(AssignMatch::assign_pos)
    }
}

#[derive(Debug, Clone, Hash, PartialEq, Eq)]
pub struct AssignMatch {
    assign_pos: AssignPos,
    len: u32,
    cigar: CigarString,
    md: Vec<u8>,
}

impl AssignMatch {
    #[allow(dead_code)]
    pub fn tid(&self) -> u32 {
        self.assign_pos.tid()
    }

    pub fn target(&self, targets: &[&str]) -> String {
        self.assign_pos.target(targets)
    }

    pub fn pos(&self) -> i64 {
        self.assign_pos.pos()
    }
    pub fn is_reverse(&self) -> bool {
        self.assign_pos.is_reverse()
    }
    pub fn len(&self) -> u32 {
        self.len
    }

    pub fn assign_pos(&self) -> &AssignPos {
        &self.assign_pos
    }

    #[allow(dead_code)]
    pub fn cigar(&self) -> CigarString {
        self.cigar.clone()
    }
    pub fn cigar_string(&self) -> String {
        self.cigar.to_string()
    }
    pub fn md(&self) -> &[u8] {
        self.md.as_slice()
    }

    pub fn is_cigar_perfect(&self) -> bool {
        let CigarString(ref v) = self.cigar;
        (v.len() == 1)
            && (match v[0] {
                Cigar::Match(_) => true,
                _ => false,
            })
    }

    pub fn is_md_perfect(&self) -> bool {
        self.md.iter().all(u8::is_ascii_digit)
    }

    pub fn header() -> String {
        "barcode\ttid\tpos\tcigar\tmd".to_string()
    }

    pub fn line(&self, bc_str: &str, targets: &[&str]) -> Result<String, failure::Error> {
        Ok(format!(
            "{}\t{}\t{}\t{:?}\t{}",
            bc_str,
            self.target(&targets),
            self.pos(),
            self.cigar_string(),
            str::from_utf8(self.md())?
        ))
    }

    pub fn field(&self) -> String {
        format!(
            "{}:{}",
            self.cigar_string(),
            String::from_utf8_lossy(self.md())
        )
    }
}
