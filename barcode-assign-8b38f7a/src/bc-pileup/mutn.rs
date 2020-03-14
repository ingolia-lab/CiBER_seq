use std::collections::HashMap;
use std::fmt::{Display, Formatter};
use std::fs::File;
use std::io::Write;
use std::path::Path;
use std::rc::Rc;
use std::str;

use aln_pos::AlnCons;
use trl::STD_CODONS;

#[derive(Debug, Clone, Hash, PartialEq, Eq, PartialOrd, Ord)]
pub struct NtMutn {
    pos: usize,
    refnt: u8,
    mutseq: Vec<u8>,
}

#[allow(dead_code)]
static NTS: [u8; 4] = [b'A', b'C', b'G', b'T'];

impl NtMutn {
    pub fn new(pos: usize, refnt: u8, mutseq: Vec<u8>) -> Self {
        NtMutn {
            pos: pos,
            refnt: refnt,
            mutseq: mutseq,
        }
    }

    pub fn pos(&self) -> usize {
        self.pos
    }
    pub fn mutseq(&self) -> &[u8] {
        &self.mutseq
    }
    pub fn refnt<T: From<u8>>(&self) -> T {
        T::from(self.refnt)
    }

    #[allow(dead_code)]
    pub fn all_substs(start: usize, len: usize, refseq: &[u8]) -> Vec<NtMutn> {
        let mut substs = Vec::new();
        for pos in start..(start + len) {
            if let Some(refnt) = refseq.get(pos) {
                for nt in NTS.iter().filter(|&nt| nt != refnt) {
                    substs.push(NtMutn::new(pos, *refnt, vec![*nt]));
                }
            }
        }
        substs
    }

    pub fn tsv(&self) -> String {
        format!(
            "{}\t{}\t{}",
            self.pos,
            self.refnt::<char>(),
            str::from_utf8(&self.mutseq).unwrap_or("???")
        )
    }
}

#[derive(Debug, Clone, Hash, PartialEq, Eq, PartialOrd, Ord)]
pub struct PeptMutn {
    pos_aa: usize,
    ref_aa: u8,
    cons_aa: u8,
}

impl PeptMutn {
    pub fn new(pos_aa: usize, ref_aa: u8, cons_aa: u8) -> Self {
        PeptMutn {
            pos_aa: pos_aa,
            ref_aa: ref_aa,
            cons_aa: cons_aa,
        }
    }

    pub fn tsv(&self) -> String {
        format!(
            "{}\t{}\t{}",
            self.pos_aa, self.ref_aa as char, self.cons_aa as char
        )
    }
}

impl Display for PeptMutn {
    fn fmt(&self, f: &mut Formatter) -> ::std::fmt::Result {
        write!(
            f,
            "{}{}{}",
            self.ref_aa as char, self.pos_aa, self.cons_aa as char
        )
    }
}

#[derive(Debug, Clone)]
pub struct MutnBarcodes {
    mut_bc: HashMap<NtMutn, Vec<Rc<String>>>,
    bc_mut: HashMap<Rc<String>, Vec<NtMutn>>,
}

impl MutnBarcodes {
    pub fn new() -> Self {
        MutnBarcodes {
            mut_bc: HashMap::new(),
            bc_mut: HashMap::new(),
        }
    }

    pub fn insert(&mut self, mutn: NtMutn, barcode: Rc<String>) {
        let mutn_vec = self
            .mut_bc
            .entry(mutn.clone())
            .or_insert_with(|| Vec::new());
        mutn_vec.push(barcode.clone());
        let bc_vec = self.bc_mut.entry(barcode).or_insert_with(|| Vec::new());
        bc_vec.push(mutn);
    }

    #[allow(dead_code)]
    pub fn count_table(&self) -> String {
        let mut buf = String::new();

        let mut mutn_barcodes: Vec<(&NtMutn, &Vec<Rc<String>>)> = self.mut_bc.iter().collect();
        mutn_barcodes.sort_by_key(|&(&ref mutn, _barcodes)| mutn);
        for (&ref mutn, &ref barcodes) in mutn_barcodes.into_iter() {
            buf.push_str(&mutn.tsv());
            buf.push_str(&format!("\t{}\n", barcodes.len()));
        }

        buf
    }

    #[allow(dead_code)]
    pub fn barcode_table(&self) -> String {
        let mut buf = String::new();

        let mut mutn_barcodes: Vec<(&NtMutn, &Vec<Rc<String>>)> = self.mut_bc.iter().collect();
        mutn_barcodes.sort_by_key(|&(&ref mutn, _barcodes)| mutn);
        for (&ref mutn, &ref barcodes) in mutn_barcodes.into_iter() {
            buf.push_str(&mutn.tsv());
            for barcode in barcodes {
                buf.push('\t');
                buf.push_str(barcode.as_str());
            }
            buf.push('\n');
        }

        buf
    }

    fn is_single_mutn(&self, barcode: &Rc<String>) -> bool {
        self.bc_mut
            .get(barcode)
            .map_or(false, |&ref muts| muts.len() == 1)
    }

    #[allow(dead_code)]
    pub fn mutn_table<'a, I>(&self, mutns: I) -> String
    where
        I: Iterator<Item = &'a NtMutn>,
    {
        let mut buf = String::new();

        for mutn in mutns {
            let barcodes = self.mut_bc.get(mutn);
            let n_total = barcodes.map_or(0, Vec::len);
            let n_single = barcodes.map_or(0, |&ref bcs| {
                bcs.iter().filter(|&ref bc| self.is_single_mutn(bc)).count()
            });
            //            let n_single = barcodes.map_or(0, |&bcs| bcs.iter().filter(&is_only).count());

            buf.push_str(&mutn.tsv());
            buf.push_str(&format!("\t{}\t{}\n", n_total, n_single));
        }

        buf
    }
}

pub struct MutnAnalysis {
    mutn_out: File,
    pept_mutn_out: File,
    exon_start: usize,
    exon_upstream: Vec<u8>,
    mutn_barcodes: MutnBarcodes,
}

impl MutnAnalysis {
    pub fn new(
        outdir: &Path,
        prefix: &str,
        exon_start: usize,
        exon_upstream: &[u8],
    ) -> Result<Self, failure::Error> {
        let mut mutn_out = File::create(
            outdir
                .to_path_buf()
                .join(prefix.to_string() + "barcode-mutations.txt"),
        )?;
        let mut pept_mutn_out = File::create(
            outdir
                .to_path_buf()
                .join(prefix.to_string() + "barcode-peptide-mutations.txt"),
        )?;
        let mutn_barcodes = MutnBarcodes::new();

        write!(mutn_out, "Barcode\tPos\tRef\tRead\n")?;
        write!(pept_mutn_out, "Barcode\tPos\tRef\tRead\n")?;

        Ok(MutnAnalysis {
            mutn_out: mutn_out,
            pept_mutn_out: pept_mutn_out,
            exon_start: exon_start,
            exon_upstream: exon_upstream.to_vec(),
            mutn_barcodes: mutn_barcodes,
        })
    }

    pub fn analyze_aln_cons(
        self: &mut Self,
        bc_str: &str,
        aln_cons: AlnCons,
    ) -> Result<(), failure::Error> {
        let mut_posn = aln_cons
            .pos_iter()
            .filter(|&(_pos, apc)| !apc.is_wildtype() && !apc.is_uncovered());
        let bc_rc = Rc::new(bc_str.to_owned());
        for (pos, apc) in mut_posn {
            let mutn = NtMutn::new(pos, apc.ref_nt(), apc.mutseq());
            write!(
                self.mutn_out,
                "{}\t{}\t{}\t{}\n",
                bc_str,
                mutn.pos(),
                mutn.refnt::<char>(),
                str::from_utf8(mutn.mutseq()).unwrap_or("???")
            )?;
            self.mutn_barcodes.insert(mutn, bc_rc.clone());
        }

        let mut cds_ref = self.exon_upstream.clone();
        let mut cds_cons = self.exon_upstream.clone();

        let mut frameshift_nt = None;
        for (pos, apc) in aln_cons
            .pos_iter()
            .filter(|&(pos, _apc)| pos >= self.exon_start)
        {
            apc.push_cons_seq(&mut cds_cons);
            cds_ref.push(apc.ref_nt());
            if apc.is_frameshift() && frameshift_nt.is_none() {
                frameshift_nt = Some(pos);
            }
        }

        let frameshift_codon =
            frameshift_nt.map(|nt| (nt + self.exon_upstream.len() - self.exon_start) / 3);

        let mut _nonsense_codon = None;

        let pept_ref = STD_CODONS.trl(&cds_ref);
        let pept_cons = STD_CODONS.trl(&cds_cons);
        for (pos, (ref_aa, cons_aa)) in pept_ref.iter().zip(pept_cons.iter()).enumerate() {
            if frameshift_codon.map_or(false, |c| pos >= c) {
                let pept_mutn = PeptMutn::new(pos, *ref_aa, b'!');
                write!(self.pept_mutn_out, "{}\t{}\n", bc_str, pept_mutn.tsv())?;
                break;
            }

            if cons_aa != ref_aa {
                let pept_mutn = PeptMutn::new(pos, *ref_aa, *cons_aa);
                if *cons_aa == b'*' {
                    _nonsense_codon = Some(pos);
                    write!(self.pept_mutn_out, "{}\t{}\n", bc_str, pept_mutn.tsv())?;
                    break;
                } else {
                    write!(self.pept_mutn_out, "{}\t{}\n", bc_str, pept_mutn.tsv())?;
                }
            }
        }

        Ok(())
    }
}
