use rust_htslib::bam;

#[derive(Debug, Clone, Hash, PartialEq, Eq)]
pub struct Depth {
    n_good: usize,
    n_lowqual: usize,
}

impl Depth {
    pub fn new(n_good: usize, n_lowqual: usize) -> Self {
        Depth {
            n_good: n_good,
            n_lowqual: n_lowqual,
        }
    }

    pub fn header() -> String {
        "barcode\tn_good\tn_lowqual".to_string()
    }

    pub fn line(&self, bc_str: &str) -> String {
        format!("{}\t{}\t{}", bc_str, self.n_good, self.n_lowqual)
    }

    pub fn n_good(&self) -> usize {
        self.n_good
    }

    #[allow(dead_code)]
    pub fn n_lowqual(&self) -> usize {
        self.n_lowqual
    }

    #[allow(dead_code)]
    pub fn n_total(&self) -> usize {
        self.n_good + self.n_lowqual
    }
}

pub fn filter_by_quality(rec_all: Vec<bam::Record>, min_qual: u8) -> (Vec<bam::Record>, Depth) {
    let n_all = rec_all.len();

    let rec_wanted = rec_all
        .into_iter()
        .filter(|r| (minimum_qual(r) >= min_qual))
        .collect::<Vec<bam::Record>>();

    let n_wanted = rec_wanted.len();
    let n_lowqual = n_all - n_wanted;

    (rec_wanted, Depth::new(n_wanted, n_lowqual))
}

fn minimum_qual(r: &bam::Record) -> u8 {
    *r.qual().iter().min().unwrap_or(&0)
}

#[allow(dead_code)]
fn median_qual(r: &bam::Record) -> u8 {
    let mut quals = r.qual().to_vec();
    quals.sort();
    quals.get(quals.len() / 2).map_or(0, |q| *q)
}
