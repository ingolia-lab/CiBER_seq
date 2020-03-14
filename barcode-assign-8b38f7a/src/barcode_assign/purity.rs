use rust_htslib::bam;

use assign::ReadAssign;

#[derive(Debug, Clone, Hash, PartialEq, Eq)]
pub struct GrnaPurity {
    assign_counts: Vec<(ReadAssign, usize)>,
}

impl GrnaPurity {
    pub fn new<'a, I>(r_iter: I) -> Result<Self, failure::Error>
    where
        I: Iterator<Item = &'a bam::Record>,
    {
        let mut cts = Vec::new();

        for r in r_iter {
            let asn = ReadAssign::new(r)?;

            let mut extant = false;
            for &mut (ref mut a, ref mut n) in cts.iter_mut() {
                if asn == *a {
                    *n += 1;
                    extant = true;
                }
            }
            if !extant {
                cts.push((asn, 1));
            }
        }

        cts.sort_by_key(|&(ref _a, ref n)| -(*n as isize));

        Ok(GrnaPurity { assign_counts: cts })
    }

    pub fn primary_assign(&self) -> ReadAssign {
        self.assign_counts
            .first()
            .map_or(ReadAssign::NoMatch, |&(ref a, ref _n)| a.clone())
    }

    /// Returns the fraction of reads in the group that match exactly
    /// the primary assignment for the group.
    pub fn assign_purity(&self) -> f64 {
        self.assign_counts
            .first()
            .map_or(0.0, |&(ref _a, ref nprim)| {
                (*nprim as f64) / (self.n_total() as f64)
            })
    }

    /// Returns the fraction of reads in the group that have the same
    /// target as the primary assignment for the group.
    pub fn target_purity(&self) -> f64 {
        let ttl = self.n_total();
        if ttl > 0 {
            ((self.n_primary() + self.n_primary_target()) as f64) / (ttl as f64)
        } else {
            0.0
        }
    }

    /// Returns the total number of alignments in the group
    pub fn n_total(&self) -> usize {
        self.assign_counts.iter().map(|&(ref _a, ref n)| *n).sum()
    }

    /// Returns the number of alignments for the primary assignment
    pub fn n_primary(&self) -> usize {
        self.assign_counts.first().map_or(0, |&(ref _a, ref n)| *n)
    }

    /// Returns the number of alignments with an assignment to the
    /// same target sequence, but different specific alignments
    pub fn n_primary_target(&self) -> usize {
        if let Some(primary_target) = self.assign_counts.first().map(|&(ref a, ref _n)| a.tid()) {
            self.assign_counts
                .iter()
                .skip(1)
                .filter(|&(ref a, ref _n)| a.tid() == primary_target)
                .map(|&(ref _a, ref n)| *n)
                .sum()
        } else {
            0
        }
    }

    /// Returns the number of unaligned, non-primary assignments. This
    /// will be the number of unaligned reads in the group -- unless
    /// the primary assignment is not aligned, in which case it will
    /// be zero.
    pub fn n_no_match(&self) -> usize {
        self.assign_counts
            .iter()
            .skip(1)
            .filter(|&(ref a, ref _n)| a.is_no_match())
            .map(|&(ref _a, ref n)| *n)
            .sum()
    }

    pub fn n_other_match(&self) -> usize {
        self.assign_counts
            .iter()
            .skip(1)
            .filter(|&(ref a, ref _n)| !a.is_no_match())
            .map(|&(ref _a, ref n)| *n)
            .sum()
    }

    /// Returns the number of alignments with an assignment to a
    /// different target sequence than the primary assignment for the
    /// group. This number doesn't include unaligned reads
    pub fn n_other_target_match(&self) -> usize {
        if let Some(primary_target) = self.assign_counts.first().map(|&(ref a, ref _n)| a.tid()) {
            self.assign_counts
                .iter()
                .skip(1)
                .filter(|&(ref a, ref _n)| a.tid() != primary_target && !a.is_no_match())
                .map(|&(ref _a, ref n)| *n)
                .sum()
        } else {
            0
        }
    }

    pub fn header() -> String {
        "barcode\ttotal\tprimary\tprimary_target\tother_match\tno_match".to_string()
    }

    pub fn line(&self, bc_str: &str) -> String {
        format!(
            "{}\t{}\t{}\t{}\t{}\t{}",
            bc_str,
            self.n_total(),
            self.n_primary(),
            self.n_primary_target(),
            self.n_other_target_match(),
            self.n_no_match()
        )
    }
}
