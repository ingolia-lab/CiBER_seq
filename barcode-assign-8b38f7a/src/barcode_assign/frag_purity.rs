use std::collections::HashMap;

use rust_htslib::bam;

use assign::AssignPos;

#[derive(Debug, Clone, Hash, PartialEq, Eq)]
pub struct FragPurity {
    primary_pos: Option<AssignPos>,
    n_primary: u64,
    n_other_aligned: u64,
    n_other_unaligned: u64,
}

impl FragPurity {
    pub fn new<'a, I>(r_iter: I) -> Self
    where
        I: Iterator<Item = &'a bam::Record>,
    {
        let poses: Vec<Option<AssignPos>> = r_iter.map(AssignPos::new).collect();

        let mut count_map = HashMap::new();
        for ap in poses.iter() {
            let ct = count_map.entry(ap).or_insert(0);
            *ct += 1;
        }
        let counts: Vec<(Option<AssignPos>, u64)> = count_map
            .into_iter()
            .map(|(ap, n)| (ap.clone(), n))
            .collect();

        // Bias towards aligned for equal counts of aligned, unaligned
        let (primary, _) = counts
            .iter()
            .max_by_key(|(ap, n)| (n, ap.is_some()))
            .unwrap();

        let mut frag_purity = FragPurity {
            primary_pos: primary.clone(),
            n_primary: 0,
            n_other_aligned: 0,
            n_other_unaligned: 0,
        };

        for (apo, n) in counts.iter() {
            if apo == primary {
                frag_purity.n_primary += n;
            } else if apo.is_some() {
                frag_purity.n_other_aligned += n;
            } else {
                frag_purity.n_other_unaligned += n;
            }
        }

        frag_purity
    }

    pub fn primary_pos(&self) -> &Option<AssignPos> {
        &self.primary_pos
    }

    pub fn n_primary(&self) -> u64 {
        self.n_primary
    }

    pub fn n_other_aligned(&self) -> u64 {
        self.n_other_aligned
    }

    pub fn n_other_unaligned(&self) -> u64 {
        self.n_other_unaligned
    }

    pub fn n_total(&self) -> u64 {
        self.n_primary() + self.n_other_aligned() + self.n_other_unaligned()
    }

    pub fn is_primary_unaligned(&self) -> bool {
        self.primary_pos.is_none()
    }

    pub fn read_purity(&self) -> f64 {
        (self.n_primary as f64) / (self.n_total() as f64)
    }

    pub fn align_purity(&self) -> f64 {
        if self.primary_pos.is_none() {
            0.0
        } else {
            (self.n_primary as f64) / ((self.n_primary + self.n_other_aligned) as f64)
        }
    }

    pub fn header() -> String {
        "primary\tn_total\tn_primary\tn_other_aligned\tn_other_unaligned".to_string()
    }

    pub fn primary_str(&self, targets: &[&str]) -> String {
        match self.primary_pos {
            Some(ref ap) => format!(
                "{}:{}{}",
                ap.target(targets),
                ap.pos(),
                if ap.is_reverse() { "-" } else { "+" }
            ),
            None => "N/A".to_string(),
        }
    }

    pub fn line(&self, targets: &[&str]) -> String {
        format!(
            "{}\t{}\t{}\t{}\t{}",
            self.primary_str(targets),
            self.n_total(),
            self.n_primary,
            self.n_other_aligned,
            self.n_other_unaligned
        )
    }
}

#[allow(dead_code)]
fn align_pos_counts<'a, I>(r_iter: I) -> Vec<(Option<AssignPos>, usize)>
where
    I: Iterator<Item = &'a bam::Record>,
{
    let mut count_map = HashMap::new();
    for r in r_iter {
        let ap = AssignPos::new(r);
        let ct = count_map.entry(ap).or_insert(0);
        *ct += 1;
    }
    count_map.into_iter().collect()
}

// 1. Identify primary assignment position
// 2. Fraction of primary, other aligned, and unaligned
// 3. Extract primary assignment rea
