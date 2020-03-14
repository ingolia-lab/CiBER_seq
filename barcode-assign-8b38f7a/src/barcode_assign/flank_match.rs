use std::cmp::*;

use bio::alphabets::dna;
use bio::pattern_matching::myers::Myers;

pub struct LibSpec {
    name: String,
    frag_matcher: FlankMatchSpec,
    barcode_matcher: FlankMatchSpec,
    barcode_rev: bool,
}

impl LibSpec {
    pub fn new(
        name: &str,
        frag_matcher: FlankMatchSpec,
        barcode_matcher: FlankMatchSpec,
        barcode_rev: bool,
    ) -> Self {
        LibSpec {
            name: name.to_string(),
            frag_matcher: frag_matcher,
            barcode_matcher: barcode_matcher,
            barcode_rev: barcode_rev,
        }
    }

    pub fn name(&self) -> &str {
        &self.name
    }

    pub fn best_match<'a>(&mut self, query: &'a [u8]) -> LibMatchOut<'a> {
        LibMatchOut {
            frag: self.frag_matcher.best_match(query),
            barcode: self.barcode_matcher.best_match(query),
            barcode_rev: self.barcode_rev,
        }
    }
}

#[derive(Debug, Clone, Hash, PartialEq, Eq)]
pub struct LibMatchOut<'a> {
    frag: FlankMatchOut<'a>,
    barcode: FlankMatchOut<'a>,
    barcode_rev: bool,
}

impl<'a> LibMatchOut<'a> {
    pub fn frag_match(&self) -> &FlankMatchOut<'a> {
        &self.frag
    }
    pub fn barcode_match(&self) -> &FlankMatchOut<'a> {
        &self.barcode
    }
    pub fn barcode_rev(&self) -> bool {
        self.barcode_rev
    }

    pub fn lib_match(&self) -> Option<LibMatch<'a>> {
        self.frag.flank_match().and_then(|frag| {
            self.barcode.flank_match().map(|barcode| LibMatch {
                frag: frag,
                barcode: barcode,
                barcode_rev: self.barcode_rev,
            })
        })
    }
}

#[derive(Debug, Clone, Hash, PartialEq, Eq)]
pub struct LibMatch<'a> {
    frag: FlankMatch<'a>,
    barcode: FlankMatch<'a>,
    barcode_rev: bool,
}

impl<'a> LibMatch<'a> {
    pub fn frag_match(&self) -> &FlankMatch<'a> {
        &self.frag
    }
    pub fn barcode_match(&self) -> &FlankMatch<'a> {
        &self.barcode
    }
    pub fn barcode_rev(&self) -> bool {
        self.barcode_rev
    }

    pub fn between(&self) -> isize {
        if self.frag.insert_start() > self.barcode.insert_end() {
            self.frag.insert_start() as isize - self.barcode.insert_end() as isize
        } else {
            self.barcode.insert_start() as isize - self.frag.insert_end() as isize
        }
    }

    pub fn barcode_actual(&self) -> String {
        if self.barcode_rev {
            String::from_utf8_lossy(&dna::revcomp(self.barcode.insert_seq())).to_string()
        } else {
            String::from_utf8_lossy(self.barcode.insert_seq()).to_string()
        }
    }
}

pub struct FlankMatchSpec {
    before_myers: Myers<u64>,
    after_myers: Myers<u64>,
    #[allow(dead_code)]
    before: Vec<u8>,
    #[allow(dead_code)]
    after: Vec<u8>,
    max_errors: u8,
}

impl FlankMatchSpec {
    pub fn new(before: &[u8], after: &[u8], max_errors: u8) -> Self {
        FlankMatchSpec {
            before_myers: Myers::<u64>::new(before),
            after_myers: Myers::<u64>::new(after),
            before: before.to_vec(),
            after: after.to_vec(),
            max_errors: max_errors,
        }
    }

    pub fn best_match<'a>(&mut self, query: &'a [u8]) -> FlankMatchOut<'a> {
        // N.B. end coordinate is not included in match
        let best_before = self
            .before_myers
            .find_all(query, self.max_errors)
            .by_ref()
            .min_by_key(|&(_, _, dist)| dist);
        let best_after = self
            .after_myers
            .find_all(query, self.max_errors)
            .by_ref()
            .min_by_key(|&(_, _, dist)| dist);

        FlankMatchOut {
            before: best_before,
            after: best_after,
            query: query,
        }
    }
}

/// Attempted alignment between flanking constant sequences and a
/// target query. Either of the flanking sequence matches may fail.
#[derive(Debug, Clone, Hash, PartialEq, Eq)]
pub struct FlankMatchOut<'a> {
    before: Option<(usize, usize, u8)>,
    after: Option<(usize, usize, u8)>,
    query: &'a [u8],
}

const DESC_LEN: usize = 10;

impl<'a> FlankMatchOut<'a> {
    /// Returns the successful match, or `None` if either the before
    /// or after match failed.
    pub fn flank_match(&self) -> Option<FlankMatch<'a>> {
        if let Some(before) = self.before {
            if let Some(after) = self.after {
                if before.1 <= after.0 {
                    Some(FlankMatch {
                        before: before,
                        after: after,
                        query: self.query,
                    })
                } else {
                    None
                }
            } else {
                None
            }
        } else {
            None
        }
    }

    pub fn before_match_desc(&self) -> String {
        if let Some((_start, end, _score)) = self.before {
            String::from_utf8_lossy(
                &self.query[max(DESC_LEN, end) - DESC_LEN..end].to_ascii_uppercase(),
            )
            .to_string()
                + &String::from_utf8_lossy(
                    &self.query[end..min(end + DESC_LEN, self.query.len())].to_ascii_lowercase(),
                )
                .to_string()
        } else {
            "N/A".to_string()
        }
    }

    pub fn after_match_desc(&self) -> String {
        if let Some((start, _end, _score)) = self.after {
            String::from_utf8_lossy(
                &self.query[max(DESC_LEN, start) - DESC_LEN..start].to_ascii_lowercase(),
            )
            .to_string()
                + &String::from_utf8_lossy(
                    &self.query[start..min(start + DESC_LEN, self.query.len())]
                        .to_ascii_uppercase(),
                )
                .to_string()
        } else {
            "N/A".to_string()
        }
    }
}

/// Successful alignment between flanking constant sequences and a
/// target query.
#[derive(Debug, Clone, Hash, PartialEq, Eq)]
pub struct FlankMatch<'a> {
    before: (usize, usize, u8),
    after: (usize, usize, u8),
    query: &'a [u8],
}

impl<'a> FlankMatch<'a> {
    /// Returns the total score (edit distance) of the alignments
    /// between the before and after flanking sequences and the query.
    pub fn score(&self) -> u8 {
        self.before.2 + self.after.2
    }

    /// Returns the starting position of the insert sequence in the
    /// query.
    pub fn insert_start(&self) -> usize {
        self.before.1
    }

    /// Returns the ending position of the insert sequence in the
    /// query, *non*-inclusive.
    pub fn insert_end(&self) -> usize {
        self.after.0
    }

    /// Returns the insert sequence.
    pub fn insert_seq(&self) -> &[u8] {
        &self.query[self.before.1..self.after.0]
    }

    /// Returns the query sequence that matched the constant sequence
    /// before the insert.
    pub fn before_seq(&self) -> &[u8] {
        &self.query[self.before.0..self.before.1]
    }

    /// Returns the query sequence that matched the constant sequence
    /// after the insert.
    pub fn after_seq(&self) -> &[u8] {
        &self.query[self.after.0..self.after.1]
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn perfect_match() {
        // 012345678901234567890123456789
        // CTAACGTACGTTTAATTAACAGTCAGTAAG

        let upstream = b"CTA";
        let before = b"ACGTACGT";
        let insert = b"TTAATTAA";
        let after = b"CAGTCAGT";
        let downstream = b"AAG";

        let query_perfect = build_query(upstream, before, insert, after, downstream);
        assert_eq!(query_perfect, b"CTAACGTACGTTTAATTAACAGTCAGTAAG");

        let insert_start = upstream.len() + before.len();
        let insert_end = insert_start + insert.len();

        let mut match_spec_a = FlankMatchSpec::new(before, after, 0);
        let match_a = match_spec_a
            .best_match(&query_perfect)
            .flank_match()
            .unwrap();
        assert_eq!(match_a.score(), 0);
        assert_eq!(match_a.insert_start(), insert_start);
        assert_eq!(match_a.insert_end(), insert_end);
        assert_eq!(
            match_a.insert_seq(),
            &query_perfect[insert_start..insert_end]
        );

        let before_mut = b"ACGTTCGT";
        let query_before_mut = build_query(upstream, before_mut, insert, after, downstream);
        let match_before_mut = match_spec_a.best_match(&query_before_mut).flank_match();
        assert_eq!(match_before_mut, None);
        let mut match_spec_b = FlankMatchSpec::new(before_mut, after, 0);
        let match_b = match_spec_b.best_match(&query_perfect).flank_match();
        assert_eq!(match_b, None);

        let after_mut = b"CACTCAGT";
        let query_after_mut = build_query(upstream, before, insert, after_mut, downstream);
        let match_after_mut = match_spec_a.best_match(&query_after_mut).flank_match();
        assert_eq!(match_after_mut, None);
        let mut match_spec_c = FlankMatchSpec::new(before, after_mut, 0);
        let match_c = match_spec_c.best_match(&query_perfect).flank_match();
        assert_eq!(match_c, None);
    }

    #[test]
    fn imperfect_match() {
        //             0123456789012345678901234567
        assert_match(
            b"AGATCTCGCGAGAATTAACGTGAGTGCC",
            b"CTCGC",
            b"CGTGA",
            0,
            9,
            18,
        );
        assert_match(
            b"AGATCTCGCGAGAATTAACGTGAGTGCC",
            b"CTAGC",
            b"CGTGA",
            1,
            9,
            18,
        );
        assert_match(
            b"AGATCTCGCGAGAATTAACGTGAGTGCC",
            b"CTCGC",
            b"CCTGA",
            1,
            9,
            18,
        );

        assert_match(
            b"AGATCTCGCGAGAATTAACGTGAGTGCC",
            b"CTTCGC",
            b"CGTGA",
            1,
            9,
            18,
        );
        assert_match(
            b"AGATCTCGCGAGAATTAACGTGAGTGCC",
            b"CTCGC",
            b"CGTGGA",
            1,
            9,
            18,
        );

        assert_no_match(b"AGATCTCGCGAGAATTAACGTGAGTGCC", b"CTAGC", b"CGTGA", 0);
        assert_no_match(b"AGATCTCGCGAGAATTAACGTGAGTGCC", b"CTCGC", b"CCTGA", 0);
        assert_no_match(b"AGATCTCGCGAGAATTAACGTGAGTGCC", b"CTTCGC", b"CGTGA", 0);
        assert_no_match(b"AGATCTCGCGAGAATTAACGTGAGTGCC", b"CTCGC", b"CGTGGA", 0);
    }

    fn assert_match(
        query: &[u8],
        before: &[u8],
        after: &[u8],
        max_errors: u8,
        insert_start: usize,
        insert_end: usize,
    ) {
        let mut match_spec = FlankMatchSpec::new(before, after, max_errors);
        let match_out = match_spec.best_match(query).flank_match().unwrap();
        assert_eq!(match_out.insert_start(), insert_start);
        assert_eq!(match_out.insert_end(), insert_end);
    }

    fn assert_no_match(query: &[u8], before: &[u8], after: &[u8], max_errors: u8) {
        let mut match_spec = FlankMatchSpec::new(before, after, max_errors);
        let match_out = match_spec.best_match(query).flank_match();
        assert_eq!(match_out, None);
    }

    fn build_query(
        upstream: &[u8],
        before: &[u8],
        insert: &[u8],
        after: &[u8],
        downstream: &[u8],
    ) -> Vec<u8> {
        let mut query = upstream.to_vec();
        query.extend_from_slice(before);
        query.extend_from_slice(insert);
        query.extend_from_slice(after);
        query.extend_from_slice(downstream);
        query
    }

}
