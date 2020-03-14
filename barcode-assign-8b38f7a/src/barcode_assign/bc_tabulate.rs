use std::io::Write;

use failure;

use counts::*;

pub struct CLI {
    pub inputs: Vec<String>,
    pub output: String,
    pub mintotal: Option<usize>,
    pub minsamples: Option<usize>,
    pub mininsample: Option<usize>,
    pub omitfile: Option<String>,
}

impl CLI {
    pub fn run(&self) -> Result<(), failure::Error> {
        let mut counts = Vec::new();

        for input in self.inputs.iter() {
            let input_counts = SampleCounts::from_file(&input)?;
            counts.push((input.to_string(), input_counts));
        }

        let total_counts: SampleCounts = counts.iter().map(|(_input, counts)| counts).sum();
        let mut barcodes: Vec<(Vec<u8>, usize)> = total_counts.into_iter().collect();
        barcodes.sort_by_key(|(_barcode, counts)| -(*counts as isize));

        let mut out = std::fs::File::create(&self.output)?;
        write!(out, "barcode")?;
        for (input, _input_counts) in counts.iter() {
            write!(out, "\t{}", input)?;
        }
        write!(out, "\n")?;

        let mut omit_out: Box<dyn Write> = match &self.omitfile {
            Some(f) => Box::new(std::fs::File::create(f)?),
            None => Box::new(std::io::sink()),
        };

        for (barcode, _counts) in barcodes.iter() {
            let count_vec = SampleCounts::barcode_count_vec(
                counts.iter().map(|(_input, counts)| counts),
                barcode,
            );

            if self.is_omitted(&count_vec) {
                write!(omit_out, "{}\n", String::from_utf8_lossy(barcode))?;
            } else {
                write!(out, "{}", String::from_utf8_lossy(barcode))?;
                for ct in count_vec.iter() {
                    write!(out, "\t{}", ct)?;
                }
                write!(out, "\n")?;
            }
        }

        Ok(())
    }

    pub fn is_omitted(&self, count_vec: &Vec<usize>) -> bool {
        if let Some(mintotal) = self.mintotal {
            let total: usize = count_vec.iter().copied().sum();
            if total < mintotal {
                return true;
            }
        }

        if let Some(minsamples) = self.minsamples {
            let nsamples: usize = count_vec.iter().filter(|&&ct| ct > 0).count();
            if nsamples < minsamples {
                return true;
            }
        }

        if let Some(mininsample) = self.mininsample {
            if count_vec.iter().max().copied().unwrap_or(0) < mininsample {
                return true;
            }
        }

        false
    }
}

#[cfg(test)]
mod tests {
    extern crate tempfile;
    extern crate rand;

    use std::collections::HashSet;
    use std::fmt::Write;
    use std::ops::Deref;
    use bc_tabulate::tests::rand::distributions::Range;
    use bc_tabulate::tests::rand::distributions::Sample;
    use bc_tabulate::tests::rand::thread_rng;
    use bc_tabulate::tests::rand::Rng;
    
    use super::*;

    const NTS_LEN: usize = 4;
    static NTS: [u8; NTS_LEN] = [b'A', b'C', b'G', b'T'];
    
    fn barcode(l: usize) -> Vec<u8> {
        let mut bc = Vec::with_capacity(l);
        for _i in 0..l {
            bc.push(*thread_rng().choose(&NTS).unwrap());
        }
        bc
    }

    const BC_LEN: usize = 20;
    
    fn distinct_barcodes(n: usize) -> Vec<Vec<u8>> {
        let mut bcset = HashSet::new();
        while bcset.len() < n {
            bcset.insert(barcode(BC_LEN));
        }
        bcset.into_iter().collect()
    }

    fn barcode_counts<B: Deref<Target=[u8]>, I: IntoIterator<Item = B>>(barcodes: I) -> Vec<(Vec<u8>, usize)> {
        let mut bc_ct = Vec::new();
        let mut ctdist = Range::new(1, 1234);
        for bc in barcodes {
            bc_ct.push((bc.deref().to_vec(), ctdist.sample(&mut thread_rng())));
        }
        bc_ct
    }

    fn counts_to_sample<'a, I: IntoIterator<Item = &'a (Vec<u8>, usize)>>(bc_cts: I) -> SampleCounts {
        bc_cts.into_iter().flat_map(|(bc, ct)| std::iter::repeat(bc.as_slice()).take(*ct)).collect()
    }
    
    #[test]
    fn tabulate_one() {
        let mut ctvec_a = barcode_counts(distinct_barcodes(200));
        let cts_a: SampleCounts = counts_to_sample(ctvec_a.iter());

        let mut count_file = tempfile::NamedTempFile::new().unwrap();
        cts_a.write(&mut count_file).unwrap();
        let count_path = count_file.into_temp_path();

        let table_path = tempfile::NamedTempFile::new().unwrap().into_temp_path();
        
        let cli = CLI { inputs: vec![ count_path.to_string_lossy().into_owned() ],
                        output: table_path.to_string_lossy().into_owned(),
                        mintotal: None,
                        minsamples: None,
                        mininsample: None,
                        omitfile: None
        };

        cli.run().unwrap();
        
        let mut exp_table = String::new();
        writeln!(exp_table, "barcode\t{}", count_path.to_string_lossy()).unwrap();
        ctvec_a.sort_by_key(|(_bc, ct)| -(*ct as isize));
        for (bc, ct) in ctvec_a.iter() {
            writeln!(exp_table, "{}\t{}", String::from_utf8_lossy(bc), *ct).unwrap();
        }
        let mut exp_lines: Vec<&str> = exp_table.lines().collect();
        exp_lines.sort();
        
        let out_table = String::from_utf8(std::fs::read(table_path).unwrap()).unwrap();
        let mut out_lines: Vec<&str> = out_table.lines().collect();
        out_lines.sort();
        
        for (exp_line, act_line) in exp_lines.iter().zip(out_lines.iter()) {
            assert_eq!(exp_line, act_line);
        }
    }

    #[test]
    fn tabulate_two() {
        const N_A: usize = 50;
        const N_B: usize = 50;
        const N_AB: usize = 100;
        let all_barcodes = distinct_barcodes(N_A + N_B + N_AB);
        let barcodes_a = &all_barcodes[0..N_A];
        let barcodes_b = &all_barcodes[N_A..(N_A + N_B)];
        let barcodes_ab = &all_barcodes[(N_A + N_B)..(N_A + N_B + N_AB)];
        let ctvec_a_only = barcode_counts(barcodes_a.iter().map(Vec::as_slice));
        let ctvec_a_both = barcode_counts(barcodes_ab.iter().map(Vec::as_slice));
        let ctvec_b_only = barcode_counts(barcodes_b.iter().map(Vec::as_slice));
        let ctvec_b_both = barcode_counts(barcodes_ab.iter().map(Vec::as_slice));

        let cts_a: SampleCounts = counts_to_sample(ctvec_a_only.iter().chain(ctvec_a_both.iter()));
        let cts_b: SampleCounts = counts_to_sample(ctvec_b_only.iter().chain(ctvec_b_both.iter()));

        let mut count_a_file = tempfile::NamedTempFile::new().unwrap();
        cts_a.write(&mut count_a_file).unwrap();
        let count_a_path = count_a_file.into_temp_path();

        let mut count_b_file = tempfile::NamedTempFile::new().unwrap();
        cts_b.write(&mut count_b_file).unwrap();
        let count_b_path = count_b_file.into_temp_path();

        let table_path = tempfile::NamedTempFile::new().unwrap().into_temp_path();
        
        let cli = CLI { inputs: vec![ count_a_path.to_string_lossy().into_owned(),
                                      count_b_path.to_string_lossy().into_owned(), ],
                        output: table_path.to_string_lossy().into_owned(),
                        mintotal: None,
                        minsamples: None,
                        mininsample: None,
                        omitfile: None
        };

        cli.run().unwrap();
        
        let mut exp_table = String::new();
        writeln!(exp_table, "barcode\t{}\t{}",
                 count_a_path.to_string_lossy(),
                 count_b_path.to_string_lossy()).unwrap();
        for (bc, ct) in ctvec_a_only.iter() {
            writeln!(exp_table, "{}\t{}\t0", String::from_utf8_lossy(bc), *ct).unwrap();
        }
        for (bc, ct) in ctvec_b_only.iter() {
            writeln!(exp_table, "{}\t0\t{}", String::from_utf8_lossy(bc), *ct).unwrap();
        }
        for ((bca, cta), (bcb, ctb)) in ctvec_a_both.iter().zip(ctvec_b_both.iter()) {
            assert_eq!(bca, bcb);
            writeln!(exp_table, "{}\t{}\t{}",
                     String::from_utf8_lossy(bca),
                     *cta, *ctb).unwrap();
        }
        let mut exp_lines: Vec<&str> = exp_table.lines().collect();
        exp_lines.sort();
        
        let out_table = String::from_utf8(std::fs::read(table_path).unwrap()).unwrap();
        let mut out_lines: Vec<&str> = out_table.lines().collect();
        out_lines.sort();
        
        for (exp_line, act_line) in exp_lines.iter().zip(out_lines.iter()) {
            assert_eq!(exp_line, act_line);
        }
    }

    #[test]
    fn tabulate_many() {
        const N_BARCODE: usize = 200;
        const N_COUNT: usize = 4;
        const P_ZERO: f64 = 0.5;
        let all_barcodes = distinct_barcodes(N_BARCODE);
        let ctvecs: Vec<Vec<(Vec<u8>, usize)>>
            = std::iter::repeat_with(||
                                     barcode_counts(all_barcodes.iter().map(Vec::as_slice))
                                     .into_iter()
                                     .map(|(bc, ct)|
                                          (bc, if thread_rng().next_f64() < P_ZERO { 0 } else { ct }))
                                     .collect())
            .take(N_COUNT)
            .collect();

        let mut count_paths = Vec::new();
        
        for ctvec in ctvecs.iter() {
            let cts: SampleCounts = counts_to_sample(ctvec.iter().filter(|(_bc, ct)| *ct > 0));
            let mut count_file = tempfile::NamedTempFile::new().unwrap();
            cts.write(&mut count_file).unwrap();
            count_paths.push(count_file.into_temp_path());
        }

        let table_path = tempfile::NamedTempFile::new().unwrap().into_temp_path();

        let cli = CLI { inputs: count_paths.iter().map(|p| p.to_string_lossy().into_owned()).collect(),
                        output: table_path.to_string_lossy().into_owned(),
                        mintotal: None,
                        minsamples: None,
                        mininsample: None,
                        omitfile: None
        };

        cli.run().unwrap();

        let mut exp_table = String::new();

        write!(exp_table, "barcode").unwrap();
        for path in count_paths.iter() {
            write!(exp_table, "\t{}", path.to_string_lossy()).unwrap();
        }
        write!(exp_table, "\n").unwrap();
        for i in 0..ctvecs[0].len() {
            if ctvecs.iter().all(|ctvec| ctvec[i].1 == 0) {
                continue;
            }
            write!(exp_table, "{}", String::from_utf8_lossy(&ctvecs[0][i].0)).unwrap();
            for ctvec in ctvecs.iter() {
                write!(exp_table, "\t{}", ctvec[i].1).unwrap();
            }
            write!(exp_table, "\n").unwrap();
        }
        let mut exp_lines: Vec<&str> = exp_table.lines().collect();
        exp_lines.sort();
        
        let out_table = String::from_utf8(std::fs::read(table_path).unwrap()).unwrap();
        let mut out_lines: Vec<&str> = out_table.lines().collect();
        out_lines.sort();
        
        for (exp_line, act_line) in exp_lines.iter().zip(out_lines.iter()) {
            assert_eq!(exp_line, act_line);
        }
    }        
}
