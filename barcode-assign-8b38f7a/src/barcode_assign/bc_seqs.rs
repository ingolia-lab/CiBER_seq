use std::collections::HashMap;
use std::fs::File;
use std::io::Write;

use bio::io::fastq;

use fastq_pair;
use neighborhood::*;

#[derive(Debug)]
pub struct Config {
    pub barcode_fastq: String,
    pub sequ_fastq: String,
    pub out_fastq: String,
    pub out_barcodes: Option<String>,
    pub out_barcode_freqs: Option<String>,
    pub neighborhood: Option<String>,
}

pub fn bc_seqs(config: Config) -> Result<(), failure::Error> {
    let barcode_reader = fastq::Reader::from_file(&config.barcode_fastq)?;
    let sequ_reader = fastq::Reader::from_file(&config.sequ_fastq)?;

    let mut fastq_writer = fastq::Writer::to_file(&config.out_fastq)?;

    let mut barcode_recs = HashMap::new();

    let pair_records =
        fastq_pair::PairRecords::new(barcode_reader.records(), sequ_reader.records());

    for pair_result in pair_records {
        let (barcode_record, sequ_record) = pair_result?;
        let barcode = barcode_record.seq().to_vec();
        let recs = barcode_recs.entry(barcode).or_insert_with(|| Vec::new());
        recs.push(sequ_record);
    }

    let bc_recs = if let Some(nbhd_filename) = config.neighborhood {
        let nbhds_raw = Neighborhood::gather_neighborhoods(barcode_recs);
        let mut nbhds: Vec<_> = nbhds_raw.into_iter().map(|n| n.into_sorted()).collect();
        nbhds.sort_unstable_by(|nbhdl, nbhdr| nbhdl.key_barcode().0.cmp(nbhdr.key_barcode().0));

        let mut nbhd_recs = Vec::new();
        let mut nbhd_counts = Vec::new();
        
        for nbhd in nbhds {
            nbhd_counts.push(nbhd.to_counts());

            let barcode = nbhd.key_barcode().0.to_vec();
            let recs = nbhd.into_barcodes().flat_map(|(_, recs)| recs.into_iter()).collect();

            nbhd_recs.push((barcode, recs));
        }

        SortedNeighborhood::write_tables(&nbhd_filename, nbhd_counts.iter())?;

        nbhd_recs
    } else {
        let mut bc_recs: Vec<(Vec<u8>, Vec<fastq::Record>)> = barcode_recs.into_iter().collect();
        bc_recs.sort_unstable_by(|(bcl, _recsl), (bcr, _recsr)| bcl.cmp(bcr));
        bc_recs
    };
    
    if let Some(barcode_filename) = config.out_barcodes {
        let barcode_writer = File::create(barcode_filename)?;
        write_barcode_table(barcode_writer, bc_recs.iter())?;
    }
    
    if let Some(freq_filename) = config.out_barcode_freqs {
        let freq_writer = File::create(freq_filename)?;
        write_freq_table(freq_writer, bc_recs.iter())?;
    }

    for (bc, recs) in bc_recs {
        let barcode = String::from_utf8(bc)?;
        
        for (recidx, rec) in recs.into_iter().enumerate() {
            let name = format!("{}_{}", barcode, recidx + 1);
            let named_record =
                fastq::Record::with_attrs(&name, None, rec.seq(), rec.qual());
            fastq_writer.write_record(&named_record)?;
        }
    }
        
    Ok(())
}

fn write_barcode_table<'i, W, I, A: 'i>(
    barcode_out: W,
    barcode_iter: I,
) -> Result<(), failure::Error>
where
    W: std::io::Write,
    I: Iterator<Item = &'i (Vec<u8>, Vec<A>)>
{
    let mut bcout = std::io::BufWriter::new(barcode_out);

    for (barcode, entries) in barcode_iter {
        bcout.write(barcode)?;
        bcout.write("\t".as_bytes())?;
        bcout.write(entries.len().to_string().as_bytes())?;
        bcout.write("\n".as_bytes())?;
    }

    Ok(())
}

fn write_freq_table<'i, W, I, A: 'i>(
    freq_out: W,
    barcode_iter: I,
) -> Result<(), failure::Error>
where
    W: std::io::Write,
    I: Iterator<Item = &'i (Vec<u8>, Vec<A>)>
{
    let mut fout = std::io::BufWriter::new(freq_out);

    let mut freq_counts = HashMap::new();

    for freq in barcode_iter {
        let freq_count = freq_counts.entry(freq.1.len()).or_insert(0);
        *freq_count += 1;
    }

    let mut freqs: Vec<usize> = freq_counts.keys().map(|&k| k).collect();
    freqs.sort();

    for freq in freqs {
        write!(fout, "{}\t{}\n", freq, freq_counts.get(&freq).unwrap_or(&0))?;
    }

    Ok(())
}

