use rust_htslib::bam;
use rust_htslib::bam::Read;

use failure::*;

pub fn read_tag(r1: &bam::Record) -> Option<&[u8]> {
    if let Some(delim_pos) = r1.qname().iter().position(|&ch| ch == b'_') {
        Some(r1.qname().split_at(delim_pos).0)
    } else {
        None
    }
}

pub struct BarcodeGroups<'a> {
    bam_reader: &'a mut bam::Reader,
    next_record: Option<bam::Record>,
    read_barcode: &'a dyn Fn(&bam::Record) -> Option<&[u8]>,
}

impl<'a> BarcodeGroups<'a> {
    pub fn new_with_read_names(bam_reader: &'a mut bam::Reader) -> Result<Self, failure::Error> {
        Self::new(&read_tag, bam_reader)
    }

    pub fn new(
        read_barcode: &'a dyn Fn(&bam::Record) -> Option<&[u8]>,
        bam_reader: &'a mut bam::Reader,
    ) -> Result<Self, failure::Error> {
        let mut bg = BarcodeGroups {
            bam_reader: bam_reader,
            next_record: None,
            read_barcode: read_barcode,
        };
        bg.next_record = bg.read_next_record()?;
        Ok(bg)
    }

    fn read_next_record(&mut self) -> Result<Option<bam::Record>, failure::Error> {
        let mut rec = bam::Record::new();
        match self.bam_reader.read(&mut rec) {
            Ok(true) => Ok(Some(rec)),
            Ok(false) => Ok(None),
            Err(e) => Err(e.into()),
        }
    }

    fn barcode_group(
        &mut self,
        curr: bam::Record,
    ) -> Result<(Vec<u8>, Vec<bam::Record>), failure::Error> {
        let curr_bc = (self.read_barcode)(&curr)
            .ok_or(err_msg("No barcode for current read"))?
            .to_vec();
        let mut bc_group = Vec::new();
        bc_group.push(curr);

        loop {
            let next = self.read_next_record()?;
            if let Some(rec) = next {
                // Ends borrow of rec to analyze barcode
                let same_group = (self.read_barcode)(&rec)
                    .ok_or(format_err!("No barcode for read {:?}", rec))?
                    == curr_bc.as_slice();

                if same_group {
                    bc_group.push(rec);
                } else {
                    self.next_record = Some(rec);
                    break;
                }
            } else {
                self.next_record = None;
                break;
            }
        }

        Ok((curr_bc, bc_group))
    }
}

impl<'a> Iterator for BarcodeGroups<'a> {
    type Item = Result<(Vec<u8>, Vec<bam::Record>), failure::Error>;

    fn next(&mut self) -> Option<Self::Item> {
        if let Some(curr) = self.next_record.take() {
            Some(self.barcode_group(curr))
        } else {
            None
        }
    }
}
