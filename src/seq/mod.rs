use serde::{Serialize, Deserialize};
use std::io;
use std::fs::File;
use std::io::prelude::*;

pub trait SeqData: Clone {
	fn len(&self) -> usize;
	fn id(&self) -> &str;
	fn seq(&self) -> Vec<char>;
	fn is_empty(&self) -> bool;
	//fn from(&self, id: String, sites: Vec<char>);
}


#[derive(Serialize, Deserialize, Debug, Clone)]
pub struct Seq {
	pub sites: Vec<char>,
	pub id: String
}

impl SeqData for Seq {
	fn len(&self) -> usize { self.sites.len() }
	fn id(&self) -> &str { &self.id }
	fn seq(&self) -> Vec<char> { self.sites.clone() }
	fn is_empty(&self) -> bool { self.sites.is_empty() }
	//fn from(&self, id: String, sites: Vec<char>) { self.id = id; self.sites = sites; }
}

impl Seq {
	pub fn new(id: String, sites: Vec<char>) -> Seq {
		Seq {
			id,
			sites
		}
	}

	pub fn serialize_fasta(&self) -> String {
		format!(">{}\n{}", self.id(), self.seq().into_iter().collect::<String>())
	}

	pub fn write(&self, file_name: &str, file_type: &str) -> Result<(), io::Error> {
        match File::create(file_name) {
            Ok(mut handle) => {
                let write_string = match file_type {
					"fasta" => self.serialize_fasta(),
                    _ => return Err(io::Error::new(io::ErrorKind::InvalidInput, "Invalid Format"))
                };
                
                handle.write_all(write_string.as_bytes())?;
            },
            Err(e) => return Err(e)
        }
        Ok(())
    }
}


#[derive(Serialize, Deserialize, Debug, Clone)]
pub struct QCSite {
    pub val: char,
    pub qual: u32
}

#[derive(Serialize, Deserialize, Debug, Clone)]
pub struct QCSeq {
	pub sites: Vec<QCSite>,
	pub id: String
}

impl SeqData for QCSeq {
	fn len(&self) -> usize { (&(self.sites)).len() }
	fn id(&self) -> &str { &self.id }
	fn seq(&self) -> Vec<char> { 
		self.sites.iter().map(|x| x.val).collect()
	}
	fn is_empty(&self) -> bool { (&(self.sites)).is_empty() }
}

#[derive(Serialize, Deserialize, Debug, Clone)]
pub struct ProbabalisticSite {
	pub val: char,
	pub prob: f64
}

#[derive(Serialize, Deserialize, Debug, Clone)]
pub struct ProbabalisticSeq {
	pub sites: Vec<Vec<ProbabalisticSite>>,
	pub id: String
}

impl ProbabalisticSeq {
	fn likely_at(&self, loc: usize) -> &ProbabalisticSite {
		let mut likely_site = &ProbabalisticSite{val: '-', prob: 0.0};
		for site in &(self.sites)[loc] {
			if site.prob > likely_site.prob {
				likely_site = site
			}
		}
		likely_site
	}
}

impl SeqData for ProbabalisticSeq {
	fn len(&self) -> usize { (&(self.sites)).len() }
	fn id(&self) -> &str { &self.id }
	fn seq(&self) -> Vec<char> {		
		let mut likely_seq: Vec<char> = Vec::new();
		for i in 0..(&(self.sites)).len() {
			likely_seq.push(self.likely_at(i).val);
		}
		likely_seq
	}
	fn is_empty(&self) -> bool { (&(self.sites)).is_empty() }
}