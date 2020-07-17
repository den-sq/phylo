//! Tools for the use of alignments, including for alignments with multiple possible residues per position
//! or other metadata (such as quality metrics) associated at the sequence level.
use crate::utils::{arithmatic_mean, arithmatic_median, arithmatic_max, Ranges};
use crate::seq::{Seq, SeqData};
use std::io;
use std::fs::File;
use std::io::prelude::*;
use std::io::BufReader;
use std::io::SeekFrom;
use std::io::Cursor;
use std::fs::OpenOptions;
use std::process::{Command, Stdio};
use serde::{Serialize, Deserialize};

pub trait Aligner {
	fn align_simple_file(&self, file_name: &str, file_type: &str) -> Result<Alignment<Seq>, io::Error>;
	//fn align_seq<T: SeqData>(&self, seq_list: Vec<T>) -> Alignment<T>;
}

pub trait GroupAligner {
	fn align_simple_files<T: SeqData>(&self, file_one: &str, file_type_one: &str,
										file_two: &str, file_type_two: &str) -> Alignment<Seq>;
	//fn align_seq_groups<T: SeqData>(&self, group_one: Vec<T>, group_two: Vec<T>) -> Alignment<T>;
}

pub struct Mafft;

impl Aligner for Mafft {
	fn align_simple_file(&self, file_name: &str, _file_type: &str) -> Result<Alignment<Seq>, io::Error> {
		let mafft_run = Command::new("mafft")
			.args(&[file_name])
			.stdout(Stdio::piped())
			.output()
			.expect("Failed to Execute Process");
		Alignment::<Seq>::read_bytes(mafft_run.stdout, "fasta")
	}
}

#[derive(Serialize, Deserialize, Debug, Clone, Default)]
/// An alignment of some kind of sequence data.
pub struct Alignment<T: SeqData> {
	pub members: Vec<T>
}
impl<T:SeqData> Alignment<T> {
	/// Creates a new, empty alignment.
    pub fn new() -> Alignment<T> {
        Alignment::<T> {
			members: Vec::new()
        }
	}

	/// Creates a new Alignment from a vector of sequences that are already aligned.
	/// 
	/// Checks that all members are the same length, returns None otherwise.
    pub fn from(sequences: Vec<T>) -> Option<Alignment<T>> {
		if sequences.is_empty() {
			Some(Alignment::new())
		} else {
			let align_len = sequences[0].len();
			for seq in sequences.iter() {
				if seq.len() != align_len {
					return None
				}
			}
			Some(
				Alignment::<T> {
					members: sequences
			})
		}
	}

	/// Gets the mean, median, and max sequence identities for all the sequences when compare to the first sequence.
	pub fn get_root_identity(&self) -> Ranges {
		let base_seq = self.members[0].seq();
		let mut identities: Vec<f64> = Vec::new();
		println!("{}",self.members.len());
		for i in 1..self.members.len() {
			let compare_seq = self.members[i].seq();
			let mut match_count = 0;
			let mut residue_count = 0;
			let mut base_started = false;
			let mut compare_started = false;
			let mut unmatched_base_residues = 0;
			let mut unmatched_compare_residues = 0;
			for j in 0..self.len() {
				match (base_seq[j], compare_seq[j]) {
					('-', '-') => {}, // Do Nothing, it's an induced alignment gap
					('-', _c) => {
						// We've stared the compare string
						compare_started = true;

						if base_started { 
							// If we're not at the start of the sequence, start measuring the unmatched residues.
							unmatched_compare_residues += 1;

							// Since we know we haven't passed the end of the compare string, add in the unmatched base residues.
							residue_count += unmatched_base_residues;
							unmatched_base_residues = 0;
						};
					},
					(_b, '-') => {
						// We've stared the base string
						base_started = true;

						if compare_started { 
							// If we're not at the start of the sequence, start measuring the unmatched residues.
							unmatched_base_residues += 1;

							// Since we know we haven't passed the end of the base string, add in the unmatched compare residues.
							residue_count += unmatched_compare_residues;
							unmatched_compare_residues = 0;
						};

					},
					(b, c) => {
						// Both have starte, clear out the bases and increment.
						base_started = true;
						compare_started = true;
						residue_count += unmatched_base_residues;
						residue_count += unmatched_compare_residues;
						unmatched_base_residues = 0;
						unmatched_compare_residues = 0;
						residue_count += 1;
						
						if b == c {
							// If we have a match, increment.
							match_count += 1;
						}
					}
				}
			}

			//println!("{}:{}", match_count, residue_count);
			if residue_count == 0 {
				
			} else {
				identities.push(match_count as f64 / residue_count as f64)
			}
		}

		Ranges {
			mean: arithmatic_mean(identities.iter()),
			median: arithmatic_median(&mut identities.clone()).unwrap(),
			max: arithmatic_max(&identities).unwrap()
		}
	}
	/// Gets the mean, median, and max sequence identities for all the sequences when compare to the first sequence.
	/// 
	pub fn get_identity(&self) -> Ranges {
		let mut identities: Vec<f64> = Vec::new();
		//println!("{}",self.members.len());
		for k in 0..self.members.len() - 1 {
			let base_seq = self.members[k].seq();
			for i in (k + 1)..self.members.len() {
				let compare_seq = self.members[i].seq();
				let mut match_count = 0;
				let mut residue_count = 0;
				let mut base_started = false;
				let mut compare_started = false;
				let mut unmatched_base_residues = 0;
				let mut unmatched_compare_residues = 0;
				for j in 0..self.len() {
					match (base_seq[j], compare_seq[j]) {
						('-', '-') => {}, // Do Nothing, it's an induced alignment gap
						('-', _c) => {
							// We've stared the compare string
							compare_started = true;

							if base_started { 
								// If we're not at the start of the sequence, start measuring the unmatched residues.
								unmatched_compare_residues += 1;

								// Since we know we haven't passed the end of the compare string, add in the unmatched base residues.
								residue_count += unmatched_base_residues;
								unmatched_base_residues = 0;
							};
						},
						(_b, '-') => {
							// We've stared the base string
							base_started = true;

							if compare_started { 
								// If we're not at the start of the sequence, start measuring the unmatched residues.
								unmatched_base_residues += 1;

								// Since we know we haven't passed the end of the base string, add in the unmatched compare residues.
								residue_count += unmatched_compare_residues;
								unmatched_compare_residues = 0;
							};

						},
						(b, c) => {
							// Both have starte, clear out the bases and increment.
							base_started = true;
							compare_started = true;
							residue_count += unmatched_base_residues;
							residue_count += unmatched_compare_residues;
							unmatched_base_residues = 0;
							unmatched_compare_residues = 0;
							residue_count += 1;
							
							if b == c {
								// If we have a match, increment.
								match_count += 1;
							}
						}
					}
				}

				//println!("{}:{}", match_count, residue_count);
				if residue_count == 0 {
					
				} else {
					identities.push(match_count as f64 / residue_count as f64)
				}
			}
		}

		Ranges {
			mean: arithmatic_mean(identities.iter()),
			median: arithmatic_median(&mut identities.clone()).unwrap(),
			max: arithmatic_max(&identities).unwrap()
		}
	}


	/// List of IDs of sequences in alignment
	pub fn id_list(&self) -> Vec<String> {
		self.members.iter().map(|x| x.id().to_string()).collect()
	}

	/// Returns sequences not shared between the two alignments based on sequence id.
	pub fn unshared(&self, other: Alignment<T>) -> Vec<String> {
		let mut other_ids = other.id_list();
		let mut unshared_ids: Vec<String> = Vec::new();

		for id in self.id_list() {
			if let Some(x) = other_ids.iter().position(|x| *x == id) {
				other_ids.remove(x);
			} else {
				unshared_ids.push(id);
			}
		}
		unshared_ids.append(&mut other_ids);
		unshared_ids
	}

	/// Returns alignment filtered by given list of ids.
	pub fn filter_by_id(&self, filter: Vec<String>) -> Alignment<T> {
		let mut new_align: Alignment<T> = Alignment::new();
		for id in filter {
			if let Some(x) = self.members.iter().position(|x| x.id() == id) {
				new_align.members.push(self.members[x].clone());
			}
		}

		new_align
	}

	/// Returns a Vec of sequences that is the alignment with gaps removed.
	/// Collapses probabalistic sequences into single-sites.
	pub fn degap(&self) -> Vec<Seq> {
		let mut seq_set: Vec<Seq> = Vec::new();
		for seq in &self.members {
			let mut compact_seq = seq.seq();
			compact_seq.retain(|x| *x != '-');
			seq_set.push(Seq::new(
				seq.id().to_string(),
				compact_seq
			));
		}
		seq_set
	}

	pub fn likely_alignment(&self) -> Alignment<Seq> {
		let mut likely_seqs: Vec<Seq> = Vec::new();
		for seq in &self.members {
			likely_seqs.push(Seq {
				sites: seq.seq(),
				id: seq.id().to_string()
			});
		}

		// TODO
		
		Alignment::from(likely_seqs).unwrap()
	}

	/// Gets length of alignment (# of characters in alignment length)
	pub fn len(&self) -> usize {
		self.members[0].len()
	}

	/// Returns whether there are no alignment members.
	pub fn is_empty(&self) -> bool {
		if self.members.is_empty() { true }
		else { self.members[0].is_empty() }
	}

	/// Creates an new Alignment from an unaligned vector of sequences.
	pub fn build_from<U: Aligner>(sequences: Vec<Seq>, aligner: U) -> Result<Alignment<Seq>, io::Error> {
		let mut seq_fasta = String::new();
		for seq in sequences {
			seq_fasta.push_str(&seq.serialize_fasta());
		}
		
		match File::create("mafft_tmp.fa") {
            Ok(mut handle) => {
				handle.write_all(seq_fasta.as_bytes())?;
            },
            Err(e) => { return Err(e); }
		}
		aligner.align_simple_file("mafft_tmp.fa", "fasta")
	}

    /// Parse a phylip file into alignment, with loose standards.
	fn parse_phylip_loose<U: Read + Seek>(mut reader: BufReader<U>) -> Result<Alignment<Seq>, io::Error> {
		let mut line = String::new();
		let mut header_line = String::new();

		// Get size parameters.
		reader.read_line(&mut header_line)?;
		let mut v = header_line.trim().split_whitespace();
		let seq_count: u32 = match v.next().unwrap().parse::<u32>() {
			Ok(count) => count,
			Err(_) => return Err(io::Error::new(io::ErrorKind::InvalidInput, 
                    format!("Invalid Phylip Format, Invalid Header Section: {}.", header_line)))
		};
		let seq_len: u32 = match v.next().unwrap().parse::<u32>() {
			Ok(length) => length,
			Err(_) => return Err(io::Error::new(io::ErrorKind::InvalidInput, 
                    format!("Invalid Phylip Format, Invalid Header Section: {}.", header_line)))
		};
		

		// Skip blank lines.
		while reader.read_line(&mut line).is_ok() {
			if line != "\n" {
				reader.seek(SeekFrom::Current(-(line.len() as i64)))?;
				break;
			}
			line = String::new();
		}

		let mut data: Vec<Seq> = Vec::new();

		// Read alignment lines.  Allow spaces between ID and sequence and between codons.
		for _i in 0..seq_count {
			let mut line = String::new();
			reader.read_line(&mut line)?;
			let id = line[0..10].trim();
			let mut seq: Vec<char> = Vec::new();
			for c in line[10..].chars() {
				match c {
					' ' | '\u{A0}' | '\n' | '\t' => continue,
					_ => seq.push(c)
				}
			}
			if seq.len() as u32 == seq_len {
				data.push(Seq::new(String::from(id), seq));
			} else {
				return Err(io::Error::new(io::ErrorKind::InvalidInput, 
                    format!("Invalid Phylip Format, Invalid Alignment Line: {}.  Sequence Length Found {} vs {} expected.", 
								line, seq.len(), seq_len)))
			}
		};

		Ok(Alignment::from(data).unwrap())
	}

    /// Parses a phylip file into alignment object, following strict rules
	fn parse_phylip_strict<U: Read>(mut reader: BufReader<U>) -> Result<Alignment<Seq>, io::Error> {
		let mut line = String::new();
		let mut header_line = String::new();

		// Get size parameters.
		reader.read_line(&mut header_line)?;
		let mut v = header_line.trim().split_whitespace();
		let seq_count: u32 = match v.next().unwrap().parse::<u32>() {
			Ok(count) => count,
			Err(_) => return Err(io::Error::new(io::ErrorKind::InvalidInput, 
                    format!("Invalid Phylip Format, Invalid Header Section: {}.", header_line)))
		};
		let seq_len: usize = match v.next().unwrap().parse::<usize>() {
			Ok(length) => length,
			Err(_) => return Err(io::Error::new(io::ErrorKind::InvalidInput, 
                    format!("Invalid Phylip Format, Invalid Header Section: {}.", header_line)))
		};

		let mut data: Vec<Seq> = Vec::new();

		// Read alignment lines.  Strict adherence to 10 id chars then alignment.
		for _i in 0..seq_count {
			reader.read_line(&mut line)?;
			data.push(Seq::new(line[0..10].trim().to_string(), line[10..(seq_len + 10)].chars().collect()))
		}

		Ok(Alignment::from(data).unwrap())
	}

	/// Parses fasta data from a buffered reader into an alignment object.
	fn parse_fasta<U: Read>(mut reader: BufReader<U>) -> Result<Alignment<Seq>, io::Error> {
		let mut line = String::new();
		let mut data: Vec<Seq> = Vec::new();
		let mut cur_seq: Vec<char> = Vec::new();
		let mut name = String::new();
		while reader.read_line(&mut line)? > 0 {
			if line.starts_with('>') {
				if !cur_seq.is_empty() {
					data.push(Seq {
						id: name,
						sites: cur_seq
					});
					cur_seq = Vec::new();
				}
				name = line.clone()[1..line.len() - 1].to_string();
			} else {
				let cap = line.len() - 1;
				cur_seq.append(&mut line[..cap].chars().collect::<Vec<char>>());
			}
			line = String::new();
		}
		data.push(Seq {
			id: name,
			sites: cur_seq
		});

		Ok(Alignment::from(data).unwrap())
	}

	/// Reads alignment from a vector of bytes into an alignment object.
	/// 
	/// _strict alignment types must precisely match file format.
	/// _loose are more forgiving but must still be parseable.
	/// 
	/// # Errors
	/// 
	/// align_type must be "phylip_loose", "phylip", "phylip_strict", or "fasta".
	pub fn read_bytes(new_align: Vec<u8>, align_type: &str) -> Result<Alignment<Seq>, io::Error> {
        let parse_func = match align_type {
            "phylip_loose" => Alignment::<Seq>::parse_phylip_loose,
			"phylip" | "phylip_strict" => Alignment::<Seq>::parse_phylip_strict,
			"fasta" => Alignment::<Seq>::parse_fasta,
            _ => return Err(io::Error::new(io::ErrorKind::InvalidInput, "Invalid Format"))
        };

    	let buf_reader = BufReader::new(Cursor::new(new_align));
        parse_func(buf_reader)
	}

	/// Reads alignment from file into an alignment object.
	/// 
	/// _strict file types must precisely match file format.
	/// _loose are more forgiving but must still be parseable.
	/// 
	/// # Errors
	/// 
	/// file_name must be valid path to an existing file that must match the file format specification.
	/// file_type must be "phylip_loose", "phylip", "phylip_strict", or "fasta".
	pub fn read(file_name: &str, file_type: &str) -> Result<Alignment<Seq>, io::Error> {
        let parse_func = match file_type {
            "phylip_loose" => Alignment::<Seq>::parse_phylip_loose,
			"phylip" | "phylip_strict" => Alignment::<Seq>::parse_phylip_strict,
			"fasta" => Alignment::<Seq>::parse_fasta,
            _ => return Err(io::Error::new(io::ErrorKind::InvalidInput, "Invalid Format"))
        };

        match File::open(file_name) {
            Ok(file) => {
                let buf_reader = BufReader::new(file);
                parse_func(buf_reader)
            }
            Err(e) => Err(e)
        }
    }

	/// Serializes alignment to a strict Phylip format string, truncating data where necessary.
	/// 
	/// Uses strict phylip format, which truncates ID to 10 characters.
	/// Returns empty string if there are no members of the alignment.
	fn serialize_phylip(&self) -> String {
		if !self.members.is_empty() {
			let mut serial = format!("{} {}\n", self.members.len(), self.members[0].len());
            for seq in &self.members {
                serial.push_str(&(seq.id().to_owned() + "          ")[0..10]);
                serial.push_str(&(seq.seq().iter().collect::<String>()));
                serial.push('\n');
            }
            serial
		} else { String::new() }
	}

	/// Serializes alignment to a loosely matched phylip format string.
	/// 
	/// Loose Phylip format means the ID may be of any length, instead of only 10 characters.
	/// Returns empty string if there are no members of the alignment.
	fn serialize_phylip_loose(&self) -> String {
		if !self.members.is_empty() {
			let mut serial = format!("{} {}\n", self.members.len(), self.members[0].len());
            for seq in &self.members {
			    serial.push_str(&format!("{} {}", seq.id(), seq.seq().iter().collect::<String>()));
		    }
            serial
		}  else { String::new() }
	}

	/// Serialzes alignment to fasta format, as a string.
	/// 
	/// Returns empty string if there are no members of the alignment.
	/// Cuts lines at 80 characters, unless there are less than 80 characters left.
	fn serialize_fasta(&self) -> String {
		let mut serial = String::new();
		for seq in &self.members {
			let mut pos = 0;			
			let residues = seq.seq();
			serial.push_str(&format!(">{}\n", seq.id()));
			while (pos + 80) < residues.len() {
				serial.push_str(&format!("{}\n", residues[pos..pos + 80].iter().collect::<String>()));
				pos += 80;
			}
			serial.push_str(&format!("{}\n",residues[pos..residues.len()].iter().collect::<String>()));
		}
		serial
	}


	/// Writes an Alignment object to a file.
	/// 
	/// Uses serialization functions to generate data to write.
	/// Uses "strict" versions of file format if not specified.
	/// 
	/// # Errors
	/// file_name must be writeable.
	/// file_type must be "phylip", "phylip_loose", "phylip_strict", or "fasta".
    pub fn write(&self, file_name: &str, file_type: &str) -> Result<(), io::Error> {
        match File::create(file_name) {
            Ok(mut handle) => {
                let write_string = match file_type {
                    "phylip" | "phylip_strict" => self.serialize_phylip(),
					"phylip_loose" => self.serialize_phylip_loose(),
					"fasta" => self.serialize_fasta(),
                    _ => return Err(io::Error::new(io::ErrorKind::InvalidInput, "Invalid Format"))
                };
                
                handle.write_all(write_string.as_bytes())?;
            },
            Err(e) => return Err(e)
        }
        Ok(())
    }

	/// Apends the alignment to a file.
	/// 
	/// Uses serialization functions to generate data to write.
	/// Uses "strict" versions of file format if not specified.
	/// 
	/// # Errors
	/// file_name must be writeable.
	/// file_type must be "phylip", "phylip_loose", "phylip_strict", or "fasta".
	pub fn append(&self, file_name: &str, file_type: &str)-> Result<(), io::Error> {
		match OpenOptions::new().append(true).create(true).open(file_name) {
            Ok(mut handle) => {
                let write_string = match file_type {
                    "phylip" | "phylip_strict" => self.serialize_phylip(),
					"phylip_loose" => self.serialize_phylip_loose(),
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
