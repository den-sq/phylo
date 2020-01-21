///! Tools for fetch, parsing, and converting DSSP Data.
use std::{io, fs::File, io::{BufRead, BufReader}};
use std::str;
use phf::{phf_map, Map};

#[derive(Debug)]
pub struct Residue {
	pub id: usize,
	pub amino_acid: char,
	pub chain: char,
	pub ss: char,
	pub ss_detail: String,
	pub bridge: Option<(usize, usize, char)>,
	pub solvent_accesibility: f64,
	pub h_bonds: [(i32, f32); 4],
	pub tco: f64,
	pub bend_angle: f64,
	pub dihedral_angle: f64,
	pub backbone: (f64, f64),
	pub alpha_carbon: (f64, f64, f64),
}

#[derive(Debug)]
/// Parent object matching the data in a DSSP file.
pub struct DSSP {
	pub id: String,
	pub classification: String,
	pub authors: Vec<String>,
	pub reference: String,
	pub residues: Vec<Residue>
}
impl DSSP {
	/// Standard maximum solvent exposure values from Tien et. al.
	/// 
	/// Only contains standard 20 amino acids.
	const TIEN: Map<&'static str, f64> = phf_map! {
        "ALA" => 129.0, "ARG" => 274.0, "ASN" => 195.0, "ASP" => 193.0, "CYS" => 167.0,
        "GLN" => 225.0, "GLU" => 223.0, "GLY" => 104.0, "HIS" => 224.0, "ILE" => 197.0,
        "LEU" => 201.0, "LYS" => 236.0, "MET" => 224.0, "PHE" => 240.0, "PRO" => 159.0,
        "SER" => 155.0, "THR" => 172.0, "TRP" => 285.0, "TYR" => 263.0, "VAL" => 174.0
    };

	/// Maps 1 character amino acid abbreviations to their 3 character equivalents.
	/// 
	/// Only contains standard 20 amino acids.
    const RES_1_TO_3: Map<char, &'static str> = phf_map! {
        'A' => "ALA", 'C' => "CYS", 'D' => "ASP", 'E' => "GLU", 'F' => "PHE",
        'G' => "GLY", 'H' => "HIS", 'I' => "ILE", 'K' => "LYS", 'L' => "LEU",
        'M' => "MET", 'N' => "ASN", 'P' => "PRO", 'Q' => "GLN", 'R' => "ARG",
        'S' => "SER", 'T' => "THR", 'V' => "VAL", 'W' => "TRP", 'Y' => "TYR"
    };

	/// Creates a new, empty DSSP object.
	pub fn bare() -> DSSP {
		DSSP {
			id: String::new(),
			classification: String::new(),
			authors: Vec::new(),
			reference: String::new(),
			residues: Vec::new()
		}
	}

	/// Gets the secondary structure all residues for the given protein chain.
	/// 
	/// Returns None if there is no amino acids for the given chain.
	pub fn ss(&self, chain: char) -> Option<Vec<char>> {
		let mut ss_list = Vec::new();
		for residue in &self.residues {
			if residue.chain == chain { ss_list.push(residue.ss) }
		}
		if !ss_list.is_empty() {
			Some(ss_list)
		} else { None }
	}

	/// Gets the % of maximum solvent accessibility for all residues for the given protein chain.
	/// 
	/// Calulates based on maximum values from Tien et. al.
	/// 
	/// Returns None if there is no amino acids for the given chain.
	pub fn sa(&self, chain: char) -> Option<Vec<f64>> {
		let mut sa_list = Vec::new();
		for residue in &self.residues {
			if residue.chain == chain { 
				sa_list.push(
					match DSSP::RES_1_TO_3.get(&residue.amino_acid) {
						Some(x) => residue.solvent_accesibility / DSSP::TIEN[x],
						None => -1.0
					});
			}
		}
		if !sa_list.is_empty() {
			Some(sa_list)
		} else { None }
	}

	/// Gets the amino acid sequence for the given protein chain.
	/// 
	/// Returns None if there is no amino acids for the given chain.
	pub fn seq(&self, chain: char) -> Option<Vec<char>> {
		let mut seq = Vec::new();
		for residue in &self.residues {
			if residue.chain == chain {
				seq.push(residue.amino_acid);
			}
		}
		if !seq.is_empty() {
			Some(seq)
		} else { None }
	}

	/// Reads a DSSP file.
	/// 
	/// 
	/// # TODO:
	/// 
	/// Following sections are ignored:
	/// COMPND
	/// SOURCE
	/// REFERE
	/// 
	/// # Errors:
	/// 
	/// Returns IO Error if file cannot be opened.
	/// 
	/// # Panics:
	/// 
	/// Will likely panic if file is not in valid DSSP format.
	pub fn read(file_name: &str) -> Result<DSSP, io::Error> {
		let mut data = DSSP::bare();
		let parse_separator_list = | field: &mut Vec<String>, separator: char, line: &str | {
            let mut initial_string = if !field.is_empty() {
                field.pop().unwrap()
            } else {
                String::new()
            };
            for slice in line.split(separator) {
                initial_string.push_str(slice);
                field.push(initial_string.clone());
                initial_string = String::new();
            }
        };
		fn atov<T: str::FromStr> (a: &str) -> T 
			where <T as std::str::FromStr>::Err: std::fmt::Debug { 
			a.trim_start().parse().unwrap() 
		};

		match File::open(file_name) {
			Ok(file) => {
				let buf_reader = BufReader::new(file);
				let mut data_section = false;
				for line_result in buf_reader.lines() {
					let line = match line_result {
						Ok(value) => value,
						Err(_) => continue
					};
					if line.len() < 6 { continue }
					match &line[0..6] {
						"==== S" => continue,
						"REFERE" => continue,
						"HEADER" => {
							data.classification = line[10..50].trim().to_string();
							// data.date = ???chrono?? line[50..59] TODO
							data.id = line[62..66].to_string();
						},
						"COMPND" => continue,
						"SOURCE" => continue,
						"AUTHOR" => parse_separator_list(&mut data.authors, ',', &line[10..line.len()]), // Could ignore breaks after a comma for faster read I think
						_ => {
							if !data_section {
								
								if &line[5..12] == "RESIDUE" {
									data_section = true;
									println!("Here");
								} else {
									println!("{}|RESIDUE", &line[5..12]);
								}
							} else {
								println!("{}", line);
								let char_list: Vec<char> = line.chars().collect();
								if char_list[13] == '!' {
									continue; // Chain break residue; switching to new chain.
								}
								data.residues.push(Residue {
									id: atov(&line[5..10]), // May contain letters? Need Sample
									chain: char_list[11],
									amino_acid: char_list[13],
									ss: char_list[15],
									ss_detail: String::from(&line[17..25]),
									bridge: Some((
										atov(&line[25..29]),
										atov(&line[29..33]),
										char_list[33]
									)),
									solvent_accesibility: atov(&line[34..38]),
									h_bonds: [
										(atov::<i32>(&line[39..45]), atov::<f32>(&line[46..50])),
										(atov::<i32>(&line[50..56]), atov::<f32>(&line[57..61])),
										(atov::<i32>(&line[61..67]), atov::<f32>(&line[68..72])),
										(atov::<i32>(&line[72..78]), atov::<f32>(&line[79..83]))
									],
									tco: atov(&line[85..91]),
									bend_angle: atov(&line[91..97]),
									dihedral_angle: atov(&line[97..103]),
									backbone: (atov::<f64>(&line[103..109]), atov::<f64>(&line[109..115])),
									alpha_carbon: (atov(&line[115..122]), atov(&line[122..129]), atov(&line[129..136]))
								});
							}
						}
					}
				}
				Ok(data)
			}, 
			Err(e) => Err(e)
		}
	}
}
