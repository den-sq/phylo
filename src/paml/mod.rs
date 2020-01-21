//! Tools for the execution of PAML and parsing the output.

use std::fs::File;
use std::collections::HashMap;
use std::collections::HashSet;
use std::io::BufReader;
use std::io::prelude::*;
use std::io;
use ndarray::Array2;
use regex::RegexSet;
use crate::alignment::Alignment;
use crate::seq::Seq;
use crate::seq::ProbabalisticSeq;
use crate::seq::ProbabalisticSite;
use crate::clade::Clade;
use crate::prot::{ProteinData};
use std::process::Command;
use std::process::Stdio;
use serde::{Serialize, Deserialize};

#[derive(Debug)]
pub struct CodonTable { }

#[derive(Debug)]
/// Holds Data Related to One Section of a PAML Ouptut File.
/// 
/// PAML Sections have vastly different data so they are held in an enum.
pub enum PAMLSection
{
	Align(Alignment<Seq>),
	AlignWithAncestral(Alignment<Seq>),
	AncestralReconstruction(Clade),
	AncestralAccuracy((Vec<f64>, Vec<f64>)),
	MarginalReconstruction(Alignment<ProbabalisticSeq>, Vec<u16>),
	JointReconstruction(Vec<(f64, Vec<(Seq, f64)>)>),
	AncestralSequences(Vec<Seq>),
	AncestralResidues(Vec<Seq>),
	PerSiteChanges(Vec<u32>, Vec<(Codon, Codon)>, HashMap<String, f64>),
	ChangeSummary(HashMap<Endpoints, BranchChangeset>),
	SitePattern(Alignment<Seq>),
	Pattern(String),
	CodonUsage(CodonTable),
	CodonPos(Array2<f64>),
	CodonUseSum(CodonTable),
	CodonPosSum(Array2<f64>),
	CodonFreq(Vec<f64>),
	Tree(Box<PAMLTree>),
	Branches(HashMap<Endpoints, BranchData>),
	Empty()
}

impl PAMLSection
{
	/// Parses the PAML Section containing the Marginal Reconstruction of Ancestral Sequences.
	fn parse_marginal<U: Read>(reader: &mut BufReader<U>) -> Result<PAMLSection, io::Error> {
		let mut site = String::new();
		let mut seq_list: Vec<ProbabalisticSeq> = Vec::new();
		let mut freqs: Vec<u16> = Vec::new();
		let mut initial_seq: usize = 0;
		let mut anc_seq: usize = 0;

		reader.read_line(&mut site)?;
		if site == "\n" {
			return Err(io::Error::new(io::ErrorKind::InvalidInput, "Unreadable PAML Section."));
		}
		// Scope out this area so site_vec drops.
		{
			let mut site_vec = site.split_whitespace();
			site_vec.next(); 												// site counter
			freqs.push(site_vec.next().unwrap().parse::<u16>().unwrap()); 	// site frequency.

			// Get through the non-ancestral sequences.
			loop {
				if site_vec.next().unwrap() == ":" {
					break;
				}
				let res: Vec<char> = site_vec.next().unwrap().chars().collect();
				initial_seq += 1;
				seq_list.push(ProbabalisticSeq {
					sites: vec![vec![ProbabalisticSite { val: res[1], prob: 1.00 }]],
					id: initial_seq.to_string()
				});
			}
			// Get through the ancestral sequences.
			loop {
				if site_vec.next() == None {
					break;
				}
				let res: char = site_vec.next().unwrap().chars().next().unwrap();
			
				anc_seq += 1;
				

				seq_list.push(ProbabalisticSeq {
					sites: vec![vec![ProbabalisticSite { 
						val: res, 
						prob: match site_vec.next().unwrap().parse::<f64>() {
							Ok(prob) => prob,
							Err(_) => 0.000 // NaN; Right now default to 0.000
						} 
					}]],
					id: (initial_seq + anc_seq).to_string()
				});

				// I don't actually know what these values are for.
				site_vec.next();
				site_vec.next();
			}
		}
		site = String::new();
		reader.read_line(&mut site)?;
		while site != "\n" {
			// Now that the vectors have all be created, we can go through and push individually.
			{
				let mut site_vec = site.split_whitespace();
				site_vec.next(); 												// site counter
				freqs.push(site_vec.next().unwrap().parse::<u16>().unwrap()); 	// site frequency.

				// Get through the non-ancestral sequences.
				for seq in seq_list.iter_mut().take(initial_seq) {
					site_vec.next(); // Codon.
					let res: Vec<char> = site_vec.next().unwrap().chars().collect();
					seq.sites.push(vec![ProbabalisticSite { val: res[1], prob: 1.00 }]);
				}
				site_vec.next().unwrap();
				// Get through the ancestral sequences.
				for seq in seq_list.iter_mut().skip(initial_seq).take(anc_seq) {
					site_vec.next();
					let res: char = site_vec.next().unwrap().chars().next().unwrap();

					seq.sites.push(vec![ProbabalisticSite { 
						val: res, 
						prob: match site_vec.next().unwrap().parse::<f64>() {
							Ok(prob) => prob,
							Err(_) => 0.000 // NaN; Right now default to 0.000
						} 
					}]);

					// I don't actually know what the next values are for.
					site_vec.next();
					site_vec.next();
				}
			}
			site = String::new();
			reader.read_line(&mut site)?;
		}
		
		Ok(PAMLSection::MarginalReconstruction(Alignment::<ProbabalisticSeq> { 
			members: seq_list
		}, freqs))
	}

	/// Parses the summary of changes across an ancestor -> dependent for a single branch.
	fn parse_change_summary<U: Read>(reader: &mut BufReader<U>) -> Result<(Endpoints, BranchChangeset), io::Error> {
		let mut line = String::new();
		let ends: Endpoints;
		let mut change_set: BranchChangeset;
		reader.read_line(&mut line)?;
		
		while line == "\n" {
			// Read through the first few lines of preface until we hit a new line.
			line = String::new();
			reader.read_line(&mut line)?;
		}
		{
			let mut meta = line.split_whitespace();
			match meta.next() {
				Some("Branch") => {}, // We're good, keep going
				Some("\nBranch") => {}, // We're good, keep going
				// Otherwise We're not in the correct section.
				Some(_) => return Err(io::Error::new(io::ErrorKind::InvalidData, "Not a Change Summary Subsection.")),
				None => return Err(io::Error::new(io::ErrorKind::InvalidData, "Not a Change Summary Subsection.")),
			}
			
			meta.next();	// #:
			let mut endpoint_part = meta.next().unwrap().split("..");
			ends = Endpoints(
				endpoint_part.next().unwrap().parse::<u32>().unwrap(),
				endpoint_part.next().unwrap().parse::<u32>().unwrap()
			);
			// This is an annoying-ass line to parse.  Examples:
			// Branch 4:   32..3  (T504135325)  (n=528.0 s=18.0)
			// Branch 3:   31..32 (n= 0.0 s= 0.0)
			// So we have to play silly buggers now.
			let next = meta.next().unwrap();
			let next_vec = next.chars().collect::<Vec<char>>();
			let mut id: String;
			let mut n: &str;
			if next_vec[0] == '(' && next_vec[next_vec.len() - 1] == ')' {
				// We have an id, so assign it and grab next value for n.
				id = next.to_string();
				id = id[1..(id.len()-1)].to_string();
				n = match meta.next() {
					Some(val) => val.split('=').last().unwrap().trim(),
					None => return Err(io::Error::new(io::ErrorKind::InvalidInput, "Abrupt Ending to Change Summary SubSection.")) // Some files just end here.
				}
			} else {
				// We have no id, so assign the generated tree #; current value is start of n/s section.
				id = ends.1.to_string();
				n = next.split('=').last().unwrap().trim();
			};
			if n == "" {
				n = meta.next().unwrap(); // We have a space between = and number.
			}

			let mut s = meta.next().unwrap().split('=').last().unwrap().trim();
			if s == "" {
				s = meta.next().unwrap(); // We have a space between = and number.
			}
			change_set = BranchChangeset {
				n: n.parse::<f64>().unwrap(),
				s: s[0..(s.len() -1)].parse::<f64>().unwrap(),
				branch: ends.clone(),
				id,
				changes: HashMap::new()
			}
		}
		reader.read_line(&mut line)?;
		line = String::new();
		reader.read_line(&mut line)?;
		while line != "\n" {
			{
				// These lines can be with or without codons.
				// Examples:
				// 144 AAT (N) 0.710 -> CAT (N) 0.433
				// 	 1 GAG (V) 0.134 -> --- (*)
				// 149 L 0.321 -> M 0.391
				let meta: Vec<&str> = line.split_whitespace().collect();
				let id: u32 = meta[0].parse::<u32>().unwrap();
				let codon_push = if meta.len() < 7 { 0 } else { 1 };
				let mut ancestor = Codon {
									res: meta[1 + codon_push].chars().collect::<Vec<char>>()[codon_push],
									prob: meta[2 + codon_push].parse::<f64>().unwrap(),
									nucs: Default::default()
								};
				if codon_push > 0 {	ancestor.nucs.copy_from_slice(&meta[1].chars().collect::<Vec<char>>()[0..3]); }
				let mut descendent = Codon {
									res: meta[4 + 2 * codon_push].chars().collect::<Vec<char>>()[codon_push],
									prob: if meta.len() < (6 + 2 * codon_push) { 1.00 } 
										else { meta[5 + 2 * codon_push].parse::<f64>().unwrap() },
									nucs: Default::default()
								};
				if codon_push > 0 {	descendent.nucs.copy_from_slice(&(meta[5].chars().collect::<Vec<char>>()[0..3])); }
				change_set.changes.insert(id, (ancestor, descendent));
			}
			line = String::new();
			reader.read_line(&mut line)?;
		}
		Ok((ends, change_set))
	}

	/// Parses the PAML RST Section containing a summary of likely changes in the marginallly recostructed sequences.
	fn parse_change_summary_section<U: Read>(reader: &mut BufReader<U>) -> Result<PAMLSection, io::Error> {
		let mut all_changes: HashMap<Endpoints, BranchChangeset> = HashMap::new();
		loop {
			match PAMLSection::parse_change_summary(reader) {
				Ok((ends, change_set)) => all_changes.insert(ends, change_set),
				Err(e) => match e.kind() {
					io::ErrorKind::InvalidData => return Ok(PAMLSection::ChangeSummary(all_changes)),
					_ => return Err(e)
				}
			};
		}
	}

	/// Parses the section holding the calculated per-branch selection data.
	#[allow(non_snake_case)]
	fn parse_dNdS_section<U: Read>(reader: &mut BufReader<U>) -> Result<PAMLSection, io::Error> {
		let mut branch_data: HashMap<Endpoints, BranchData> = HashMap::new();
		let mut line = String::new();
		let mut eol_check = reader.read_line(&mut line)?;
		while (line != "\n") && (eol_check > 0) {
			{
				let mut meta = line.split_whitespace();
				let mut endpoint_part = meta.next().unwrap().split("..");
				let ends = Endpoints(
					endpoint_part.next().unwrap().parse::<u32>().unwrap(),
					endpoint_part.next().unwrap().parse::<u32>().unwrap()
				);
				branch_data.insert(ends, BranchData {
					t: meta.next().unwrap().parse::<f64>().unwrap(),
					N: meta.next().unwrap().parse::<f64>().unwrap(),
					S: meta.next().unwrap().parse::<f64>().unwrap(),
					dNdS: meta.next().unwrap().parse::<f64>().unwrap(),
					dN: meta.next().unwrap().parse::<f64>().unwrap(),
					dS: meta.next().unwrap().parse::<f64>().unwrap(),
					NdN: meta.next().unwrap().parse::<f64>().unwrap(),
					SdS: meta.next().unwrap().parse::<f64>().unwrap()
				});
			}
			line = String::new();
			eol_check = reader.read_line(&mut line)?;
		}
		Ok(PAMLSection::Branches(branch_data))
	}

	/// Fetches and parses the next section of a PAML output file.
	pub fn get_next_section<U: Read>(reader: &mut BufReader<U>) -> Result<PAMLSection, io::Error> {
		let mut line = String::new();
		// Create early and once.
		let set = RegexSet::new(&[
						r"^CODONML.*",
						r"^\nCODONML.*",		
						r"^TREE # *[0-9]*:.*\n",
						r"^\nTREE # *[0-9]*:.*\n",
						r"^Supplemental results for CODEML.*\n",
						r"^ *[0-9]+ *[0-9]+ ?\n",
						r"^\n *[0-9]+ *[0-9]+ ?\n",					// Start of Phylip Section.
						r"^\nTime used: .*\n",
						r"Time used: .*\n",
					// Codon Table
					]).unwrap();

		while reader.read_line(&mut line)? > 0 {
			// First, do the quick string comparsion.  Avoid regex if possible.
			match line.as_str() {
				//"\nPrinting out site pattern counts.\n" => return Ok(PAMLSection::Empty()),
				//"\nCodon position x base (3x4) table for each sequence.\n" => return Ok(PAMLSection::Empty()),
				//"\nSums of codon usage counts.\n" => return Ok(PAMLSection::Empty()),
				//"\nCodon position x base (3x4) table, overall.\n" => return Ok(PAMLSection::Empty()),
				//"\nCodon frequencies under model, for use in evolver (TTT TTC TTA TTG ... GGG):\n" => return Ok(PAMLSection::Empty()),
				//"\nDetailed output identifying parameters\n" => return Ok(PAMLSection::Empty()),
				"\ndN & dS for each branch\n" | "dN & dS for each branch\n" => {
					let mut blank = String::new();
					for _i in 0..3 {
						reader.read_line(&mut blank)?;
					}

					return PAMLSection::parse_dNdS_section(reader);
				}
				"\nAncestral reconstruction by CODONML.\n" | "Ancestral reconstruction by CODONML.\n" => {
					let mut tree = String::new();
					let sec: PAMLSection;

					// Reads the first tree, with IDs & branch lengths.
					if reader.read_line(&mut tree).is_ok() {
						reader.read_line(&mut tree)?;
						match Clade::parse_newick_string(tree) {
							Ok(ac_tree) => sec = PAMLSection::AncestralReconstruction(ac_tree),
							Err(e) => return Err(e)
						}
					} else {
						return Err(io::Error::new(io::ErrorKind::InvalidInput, "Unreadable Ancestral Reconstruction Section."))
					}
					tree = String::new();
					// Ignore the rest of the trees for now.
					for _i in 0..10 {
						reader.read_line(&mut tree)?;
						// 0,2,4,7 empty.
						// 1 numbered labels
						// 3 ancestral tree references.
						// 4-5 Rod Page Treeview
						// 8-9 ID ancestral nodes.
					}
					return Ok(sec)
				},
				"\n(1) Marginal reconstruction of ancestral sequences\n" => {
					let mut blank = String::new();
					for _i in 0..7 {
						reader.read_line(&mut blank)?;
					}

					return PAMLSection::parse_marginal(reader);
				}
				"\nSummary of changes along branches.\n" | "Summary of changes along branches.\n" => {
					let mut blank = String::new();
					for _i in 0..2 {
						reader.read_line(&mut blank)?;
					}

					return PAMLSection::parse_change_summary_section(reader);
				},
				//"\nList of extant and reconstructed sequences\n" => return Ok(PAMLSection::Empty()),
				//"\nOverall accuracy of the 28 ancestral sequences:\n" => return Ok(PAMLSection::Empty()),
				//"\nAmino acid sequences inferred by codonml.\n" => return Ok(PAMLSection::Empty()),
				//"\nCounts of changes at sites, listed by site\n" => return Ok(PAMLSection::Empty()),
				//"\n(2) Joint reconstruction of ancestral sequences\n" => return Ok(PAMLSection::Empty()),

				"\n" => continue,
				this_line => {
					match set.matches(this_line).iter().next() {
						Some(0) => (),
						Some(1) => (),
						Some(2) => (), // Handle Multiple Trees in the future?
						Some(3) => (),
						Some(8) | Some (9) => return Ok(PAMLSection::Empty()),
						Some(_) => (),
						None => ()
					}
				}
			}
			line = String::new();
		}

		Ok(PAMLSection::Empty())
	}
}

/// A single codon in a protein sequence, with associated nucleotides.
#[derive(Debug)]
pub struct Codon
{
	res: char,
	nucs: [char; 3],
	prob: f64
}

/// Set of changes along a branch.
#[derive(Debug)]
pub struct BranchChangeset
{
	n: f64,
	s: f64,
	branch: Endpoints,
	id: String,
	changes: HashMap<u32,(Codon, Codon)>
}

/// Pair of Branch Endpoints.
#[derive(Debug,Eq,PartialEq,Hash,Clone)]
pub struct Endpoints(u32, u32);

/// Branch Selection Data.
#[derive(Debug)]
#[allow(non_snake_case)]
pub struct BranchData 
{
	t: f64,
	N: f64,
	S: f64,
	dNdS: f64,
	dN: f64,
	dS: f64,
	NdN: f64,
	SdS: f64
}

/// PAML Tree Section.
#[derive(Debug)]
#[allow(non_snake_case)]
pub struct PAMLTree
{
	tree: Clade,
	mp: i32,
	lNL: String,
	tree_length: f64,
	kappa: f64, //ϰ: f64,
	omega: f64, //ω: f64,
	branch_data: HashMap<Endpoints, BranchData>
}


#[derive(Serialize, Deserialize, Debug, Clone)]
#[allow(non_snake_case)]
/// Associated data for a sequence in a PAML file, with selection data relative to its ancestor.
pub struct PAMLRecord { 
	pub name: String, 
	pub id: u32,
	pub ancestor: u32,
	pub descendents: HashSet<u32>,
	pub S: f64,
	pub dS: f64,
	dNdS: f64
}

#[allow(non_snake_case)]
impl PAMLRecord {
	/// Identifies whether this is a PAML-generated ancestor record.
	pub fn is_paml_generated(&self) -> bool {
		(self.name == String::new()) || (self.name == self.id.to_string())
	}

	/// Gets dNdS for the record.
	pub fn dNdS(&self) -> f64 { self.dNdS }

	/// Creates a blank PAML Record with Default Values.
	pub fn bare() -> PAMLRecord {
		PAMLRecord {
			name: String::new(),
			id: 0,
			ancestor: 0,
			descendents: HashSet::new(),
			S: 0.000,
			dS: -1.000,
			dNdS: -1.000
		}
	}
}

/// Record for handling the Phylogenetic Analysis by Maximum Likelihood tool.
#[derive(Debug)]
pub struct PAML 
{
	parameters: Vec<String>,
	pub records: Vec<PAMLRecord>,
	pub sections: Vec<PAMLSection>
}


impl PAML {
	/// Creates an empty PAML Record.
	pub fn bare() -> PAML {
		PAML {
			parameters: Vec::new(),
			records: Vec::new(),
			sections: Vec::new()
		}
	}

	/// Will generate identified branch records from parsed PAML Sections.
	pub fn regenerate_records(&mut self) {
		let mut records: Vec<PAMLRecord> = Vec::new();
		//println!("Section Count:{}", self.sections.len());
		for section in &self.sections {
			match section {/*
				PAMLSection::MarginalReconstruction(alignment) => {
					println!("In Reconstruction.");
					if records.len() < alignment.members.len() {
						records.resize(alignment.members.len(), PAMLRecord::bare());
					}
					for i in 0..alignment.members.len() {
						println!("NameSet:{}", alignment.members[i].id.clone());
						records[i].name = alignment.members[i].id.clone()
					}
				},*/
				PAMLSection::ChangeSummary(changesets) => {
					for (key, value) in changesets {
						let index = key.1 as usize - 1;

						if records.len() < key.1 as usize {
							records.resize(key.1 as usize, PAMLRecord::bare());
						}
						records[index].id = key.1;
						if records[index].name == "" {
							records[index].name = value.id.clone();
						}
						records[index].ancestor = key.0;
						if records.len() < key.0 as usize {
							records.resize(key.0 as usize, PAMLRecord::bare());
						}
						records[key.0 as usize - 1].id = key.0;
						records[key.0 as usize - 1].descendents.insert(key.1);
					}
				},
				PAMLSection::Branches(branches) => {
					for (key, value) in branches {
						let index = key.1 as usize - 1;
						if records.len() < key.1 as usize {
							records.resize(key.1 as usize, PAMLRecord::bare());
						}
						records[index].id = key.1;
						records[index].ancestor = key.0;
						records[index].S = value.S;
						records[index].dS = value.dS;
						records[index].dNdS = value.dNdS;
						if records.len() < key.0 as usize {
							records.resize(key.0 as usize, PAMLRecord::bare());
						}
						records[key.0 as usize - 1].id = key.0;
						records[key.0 as usize - 1].descendents.insert(key.1);
					}
				},
				_ => continue
			}
		}
		self.records = records
	}

	/// Appends a new file to the PAML data.
	pub fn append(&mut self, file_name: &str) -> Option<u32> {
		let mut counter: u32 = 0;
		match File::open(file_name) {
            Ok(file) => {
                let mut buf_reader = BufReader::new(file);
                loop {
					match PAMLSection::get_next_section(&mut buf_reader) {
						Ok(PAMLSection::Empty()) => break,
						Ok(section) => { counter += 1; self.sections.push(section) },
						Err(e) => {
							eprintln!("Error Reading Section: {}", e);
							continue
						}
					}
				};
				Some(counter)
            }
            Err(_e) => None
		}
	}

	/// Appends PDB Protein Structure data to the tree data held by PAML.
	pub fn append_prot_structure<T: ProteinData>(&mut self, prot: &T) {

		// Find and write alignment.
		Alignment::from(vec![prot.get_sequence()]).unwrap().write("pdb_a.fasta", "fasta").unwrap();
		if let Some(align) = self.get_ancestral_alignment() {
			align.write("temp_a.fasta", "fasta").unwrap();
		} else {
			// Awoooga
		}

		let mafft_run = Command::new("mafft-profile")
							.args(&["pdb_a.fasta", "temp_a.fasta"])
    						.stdout(Stdio::piped())
							.output()
							.expect("Failed to Execute Process");
		match Alignment::<Seq>::read_bytes(mafft_run.stdout, "fasta") {
			Ok(align) => {
				align.write("mafft.fasta", "fasta").unwrap();
				self.sections.insert(0, PAMLSection::AlignWithAncestral(align))
			}
			Err(e) => eprintln!("MAFFT Run Failed: {:?}", e)
		}
	}

	/// Fetchs the full alignment of original and ancestral sequences.
	pub fn get_ancestral_alignment(&self) -> Option<Alignment<Seq>> {
		for section in &self.sections {
			match section {
				PAMLSection::MarginalReconstruction(prob_align, _site_freqs) => return Some(prob_align.likely_alignment()),
				PAMLSection::AlignWithAncestral(align) => return Some((*align).clone()),
				_ => continue
			}
		}
		None
	}

	/// Reads a file into a new PAML section.
	pub fn read(file_name: &str) -> Result<PAML, io::Error> {
		let mut data: PAML = PAML::bare();

		match File::open(file_name) {
            Ok(file) => {
                let mut buf_reader = BufReader::new(file);
                loop {
					match PAMLSection::get_next_section(&mut buf_reader) {
						Ok(PAMLSection::Empty()) => break,
						Ok(section) => data.sections.push(section),
						Err(e) => return Err(e)
					}
				};
				Ok(data)
            }
            Err(e) => Err(e)
        }
	}
}
