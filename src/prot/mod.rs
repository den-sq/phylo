//! Tools for the analysis of Protein Databank Information
use crate::seq::{Seq, SeqData};
use crate::alignment::Alignment;
use crate::paml::PAMLRecord;
use serde::{Serialize, Deserialize};

/// Holds simple protein structure data from a parsed DSSP file.
#[derive(Serialize, Deserialize, Debug, Clone, Default)]
pub struct ProteinStructure {
	pub residues: Vec<char>,
	pub ss_simple: Vec<char>,
	solvent_accessibility: Vec<String>
}

impl ProteinStructure {
	pub fn new() -> ProteinStructure {
		ProteinStructure {
			residues: Vec::new(),
			ss_simple: Vec::new(),
			solvent_accessibility: Vec::new()
		}
	}
}

pub trait ProteinData {
	fn get_solvent_accessibility (&self, pos: usize) -> Option<f64>;
	fn get_sequence(&self) -> Seq;
	fn missing_residue_check(&self) -> bool;
}

impl ProteinData for ProteinStructure {
	/// Fetches the solvent accessibility at a given position.
	fn get_solvent_accessibility (&self, pos: usize) -> Option<f64> {
		if pos < self.solvent_accessibility.len() {
			match self.solvent_accessibility[pos].parse::<f64>() {
				Ok(value) => Some(value),
				Err(_) => None
			}
		} else { None }
	}

	/// Gets the amino acid sequence as a Seq object.
	fn get_sequence(&self) -> Seq {
		Seq::new(String::from("pdb"), self.residues.clone())
	}

	/// Checks if the residue at a given position is missing.
	fn missing_residue_check(&self) -> bool {
		for c in &self.residues {
			if *c == '-' { return true }
		}
		false
	}
}


#[derive(Serialize, Deserialize, Debug, Clone, Default)]
pub struct ExtendedProteinStructure {
	pub base: ProteinStructure,
	families: Vec<(Vec<PAMLRecord>, Alignment<Seq>, String)>,
	pdb_id: String,
}

impl ExtendedProteinStructure {
	pub fn extend(prot: ProteinStructure, pdb_id: String) -> ExtendedProteinStructure {
		ExtendedProteinStructure {
			base: prot,
			families: Vec::new(),
			pdb_id,
		}
	}
	pub fn new() -> ExtendedProteinStructure {
		ExtendedProteinStructure {
			base: ProteinStructure::new(),
			families: Vec::new(),
			pdb_id: String::new(),	
		}
	}

	pub fn add_family(&mut self, record_set: Vec<PAMLRecord>, align: Alignment<Seq>, from: String) {
		self.families.push((record_set,align, from));
	}

	pub fn get_families(&self) -> &Vec<(Vec<PAMLRecord>, Alignment<Seq>, String)> {
		&self.families
	}
}


impl SeqData for ExtendedProteinStructure {
	fn len(&self) -> usize { self.base.residues.len() }
	fn id(&self) -> &str { &self.pdb_id }
	fn seq(&self) -> Vec<char> { Seq::new(String::from("pdb"), self.base.residues.clone()).seq() }
	fn is_empty(&self) -> bool { self.base.residues.is_empty() }
}


impl ProteinData for ExtendedProteinStructure {
	/// Fetches the solvent accessibility at a given position.
	fn get_solvent_accessibility (&self, pos: usize) -> Option<f64> {
		self.base.get_solvent_accessibility(pos)
	}

	/// Gets the amino acid sequence as a Seq object.
	fn get_sequence(&self) -> Seq {
		Seq::new(String::from("pdb"), self.base.residues.clone())
	}

	/// Checks if the residue at a given position is missing.
	fn missing_residue_check(&self) -> bool {
		for c in &self.base.residues {
			if *c == '-' { return true }
		}
		false
	}
}