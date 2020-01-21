//! Tools for handling phenotypic data.

use crate::clade::Clade;
use serde::{Serialize, Deserialize};
use std::collections::HashMap;

#[derive(Serialize, Deserialize)]
pub struct PhylogeneticSet {
	pub species: HashMap<String, PhylogeneticData>,
    pub species_list: Vec<String>
}



#[derive(Serialize, Deserialize)]
pub struct PhylogeneticData {
	pub files: HashMap<String, PhylogeneticFiles>,
	pub pathways: HashMap<String, Vec<String>>,
	pub genes: HashMap<String, Vec<String>>,
}

#[derive(Serialize, Deserialize)]
pub struct PhylogeneticFiles {
	pub tree: String,
	pub paml: String,
	pub pdb: String,
	pub gene_name: String
}

#[derive(Serialize, Deserialize)]
pub struct ProcessedData {
	pub selection: HashMap<String, Vec<f64>>,
	pub related_selection: HashMap<String, f64>,
	pub ancestor_pairs: HashMap<String, (u32, u32)>,
	//trees: HashMap<String, Vec<Clade>>,
	pub related_arctic_species:Vec<String>,
	pub related_nonarctic_species: Vec<String>,
	pub pathways: HashMap<String, Vec<String>>,
	pub genes: HashMap<String, Vec<String>>,
}

#[derive(Serialize, Deserialize)]
pub struct PhenotypicRecord {
    pub species: String,
    pub winter_coat: String,
    pub summer_coat: String,
    pub concealing_coat: String,
    pub seasonal_coat: String,
    pub food_hoarding: String,
    pub body_mass: f64,
    pub body_mass_seasonal: String,
    _ignore: String
}

impl ProcessedData {
	pub fn from(pathways: HashMap<String, Vec<String>>, genes: HashMap<String, Vec<String>>) -> ProcessedData {
		ProcessedData {
			selection: HashMap::new(),
			related_selection: HashMap::new(),
			ancestor_pairs: HashMap::new(),
			//trees: HashMap::new(),
			related_arctic_species: Vec::new(),
			related_nonarctic_species: Vec::new(),
            pathways,
            genes
		}
	}
}

/// Gets all entries in the phylogenetic tree that mach a name in the species list.
/// 
/// # Arguments
/// * 'tree' - Reference to root node of tree to search.
/// * 'name_check' - Reference to vector of names to search
/// 
/// # TODO
/// 
/// Switch vector to generic slice?
pub fn get_matching_entries(tree: &Clade, name_check: &[String]) -> Vec<String> {
	let mut names: Vec<String> = Vec::new();
	for child in tree.children.iter() {
        names.append(&mut get_matching_entries(child, name_check));
	}
    let node_name = tree.name.split('#').last().unwrap().to_lowercase();
	if name_check.contains(&node_name) {
        names.push(node_name);
	}
	names
}

/// Gets all entries in the phylogenetic tree that don't mach a name in the species list.
/// 
/// # Arguments
/// * 'tree' - Reference to root node of tree to search.
/// * 'name_check' - Reference to vector of names to search
/// 
/// # TODO
/// 
/// Switch vector to generic slice?
pub fn get_divergent_entries(tree: &Clade, name_check: &[String]) -> Vec<String> {
	let mut names: Vec<String> = Vec::new();
	for child in tree.children.iter() {
		names.append(&mut get_divergent_entries(child, name_check));
	}

    let node_name = tree.name.split('#').last().unwrap().to_lowercase();
	if !name_check.contains(&node_name) {
		names.push(node_name);
	}
	names
}