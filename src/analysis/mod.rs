extern crate rand;
extern crate num;
extern crate indexmap;

use std::collections::{HashMap, HashSet};
use std::hash::Hash;
use std::fmt;
use std::iter::{Iterator, FromIterator};
use crate::seq::SeqData;
use crate::prot::{ExtendedProteinStructure, ProteinData};
use crate::utils::{arithmatic_mean};
use self::rand::prelude::*;
use self::rand::distributions::Uniform;
use self::indexmap::IndexMap;
use serde::{Serialize, Deserialize};

pub type StatGroupExport = HashMap<String, HashMap<String, (f64, f64, f64)>>;
pub type CountGroupExport = HashMap<String, HashMap<String, Vec<f64>>>;

#[derive(Serialize, Deserialize, Debug, Clone)]
/// A simple statistical value with supporting information (name, standard deviation, and p-value).
pub struct Statistic {
	field: String,
	pub value: f64,
	std_dev: f64,
	p: f64
}

impl Statistic {
	/// Create a new statistic, given name, value, and stddev.
	pub fn new(field: String, value: f64, std_dev: f64) -> Statistic {
		Statistic {
			field,
			value,
			std_dev,
			p: 0.00
		}
	}

	/// Computes the distance between two statistics (absolute value of difference).
	pub fn distance(&self, other: &Statistic) -> f64 { (self.value - other.value).abs() }

}

impl fmt::Display for Statistic {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
		write!(f, " \"{}\": \"{} +/- {}\"", self.field, self.value, self.std_dev)
	}
}

#[derive(Serialize, Deserialize, Debug, Clone)]
pub struct StatGroup {
	pub name: String,
	pub stats: HashMap<String, Statistic>,
	n: u64
}

impl StatGroup {
	pub fn new(name: String, n: u64) -> StatGroup {
		StatGroup {
			name,
			n,
			stats: HashMap::new(),
		}
	}

	pub fn export(&self, format_map: &HashMap<char, String>) -> StatGroupExport {
		let mut formatted_stats: HashMap<String,  HashMap<String, (f64, f64, f64)>> = HashMap::new();
		for stat in self.stats.values() {
			let first = stat.field.clone().chars().next().unwrap();
			let extra = if stat.field.len() > 1 { stat.field[1..stat.field.len()].to_string() } else { "-".to_string() };
			//eprintln!("{}:{}", first, extra);
			formatted_stats.entry(format_map.get(&first).unwrap().clone()).or_insert_with(HashMap::new)
				.insert(extra, (stat.value, stat.std_dev, stat.p));
			//formatted_stats.insert(format_map.get(&stat.field.clone().chars().next().unwrap()).unwrap().clone(), format!("{:.3} +/- {:.3}", stat.value, stat.std_dev));
		}
		formatted_stats
	}

}

impl fmt::Display for StatGroup {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
		let mut contents = String::new();
		
		for stat in self.stats.values() {
			contents.push_str(&format!("{}: {:.3} +/- {:.3}", stat.field, stat.value, stat.std_dev));
		}
		let x = contents.len(); contents.truncate(x - 1);
		write!(f, "\"{}\": {{ {} }}", self.name, contents)
	}
}


// Secondary Structure and Solvent Accessibility Pair.
#[derive(Serialize, Deserialize, Debug, Eq, PartialEq, Hash, Clone, Copy)]
pub struct ResFeature { pub ss: char, pub sa: char }

impl fmt::Display for ResFeature {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
		write!(f, "{}{}", self.ss, self.sa)
	}
}

impl From<String> for ResFeature {
	fn from(source: String) -> Self {
		let mut char_list = source.chars();
		ResFeature {
			ss: char_list.next().unwrap(),
			sa: char_list.next().unwrap()
		}
	}
}

#[derive(Serialize, Deserialize, Debug, Default)]
#[allow(non_snake_case)]
pub struct ProteinAnalytics {
	pub proteins: Vec<ExtendedProteinStructure>,
	pub protein_stats: Vec<HashMap<String, StatGroup>>,
	pub stats: HashMap<String, StatGroup>,
	pub bad_datasources: HashMap<String, HashSet<String>>,	
	pub dNdS_distribution: Vec<(f64, f64, f64)>,
	protein_counts: Vec<Vec<HashMap<String, HashMap<char, u32>>>>,
	dist_protein_counts: Vec<Vec<HashMap<String, HashMap<ResFeature, u32>>>>,
	net_counts: HashMap<String, HashMap<char, u32>>,
	dist_net_counts: HashMap<String, HashMap<ResFeature, u32>>,
	bad_counts: HashMap<String, u64>,
	wilke: HashMap<String, f64>,
	aa_1_to_3: HashMap<char, String>,
	format_map: HashMap<char, String>
}



impl fmt::Display for ProteinAnalytics {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        //write!(f, "({}, {})", self.x, self.y)
		
		for (name, stat) in &self.stats {
			let name_vec = ProteinAnalytics::long_stat_categorization(name);
			match write!(f, "{} ({}) ({}): {}", name_vec[0], name_vec[1], name_vec[2], stat) {
				Ok(_) => (),
				Err(e) => return Err(e)
			}
		}
		Ok(())
    }
}

impl ProteinAnalytics {
	fn long_stat_categorization(short: &str) -> Vec<String> {
		let mut long_cat: Vec<String> = Vec::new();
		let name_vec: Vec<char> = short.chars().collect();

		long_cat.push(match name_vec[0] {
			'l' => "Changed on a Lineage".to_string(),
			'p' => "Changed Per-Protein".to_string(), _ => String::new()
		});
		long_cat.push(match name_vec[1] {
			'p' => "Positively Selected Lineages".to_string(),
			'n' => "Negatively Selected Lineages".to_string(), _ => String::new()
		});
		long_cat.push(match name_vec[2] {
			'c' => "Substituted Sites".to_string(),
			'n' => "Non-Substituted Sites".to_string(), 
			'a' => "All Sites".to_string(), 
			'g' => "Randomly Allocated Sites".to_string(), 
			'r' => "Rate of Substitution".to_string(), _ => String::new()
		});
		long_cat.push(if name_vec.len() > 3 {
			match name_vec[3] {
				'e' => "(Split)".to_string(),
				_ => String::new()
			}
		} else {
			"(All)".to_string()
		});

		long_cat
	}

	pub fn export_raw_counts(&mut self) -> HashMap<String, HashMap<String, CountGroupExport>> {
		let mut count_hash: HashMap<String, HashMap<String, CountGroupExport>> = HashMap::new();

		for prot in self.protein_counts.iter_mut() {
			for prot_count in prot.iter_mut() {
				for (stat_type, stats) in prot_count {
					let stat_group = ProteinAnalytics::long_stat_categorization(stat_type);
					let hash = count_hash.entry(stat_group[0].clone()).or_insert_with(HashMap::new)
						.entry(stat_group[1].clone()).or_insert_with(HashMap::new)
						.entry(stat_group[2].clone()).or_insert_with(HashMap::new);

					// Calc total; if 0 ignore this set.
					let total_residues = stats[&'F'] + stats[&'P'];
					if total_residues == 0 { continue };

					for (stat, name) in &self.format_map {
						let value = stats.entry(*stat).or_insert(0);
						hash.entry(name.clone())
							.and_modify(|x| { x.push(*value as f64 / total_residues as f64) })
							.or_insert_with(|| vec![*value as f64 / total_residues as f64]);
					}
				}
			}
		}

		count_hash
	}

	// Exports accumulated statistical data.
	pub fn export_stats(&self) -> HashMap<String, HashMap<String, HashMap<String, StatGroupExport>>> {
		let mut stat_hash: HashMap<String, HashMap<String, HashMap<String, StatGroupExport>>> = HashMap::new();
		for (name, stat_group) in &self.stats {
			let cats = ProteinAnalytics::long_stat_categorization(name);
			let data = stat_group.export(&self.format_map);

			stat_hash.entry(cats[0].clone()).or_insert_with(HashMap::new)
				.entry(cats[1].clone()).or_insert_with(HashMap::new)
				.entry(cats[2].clone()).and_modify(|x| { for (key, value) in data.clone() { x.insert(key, value); } })
				.or_insert(data);
		}
		stat_hash
	}

	// Exports Relevant Counts of Values.
	pub fn export_counts(&self) -> HashMap<String, u64> {
		let mut issue_hash: HashMap<String, u64> = self.bad_counts.clone();
		for (key, records) in &self.bad_datasources {
			issue_hash.insert(key.clone(), records.len() as u64);
		}
		let residue_total =
			*issue_hash.entry("Skipped Residues".to_string()).or_insert(0) + 
			*issue_hash.entry("Unidentified Residues".to_string()).or_insert(0) + 
			*issue_hash.entry("Compared Residues".to_string()).or_insert(0);
		issue_hash.insert("Total Residues Scanned".to_string(), residue_total);
		let residue_total =
			*issue_hash.entry("Per-Protein Skipped Residues".to_string()).or_insert(0) + 
			*issue_hash.entry("Per-Protein Unidentified Residues".to_string()).or_insert(0) + 
			*issue_hash.entry("Per-Protein Compared Residues".to_string()).or_insert(0);
		issue_hash.insert("Per-Protein Total Residues Scanned".to_string(), residue_total);
		issue_hash.insert("Protein Count".to_string(), self.proteins.len() as u64);
		issue_hash.insert("Branch Count".to_string(), self.stats["lpne"].n);
		issue_hash
	}

	// Creates a new Protein Analytics object, with default values.
	pub fn new() -> ProteinAnalytics {
		ProteinAnalytics {
			proteins: Vec::new(),
			wilke: [(String::from("ALA"), 129.0), (String::from("ARG"), 274.0), (String::from("ASN"), 195.0), 
					(String::from("ASP"), 193.0), (String::from("CYS"), 167.0), (String::from("GLN"), 225.0), 
					(String::from("GLU"), 223.0), (String::from("GLY"), 104.0), (String::from("HIS"), 224.0), 
					(String::from("ILE"), 197.0), (String::from("LEU"), 201.0), (String::from("LYS"), 236.0), 
					(String::from("MET"), 224.0), (String::from("PHE"), 240.0), (String::from("PRO"), 159.0),
					(String::from("SER"), 155.0), (String::from("THR"), 172.0), (String::from("TRP"), 285.0), 
					(String::from("TYR"), 263.0), (String::from("VAL"), 174.0)].iter().cloned().collect(),
			aa_1_to_3: [('A', String::from("ALA")), ('C', String::from("CYS")), ('D', String::from("ASP")), 
						('E', String::from("GLU")), ('F', String::from("PHE")), ('G', String::from("GLY")),
						('H', String::from("HIS")), ('I', String::from("ILE")), ('K', String::from("LYS")), 
						('L', String::from("LEU")), ('M', String::from("MET")), ('N', String::from("ASN")),
						('P', String::from("PRO")), ('Q', String::from("GLN")), ('R', String::from("ARG")),
						('S', String::from("SER")), ('T', String::from("THR")), ('V', String::from("VAL")), 
						('W', String::from("TRP")), ('Y', String::from("TYR"))].iter().cloned().collect(),
			format_map: [('A', "Helix".to_string()), ('*', "All Sites".to_string()),
						('G', "3₁₀ Helix".to_string()), ('H', "α-Helix".to_string()), ('I', "π-Helix".to_string()), 
						('T', "Turn".to_string()), ('S', "Bend".to_string()), ('C', "Coil".to_string()), 
						('B', "β-Bridge".to_string()), ('E', "β-Sheet".to_string()), (' ', "Unidentified".to_string()), 
						('P', "Buried".to_string()), ('F', "Exposed".to_string())].iter().cloned().collect(),
			dNdS_distribution: Vec::new(),
			bad_counts: HashMap::new(),
			bad_datasources: HashMap::new(),
			protein_counts: Vec::new(),
			dist_protein_counts: Vec::new(),
			net_counts: HashMap::new(),
			dist_net_counts: HashMap::new(),
			protein_stats: Vec::new(),
			stats: HashMap::new()
		}
	}

	// Extends analytics by adding a new protein.
	pub fn extend(&mut self, ps: ExtendedProteinStructure) {
		self.proteins.push(ps)
	}

}

impl ProteinAnalytics {

	// Calculates % of how residues are distributed amongst features (secondary structure / solvent accessibility)
	fn distribution_calc<T: Eq + PartialEq + Hash + fmt::Display>
			(total: f64, counts: &mut HashMap<T, u32>, name: String, n: u64) -> StatGroup {
		let mut stat_group = StatGroup::new(name, n);
		for (stat, count) in counts {
			stat_group.stats.insert(stat.to_string(), Statistic::new(stat.to_string(), *count as f64 / total, 0.00));
		}
		stat_group
	}

	// Calculates the rate of substitution.
	fn rate_calc<T: Eq + PartialEq + Hash + fmt::Display>
			(nosub: &HashMap<T, u32>, sub: &HashMap<T, u32>, name: String, n: u64) -> StatGroup {
		let mut stat_group = StatGroup::new(name, n);
		for (stat, count) in sub {

			let total = match nosub.get(stat) { Some(x) => *x, None => 0 } + *count;
			if total == 0 { continue; }

			stat_group.stats.insert(stat.to_string(), Statistic::new(stat.to_string(), 
				*count as f64 / total as f64, 0.00));
		}
		for (stat, count) in nosub {
			if *count == 0 { continue; }
			if !sub.contains_key(stat) {
				stat_group.stats.insert(stat.to_string(), Statistic::new(stat.to_string(), 
					0.00, 0.00));
			}
		}
		stat_group
	}

	#[allow(dead_code)]
	fn rate_calc_pretotaled<T: Eq + PartialEq + Hash + fmt::Display>
			(total: &HashMap<T, u32>, sub: &HashMap<T, u32>, name: String, n: u64) -> StatGroup {
		let mut stat_group = StatGroup::new(name, n);
		for (stat, count) in total {
			if *count == 0 { continue; }
			if !sub.contains_key(stat) {
				stat_group.stats.insert(stat.to_string(), Statistic::new(stat.to_string(), 
					0.00, 0.00));
			} else {
				stat_group.stats.insert(stat.to_string(), Statistic::new(stat.to_string(), 
					sub[stat] as f64 / *count as f64, 0.00));
			}
		}
		stat_group
	}

	fn rate_calc_for_hash<T: Eq + PartialEq + Hash + fmt::Display + Copy>
			(total: &HashMap<T, u32>, sub: &mut HashMap<T, u32>) -> HashMap<T, f64> {
		let mut stat_group: HashMap<T, f64> = HashMap::new();
		for (stat, count) in total {
			if *count == 0 { continue; }
			stat_group.insert(*stat, *sub.entry(*stat).or_insert(0) as f64 / *count as f64);
		}
		stat_group
	}

	pub fn generate_bootstrap(totals: HashMap<String, IndexMap<ResFeature, u64>>, 
			sample_size: HashMap<String, u64>, sample_range: HashMap<String, u64>, 
			num_groups: u32, bases: Vec<String>) -> HashMap<String, Vec<StatGroup>> {
		let mut bootstrap_data: HashMap<String, Vec<StatGroup>> = HashMap::new();

		for base in bases.iter() {
			
			let mut generated_vals: Vec<StatGroup> = Vec::new();
	
			let sample = Uniform::from(0..sample_range[base]);
			let mut sampler = rand::thread_rng();
			println!("{:?}|{}|{}", totals[base], sample_range[base], sample_size[base]);
			for _i in 0..num_groups {
				let mut random_sample = HashMap::new();
				for _j in 0..sample_size[base] {
					let res_num =  sample.sample(&mut sampler);
					for (key, value) in totals[base].iter() {
						if res_num < *value {
							random_sample.entry(*key).and_modify(|x| { *x += 1;}).or_insert(1);
							break;
						}
					};
				}
				let h = *random_sample.entry(ResFeature { ss: 'G', sa: 'F' }).or_insert(0) + 
							*random_sample.entry(ResFeature { ss: 'H', sa: 'F' }).or_insert(0) + 
							*random_sample.entry(ResFeature { ss: 'I', sa: 'F' }).or_insert(0);
				random_sample.insert(ResFeature { ss: 'A', sa: 'F' }, h);
				let h = *random_sample.entry(ResFeature { ss: 'G', sa: 'P' }).or_insert(0) + 
							*random_sample.entry(ResFeature { ss: 'H', sa: 'P' }).or_insert(0) + 
							*random_sample.entry(ResFeature { ss: 'I', sa: 'P' }).or_insert(0);
				random_sample.insert(ResFeature { ss: 'A', sa: 'P' }, h);

				// Calculate Distributions.
				let temp = ProteinAnalytics::distribution_calc(
					sample_size[base] as f64, &mut random_sample, base.to_string(), 0);
				//let mut x = 0.0;
				//for (_stat_name, stat) in temp.stats.iter() {
				//	x += stat.value;
				//}
				//eprintln!("{}:{}:{}", x, temp.stats[&ResFeature { ss: 'A', sa: 'P' }.to_string()], temp.stats[&ResFeature { ss: 'A', sa: 'F' }.to_string()]);
				generated_vals.push(temp);
			}
			bootstrap_data.insert(base.clone(), generated_vals);
		}
		bootstrap_data
	}

	pub fn generate_sample(&mut self) {
		let bases = ["lp", "ln", "pp", "pn"];
		let structure = ['G', 'H', 'I', 'E', 'B', 'T', 'S', ' '];
		let num_groups = 20000;
		//for (stat_type, stats) in &self.net_counts {
		for base in bases.iter() {
			let mut totals: IndexMap<ResFeature, u32> = IndexMap::from_iter(self.dist_net_counts[&format!("{}n", base)].clone());
			let mut summed: HashMap<ResFeature, u32> = self.dist_net_counts[&format!("{}n", base)].clone();
			let mut generated_samples: Vec<HashMap<ResFeature, u32>> = Vec::new();
			let mut generated_vals: Vec<StatGroup> = Vec::new();
			let mut generated_rate_of_sub: Vec<HashMap<ResFeature, f64>> = Vec::new();
			//println!("{:?}\n{:?}\n{:?}", totals, self.dist_net_counts[&format!("{}n", base)], self.dist_net_counts[&format!("{}c", base)]);


			// Sample size is # of substituted residues, so we get the sum minus the (duplicate) helixes.
			let sample_size: u32 = self.dist_net_counts[&format!("{}c", base)].values().sum::<u32>() - 
								self.dist_net_counts[&format!("{}c", base)][&ResFeature { ss: 'A', sa: 'F' }] -
								self.dist_net_counts[&format!("{}c", base)][&ResFeature { ss: 'A', sa: 'P' }];


			// Get rid of Helixes.  Get Total.
			totals.remove(&ResFeature { ss: 'A', sa: 'F' });
			totals.remove(&ResFeature { ss: 'A', sa: 'P' });
			let sample_range: u32 = totals.values().sum::<u32>() + sample_size;

			// We add up the incremented totals for use with the random number generator.
			let mut inc: u32 = 0;
			for (key, count) in totals.iter_mut() {
				*count += inc + self.dist_net_counts[&format!("{}c", base)][key];
				inc = *count;
				summed.entry(*key).and_modify(|x| { *x += self.dist_net_counts[&format!("{}c", base)][key]; }).or_insert(0);
			}
								
			let sample = Uniform::from(0..sample_range);
			let mut sampler = rand::thread_rng();
			println!("{:?}|{}|{}", totals, sample_range, sample_size);
			for _i in 0..num_groups {
				let mut random_sample = HashMap::new();
				for _j in 0..sample_size {
					let res_num =  sample.sample(&mut sampler);
					for (key, value) in totals.iter() {
						if res_num < *value {
							random_sample.entry(*key).and_modify(|x| { *x += 1;}).or_insert(1);
							break;
						}
					};
				}
				let h = *random_sample.entry(ResFeature { ss: 'G', sa: 'F' }).or_insert(0) + 
							*random_sample.entry(ResFeature { ss: 'H', sa: 'F' }).or_insert(0) + 
							*random_sample.entry(ResFeature { ss: 'I', sa: 'F' }).or_insert(0);
				random_sample.insert(ResFeature { ss: 'A', sa: 'F' }, h);
				let h = *random_sample.entry(ResFeature { ss: 'G', sa: 'P' }).or_insert(0) + 
							*random_sample.entry(ResFeature { ss: 'H', sa: 'P' }).or_insert(0) + 
							*random_sample.entry(ResFeature { ss: 'I', sa: 'P' }).or_insert(0);
				random_sample.insert(ResFeature { ss: 'A', sa: 'P' }, h);


				// Calculate Rate of Substitution & Insert
				// Doing this first fills in games in random_sample.
				generated_rate_of_sub.push(ProteinAnalytics::rate_calc_for_hash(&summed, &mut random_sample));

				// Calculate Distributions.
				generated_vals.push(ProteinAnalytics::distribution_calc(
					sample_size as f64, &mut random_sample, format!("{}g", base), 0));
				
				//println!("{:?}", random_sample);

				generated_samples.push(random_sample);
			}


			// Add in the All-Cases before distance calculation.
			for ss in ['A', 'G', 'H', 'I', 'E', 'B', 'T', 'S', ' '].iter() {
				let ss_all = ResFeature { ss: *ss, sa: '*' };
				let exposed = ResFeature { ss: *ss, sa: 'F' };
				let buried = ResFeature { ss: *ss, sa: 'P' };
				
				// self.stats.get(&format!("{}ce", base)) (Base Residue Sub Distribution)
				for sub_group in ["ce", "ne", "ae"].iter() {
					if let Some(stat_group) = self.stats.get_mut(&format!("{}{}", base, sub_group)) {
						let total: f64 = stat_group.stats[&exposed.to_string()].value + stat_group.stats[&buried.to_string()].value;
						stat_group.stats.insert(ss_all.to_string(), Statistic::new(ss_all.to_string(), total, 0.00));
					}
				}

				// self.stats[&format!("{}re", base)] (Base Rate of Substitution)
				let mut numerator: u32 = 0; let mut denominator: u32 = 0;
				for sa in ['F', 'P'].iter() {
					let res = ResFeature { ss: *ss, sa: *sa };
					numerator += self.dist_net_counts[&format!("{}c", base)][&res];
					denominator += self.dist_net_counts[&format!("{}n", base)][&res];
				}						
				if let Some(stat_group) = self.stats.get_mut(&format!("{}re", base)) {
					stat_group.stats.insert(ss_all.to_string(), Statistic::new(ss_all.to_string(), numerator as f64 / (numerator + denominator) as f64, 0.00));
				}
				
				// generated_vals (Random Residue Sub Distribution)
				for sub_dist in generated_vals.iter_mut() {
					//eprintln!("{:?}:{}:{}:{}", sub_dist.stats.keys(), ss_all, exposed, buried);
					let total: f64 = sub_dist.stats[&exposed.to_string()].value + sub_dist.stats[&buried.to_string()].value;
					sub_dist.stats.insert(ss_all.to_string(), Statistic::new(ss_all.to_string(), total, 0.00));
				}

				// generated_rate_of_sub (Random Rate of Substition)
				for i in 0..generated_rate_of_sub.len() {
					let mut numerator: u32 = 0; let mut denominator: u32 = 0;
					for sa in ['F', 'P'].iter() {
						let res = ResFeature { ss: *ss, sa: *sa };
						numerator += generated_samples[i][&res];
						denominator += self.dist_net_counts[&format!("{}n", base)][&res];
					}
					generated_rate_of_sub[i].insert(ss_all, numerator as f64 / (numerator + denominator) as f64);
				}
			}
			for sa in ['F', 'P'].iter() {
				let sa_all = ResFeature { ss: '*', sa: *sa };

				// self.stats.get(&format!("{}ce", base)) (Base Residue Sub Distribution)
				for sub_group in ["ce", "ne", "ae"].iter() {
					if let Some(stat_group) = self.stats.get_mut(&format!("{}{}", base, sub_group)) {
						let mut total: f64 = 0.00;
						for ss in structure.iter() {
							total += stat_group.stats[&format!("{}{}", ss, sa)].value;
						}
						stat_group.stats.insert(sa_all.to_string(), Statistic::new(sa_all.to_string(), total, 0.00));
					}
				}

				// self.stats[&format!("{}re", base)] (Base Rate of Substitution)
				let mut substitutions: u32 = 0; let mut fixed: u32 = 0;
				for ss in structure.iter() {
					let res = ResFeature { ss: *ss, sa: *sa };
					substitutions += self.dist_net_counts[&format!("{}c", base)][&res];
					fixed += self.dist_net_counts[&format!("{}n", base)][&res];
				}
						
				if let Some(stat_group) = self.stats.get_mut(&format!("{}re", base)) {
					stat_group.stats.insert(sa_all.to_string(), Statistic::new(sa_all.to_string(), substitutions as f64 / (fixed + substitutions) as f64, 0.00));
				}

				// generated_vals (Random Residue Sub Distribution)
				for sub_dist in generated_vals.iter_mut() {
					let mut total: f64 = 0.00;
					for ss in structure.iter() {
						total += sub_dist.stats[&format!("{}{}", ss, sa)].value;
					}
					sub_dist.stats.insert(sa_all.to_string(), Statistic::new(sa_all.to_string(), total, 0.00));
				}

				// generated_rate_of_sub (Random Rate of Substition)
				for i in 0..generated_rate_of_sub.len() {
					let mut substitutions: u32 = 0; let mut denominator: u32 = 0;
					for ss in structure.iter() {
						let res = ResFeature { ss: *ss, sa: *sa };
						substitutions += generated_samples[i][&res];
						denominator += self.dist_net_counts[&format!("{}n", base)][&res];
					}
					generated_rate_of_sub[i].insert(sa_all, substitutions as f64 / (denominator + substitutions) as f64);
				}
			}

			// Calculate distance between all residue distribution and positively selected residues for distribution.
			let mut gaps: HashMap<String, f64> = HashMap::new();
			let bare_stat = Statistic::new("bare".to_string(), 0.00, 0.00);
			{ 	// Braces so borrow is dropped at the end of them and we can do a mutable borrow later.
				for (name, stat) in &self.stats[&format!("{}ae", base)].stats {
					gaps.insert(name.clone(), stat.distance(
						if self.stats.contains_key(&format!("{}ce", base)) {
							&self.stats.get(&format!("{}ce", base)).unwrap().stats.get(&name.clone()).unwrap()
						} else { &bare_stat }
					));
				}
			}

			// Calculate Rate of Substituion Distances
			let mut rate_distance: HashMap<ResFeature, f64> = HashMap::new();
			let mean: f64;
			{	
				let mut ss_vals: Vec<f64> = Vec::new();
				
				for (name, stat) in &self.stats[&format!("{}re", base)].stats {
					if name.find('*') != None {
						ss_vals.push(stat.value);
					}
				}
				mean = self::arithmatic_mean(ss_vals.iter());
				for (name, stat) in &self.stats[&format!("{}re", base)].stats {
					rate_distance.insert(ResFeature::from(name.clone()), (stat.value - mean).abs());
				}
			}
			eprintln!("{}", mean);
			eprintln!("{:?}", self.stats.keys());
			eprintln!("{}", base);

			// Calculate distance from all-residues in random samples.
			let mut random_gaps: Vec<HashMap<String, f64>> = Vec::new();
			for stat_group in generated_vals {
				let mut random_gap_set: HashMap<String, f64> = HashMap::new();
				for (name, stat) in stat_group.stats {
					random_gap_set.insert(name.clone(), 
						stat.distance(&self.stats.get(&format!("{}ae", base)).unwrap()
							.stats.get(&name.clone()).unwrap()));
				}
				random_gaps.push(random_gap_set);
			}

			// As Same for Rate of Sub
			let mut random_rate_distance: Vec<HashMap<ResFeature, f64>> = Vec::new();
			for stat_group in generated_rate_of_sub {
				let mut distances: HashMap<ResFeature, f64> = HashMap::new();
				for (name, result) in stat_group.clone() {
					distances.insert(name, (result - mean).abs());
				}
				//eprintln!("{:?}\n{:?}", distances, stat_group);

				random_rate_distance.push(distances);
			}

			// Mark where the distances are greater
			let mut gap_counts: HashMap<String, u64> = HashMap::new();
			for random_gap_set in random_gaps {
				for (name, random_distance) in random_gap_set {
					//println!("{}:{}:{}", name, random_distance, gaps[name.as_str()]);
					gap_counts.entry(name.clone())
						.and_modify(|x| { *x += (random_distance > gaps[name.as_str()]) as u64 })
						.or_insert((random_distance > gaps[name.as_str()]) as u64);
				}
			}

			// Same for Rate of Sub
			let mut rate_distance_counts: HashMap<ResFeature, u64> = HashMap::new();
			for distances in random_rate_distance {
				for (name, random_distance) in distances {
					//println!("{}:{}:{}", name, random_distance, rate_distance[&name]);
					rate_distance_counts.entry(name)
						.and_modify(|x| { *x += (random_distance > rate_distance[&name]) as u64 })
						.or_insert((random_distance > rate_distance[&name]) as u64);
				}
			}

			eprintln!("{:?}", gap_counts);
			eprintln!("{:?}", rate_distance_counts);

			// Generate & assign P values.
			for (name, stat) in self.stats.get_mut(&format!("{}ce", base)).unwrap().stats.iter_mut() {
				stat.p = *gap_counts.entry(name.clone()).or_insert(0) as f64 / num_groups as f64;
			}
			for (name, stat) in self.stats.get_mut(&format!("{}re", base)).unwrap().stats.iter_mut() {
				stat.p = *rate_distance_counts.entry(ResFeature::from(name.clone())).or_insert(0) as f64 / num_groups as f64;
			}
		}
	}

	#[allow(non_snake_case)]
	pub fn dNdS_record(&mut self) {
		let mut dNdS_v_dS: Vec<(f64, f64, f64)> = Vec::new();
		for prot in &self.proteins {
			for (records, _align, _taed_id) in prot.get_families() {
				for record in records {
					if record.is_paml_generated() && record.dNdS() >= 1.00 {
						dNdS_v_dS.push((record.dNdS(), record.dS, record.S));
					}
				}
			}
		}
		self.dNdS_distribution = dNdS_v_dS;
	}


	#[allow(illegal_floating_point_literal_pattern)]
	pub fn total(&mut self) {
		fn increment<T: Hash + Eq + PartialEq + Clone>(count: &mut HashMap<String, HashMap<T, u32>>, stat_group: String, prop: T) {
			count.entry(stat_group).or_insert_with(HashMap::new)
				.entry(prop)
				.and_modify(|e| { *e += 1 }).or_insert(1);
		};
		fn build_hash<T: Hash + Eq + PartialEq + Clone>(cat: &str) -> HashMap<String, HashMap<T, u32>> {
			[	(format!("{}pc", cat), HashMap::new()), 
				(format!("{}pn", cat), HashMap::new()),
				(format!("{}nc", cat), HashMap::new()), 
				(format!("{}nn", cat), HashMap::new())].iter().cloned().collect()
		}
		let calc_helix = |x: &mut HashMap<ResFeature, u32>, sa| { 
			*x.entry(ResFeature {ss: 'H', sa}).or_insert(0) + 
			*x.entry(ResFeature {ss: 'G', sa}).or_insert(0) + 
			*x.entry(ResFeature {ss: 'I', sa}).or_insert(0) };

		let mut dist_full_total: HashMap<String, HashMap<ResFeature, u32>> = build_hash("l");
		let mut dist_per_protein_total: HashMap<String, HashMap<ResFeature, u32>> = build_hash("p");
		let mut n: u64 = 0;

		for prot in &self.proteins {
			let mut dist_protein_count = Vec::new();
			let mut pos_protein_index: Vec<Option<bool>> = vec![None; prot.base.residues.len()];
			let mut neg_protein_index: Vec<Option<bool>> = vec![None; prot.base.residues.len()];
			for (records, align, taed_id) in prot.get_families() {
				let mut dist_family_count: Vec<Option<HashMap<String, HashMap<ResFeature, u32>>>> = vec![ 
					Some(build_hash("l")); records.len()];
				let mut pdb_pos: usize = 0;
				for align_pos in 0..align.len() as usize {
					if pdb_pos >= prot.base.residues.len() {
						self.bad_counts.entry("Skipped Residues".to_string()).and_modify(|x| {*x += 1}).or_insert(1);
						continue;
					}
					if align.members[0].sites[align_pos] == '-' {
						self.bad_counts.entry("Skipped Residues".to_string()).and_modify(|x| {*x += 1}).or_insert(1);
						continue
					} else if self.aa_1_to_3.get(&prot.base.residues[pdb_pos]) == None {
						self.bad_counts.entry("Unidentified Residues".to_string()).and_modify(|x| {*x += 1}).or_insert(1);
						pdb_pos += 1;
						continue
					} else {
						self.bad_counts.entry("Compared Residues".to_string()).and_modify(|x| {*x += 1}).or_insert(1);
					}

					let exposure: char = 
						if prot.get_solvent_accessibility(pdb_pos).unwrap() / self.wilke[&self.aa_1_to_3[&prot.base.residues[pdb_pos]]] > 0.20 
							{ 'F' } else { 'P' };
					let mut k = 0;

					for record in records {
						if record.is_paml_generated() {
							dist_family_count[k] = None;
							k += 1;
							continue
						} else if record.ancestor as usize >= align.members.len() {
							if !self.bad_datasources.entry("Alignment Build Issues".to_string()).or_insert_with(HashSet::new)
								.contains(taed_id) {
								eprintln!("Mismatch Between Alignment Size [{}] and Ancestor Index [{}] for {}.", 
										align.members.len(), record.ancestor, taed_id);
								self.bad_datasources.entry("Alignment Build Issues".to_string()).or_insert_with(HashSet::new)
									.insert(taed_id.clone());
							}
							k += 1;
							continue
						};

						let ancestor = &align.members[record.ancestor as usize].sites[align_pos];
						let descendent = &align.members[record.id as usize].sites[align_pos];
						if *ancestor == '-' || *descendent == '-' { continue }
						let is_changed = *ancestor != *descendent;

						let mut selection = match record.dNdS() {
							0.000_000_1..=0.500_000_000 => { 
								neg_protein_index[pdb_pos] = match neg_protein_index[pdb_pos] {
									Some(x) => Some(x || is_changed),
									None => Some(is_changed)
								};
								"ln".to_string() 
							},
							1.000_000_0..=1000.0 => { 
								pos_protein_index[pdb_pos] = match pos_protein_index[pdb_pos] {
									Some(x) => Some(x || is_changed),
									None => Some(is_changed)
								};
								"lp".to_string()
							},
							_ => { continue } // eprint!("^{}", x);
						};

						selection.push(if is_changed { 'c' } else { 'n' });

						// Increment both this-family and overall totals.  Handle both
						// Secondary Structure and Solvent Exposure.
						if let Some(Some(ref mut branch)) = dist_family_count.get_mut(k) {
							increment(branch, selection.clone(), ResFeature { ss: prot.base.ss_simple[pdb_pos], sa: exposure });
						}
						increment(&mut dist_full_total, selection.clone(), ResFeature { ss: prot.base.ss_simple[pdb_pos], sa: exposure });
						k += 1;
					}					
					
					// Move the PDB positiion forward as we've dealt with it and haven't skipped it due to misalignment.
					pdb_pos += 1;
				}
				for branch in dist_family_count {
					if let Some(branch_counts) = branch { n += 1; dist_protein_count.push(branch_counts); }
				}
			}
			let mut dist_per_protein_count: HashMap<String, HashMap<ResFeature, u32>> = build_hash("p");

			for pdb_pos in 0..pos_protein_index.len() {
				if prot.base.residues[pdb_pos] == '-' {
					self.bad_counts.entry("Per-Protein Skipped Residues".to_string()).and_modify(|x| {*x += 1}).or_insert(1);
					continue
				} else if self.aa_1_to_3.get(&prot.base.residues[pdb_pos]) == None {
					self.bad_counts.entry("Per-Protein Unidentified Residues".to_string()).and_modify(|x| {*x += 1}).or_insert(1);
					continue
				} else {
					self.bad_counts.entry("Per-Protein Compared Residues".to_string()).and_modify(|x| {*x += 1}).or_insert(1);
				}

				let exposure: char =
					if prot.get_solvent_accessibility(pdb_pos).unwrap() / self.wilke[&self.aa_1_to_3[&prot.base.residues[pdb_pos]]] > 0.20
						{ 'F' } else { 'P' };
				let mut increment_list = Vec::new();
				if let Some(changed) = pos_protein_index[pdb_pos] 
					{ increment_list.push(if changed { "ppc".to_string() } else { "ppn".to_string() }) }
				if let Some(changed) = neg_protein_index[pdb_pos] 
					{ increment_list.push(if changed { "pnc".to_string() } else { "pnn".to_string() }) }
				
				for prefix in increment_list {
					increment(&mut dist_per_protein_count, prefix.clone(), ResFeature { ss: prot.base.ss_simple[pdb_pos], sa: exposure });
					increment(&mut dist_per_protein_total, prefix.clone(), ResFeature { ss: prot.base.ss_simple[pdb_pos], sa: exposure });
				}
			}

			dist_protein_count.push(dist_per_protein_count);
			self.dist_protein_counts.push(dist_protein_count);
		}

		// Inserts Mean value for per-lineage.
		for (stat_type, dist) in dist_full_total.iter_mut() {
			let total: u32 = dist.values().sum();
			let x = calc_helix(dist, 'F'); dist.insert(ResFeature {ss: 'A', sa: 'F'}, x);
			let x = calc_helix(dist, 'P'); dist.insert(ResFeature {ss: 'A', sa: 'P'}, x);
			let name = stat_type;
			self.stats.insert(format!("{}e", stat_type), ProteinAnalytics::distribution_calc(
				total as f64, 
				dist, name.to_string(), n));
			self.dist_net_counts.insert(stat_type.clone(), dist.clone());
		}

		// Inserts Mean value for per-protein.
		for (stat_type, dist) in dist_per_protein_total.iter_mut() {
			let total: u32 = dist.values().sum();
			let x = calc_helix(dist, 'F'); dist.insert(ResFeature {ss: 'A', sa: 'F'}, x);
			let x = calc_helix(dist, 'P'); dist.insert(ResFeature {ss: 'A', sa: 'P'}, x);
			let name = stat_type;
			self.stats.insert(format!("{}e", stat_type), ProteinAnalytics::distribution_calc(
				total as f64, 
				dist, name.to_string(), self.proteins.len() as u64));
			self.dist_net_counts.insert(stat_type.clone(), dist.clone());
		}
		
		// Inserts Rate Calcs
		let calc_list = vec![(&dist_full_total["lpc"], &dist_full_total["lpn"], "lpre", n),
							(&dist_full_total["lnc"], &dist_full_total["lnn"], "lnre", n),
							(&dist_per_protein_total["ppc"], &dist_per_protein_total["ppn"], "ppre", self.proteins.len() as u64),
							(&dist_per_protein_total["pnc"], &dist_per_protein_total["pnn"], "pnre", self.proteins.len() as u64)];
		for (sub, no_sub, rate, n) in calc_list {
			self.stats.insert(rate.to_string(), ProteinAnalytics::rate_calc(
				no_sub, sub, rate.to_string(), n));
		}

		// Inserts Totals
		for stat_type in ["lp", "ln", "pp", "pn"].iter() {
			let mut totals: HashMap<ResFeature, u32> = self.dist_net_counts[&format!("{}n", stat_type)].clone();
			for (key, value) in totals.iter_mut() {	
				*value += *self.dist_net_counts.entry(format!("{}p", stat_type)).or_insert_with(HashMap::new).entry(key.clone()).or_insert(0);		
			}
			let total: u32 = totals.values().sum::<u32>() - totals[&ResFeature {ss: 'A', sa: 'F'}] - totals[&ResFeature {ss: 'A', sa: 'P'}];
			let mut new_stat = (*stat_type).to_string();
			new_stat.push_str("ae");
			self.stats.insert(new_stat.clone(), ProteinAnalytics::distribution_calc(
				total as f64, 
				&mut totals, new_stat.clone(), 
				if new_stat.starts_with('l') { n } else { self.proteins.len() as u64 }));
		}
	}

	pub fn build_stats(&mut self) {
		// Build std deviations and per-protein stats
		let mut index: usize = 0;
		let mut counter: HashMap<String, HashMap<String, f64>> = HashMap::new();

		for prot in self.dist_protein_counts.iter_mut() {
			let mut j = 0;
			for prot_count in prot.iter_mut() {
				let mut prot_stat: HashMap<String, StatGroup> = HashMap::new();

				let mut name = self.proteins[index].id().to_string(); name.push('_'); name.push_str(&j.to_string());

				// Sum setup for standard deviation.
				let calc_dev = |mean: &StatGroup, item: &StatGroup, sums: &mut HashMap<String, f64>| {
					for name in mean.stats.keys() {
						let counter = sums.entry(name.to_string()).or_insert(0.00);
						let cmp_val = match item.stats.get(name) { Some(x) => x.value as f64, None => 0.00 };
						*counter += (mean.stats[name].value - cmp_val).powi(2);
					}
				};

				// Uses raw counts to caculate the % distributions.
				for (stat_type, dist) in prot_count.iter_mut() {
					// Calc total; if 0 ignore this set.
					let total_residues: u32 = dist.values().sum();
					if total_residues == 0 { continue };

					// Sums different helix types.
					let x = *dist.entry(ResFeature {ss: 'H', sa: 'F'}).or_insert(0) + 
							*dist.entry(ResFeature {ss: 'G', sa: 'F'}).or_insert(0) + 
							*dist.entry(ResFeature {ss: 'I', sa: 'F'}).or_insert(0); dist.insert(ResFeature {ss: 'A', sa: 'F'}, x);
					let x = *dist.entry(ResFeature {ss: 'H', sa: 'P'}).or_insert(0) + 
							*dist.entry(ResFeature {ss: 'G', sa: 'P'}).or_insert(0) + 
							*dist.entry(ResFeature {ss: 'I', sa: 'P'}).or_insert(0); dist.insert(ResFeature {ss: 'A', sa: 'P'}, x);

					// Calculates raw value.
					if stat_type.len() > 3 { eprintln!("!!{}", stat_type); }
					let extended_stat_type = format!("{}e", stat_type);
					prot_stat.insert(extended_stat_type.clone(), ProteinAnalytics::distribution_calc(
						total_residues as f64, 
						dist, name.clone(), 0));
					
					// Adds difference to std_dev calculation
					calc_dev(&self.stats[&extended_stat_type], &prot_stat[&extended_stat_type], counter.entry(extended_stat_type).or_insert_with(HashMap::new));
				}

				// Does substitution rate calculations for std_dev.
				let prefix = if prot_count.contains_key("lpc") { "l".to_string() } else { "p".to_string() };
				let set = vec![(format!("{}pc", prefix), format!("{}pn", prefix), format!("{}pre", prefix)), 
								 (format!("{}nc", prefix), format!("{}nn", prefix), format!("{}nre", prefix))];
	
				//println!("!!!!{:?}\n!!!{:?}\n!!{:?}\n!{:?}", self.stats.keys(), prot_stat.keys(), counter.keys(), prot_count.keys());
				for (sub, no_sub, rate) in set {
					// We've already got the Helix values so we just insert the local stats and then calculate the group stats.
					prot_stat.insert(rate.clone(), 
						ProteinAnalytics::rate_calc(prot_count.get(&no_sub).unwrap(), prot_count.get(&sub).unwrap(), name.clone(), 0));
					
					calc_dev(&self.stats[&rate], &prot_stat[&rate], counter.entry(rate.clone()).or_insert_with(HashMap::new));
				}
				j += 1

			}
			index += 1; 
		}

		// Assign Standard Deviations to existing stats.
		for (stat_type, counter_group) in counter.iter_mut() {
			let stat_group = self.stats.entry(stat_type.clone()).or_insert_with(|| StatGroup::new(stat_type.clone(), 0));
			for (_name, stat) in stat_group.stats.iter_mut() {
				//println!("{:?}:{}", counter_group.keys(), stat.field);
				stat.std_dev = (counter_group[&stat.field] / stat_group.n as f64).sqrt();
			}
		}
	}
}