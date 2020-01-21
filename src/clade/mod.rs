#![allow(dead_code)]
///! Data structures for use with phylogenetic trees.
use std::fs::File;
use std::collections::HashMap;
use std::io::BufReader;
use std::io::prelude::*;
use std::io;
use std::fmt;
use std::error::Error as StdError;
use serde::{Serialize, Deserialize};
use indexmap::IndexSet;

/// Errors representing possible parsing and display issues.
#[derive(Debug)]
pub enum Error {
    MissingMetadata(String),
    MalformedMetadata(String, CladeMeta)
}

impl fmt::Display for Error {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match self {
            Error::MissingMetadata(metadata) => 
                f.write_str(&format!("Metadata {} Not Found but Access Attempted.", metadata)),
            Error::MalformedMetadata(metaname, metavalue) => 
                f.write_str(&format!("Metadata {} Pointing to Incorrect Data Block {}", metaname, metavalue))
        }
    }
}
impl StdError for Error {
    fn description(&self) -> &str {
        match self {
            Error::MissingMetadata(_) => "Metadata Not Found but Access Attempted.",
            Error::MalformedMetadata(_, _) => "Metadata Pointing to Incorrect Data Block"
        }
    }
}



#[derive(Debug,Clone,Serialize,Deserialize)]
/// Handles metadata for tree nodes, matching standard NHX/Newick Metadata Information
pub enum CladeMeta {
	Name(String),
	Color(String),
	BranchLength(f64),
	Bootstrap(i64),
	Species(String),
	Taxonomy(i64),
	EC(String),
	IsDuplication(bool),
	Orthologous(i64),
	SuperOrthologous(i64),
	LogLikelihood(f64),
	SignificantlyWorse(bool),
	Collapse(bool),
	DNDS(f64),
	Type(char),
	UserDef(String),
}

impl std::default::Default for CladeMeta {
	fn default() -> CladeMeta { CladeMeta::UserDef("".to_string()) }
}

impl fmt::Display for CladeMeta {
    /// Displays Contained Metadata in Simplke Format, Y/N for booleans.
	fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match *self {
			CladeMeta::IsDuplication(v) | CladeMeta::SignificantlyWorse(v) |
				CladeMeta::Collapse(v)										=> if v { f.write_str("Y") } else { f.write_str("N") },
			CladeMeta::Name(ref v) | CladeMeta::Species(ref v) | 
				CladeMeta::EC(ref v) | CladeMeta::UserDef(ref v) 			=> f.write_str(v),
			CladeMeta::BranchLength(v) | CladeMeta::LogLikelihood(v) 		=> f.write_str(&format!("{:6}", v).trim()),
			CladeMeta::DNDS(v)												=> f.write_str(&format!("{:4}", v).trim()),
			CladeMeta::Bootstrap(v) | CladeMeta::Taxonomy(v) |
				CladeMeta::Orthologous(v) | CladeMeta::SuperOrthologous(v)	=> f.write_str(&v.to_string()),
			CladeMeta::Type(v)												=> f.write_str(&v.to_string()),
			CladeMeta::Color(_)												=> Err(fmt::Error)
        }
    }
}


/// Recursively defined tree structure.
#[derive(Debug,Serialize,Deserialize,Clone)]
pub struct Clade
{
    pub children: Vec<Clade>,
    pub name: String,
    branch_length: f64,
    pub tags: HashMap<String, CladeMeta>
}

impl Clade {
	/// Creates a new Clade from all data.
    pub fn new(branch_length: f64, name: String, clades: Vec<Clade>, tags: HashMap<String, CladeMeta>) -> Clade {
        Clade {
            branch_length,
            name,
            children: clades,
            tags
        }
    }

	/// Creates a new Clade using empty / default values.
    pub fn bare() -> Clade {
        Clade {
            branch_length: 0.0,
            name: String::from(""),
            children: Vec::new(),
            tags: HashMap::new()
        }
    }

	/// Adds a single child node.
    pub fn add_child(&mut self, child: Clade) {
        self.children.push(child)
    }

	/// Adss multiple child nodes from a vector.
    pub fn add_children(&mut self, mut new_children: Vec<Clade>) {
        self.children.append(&mut new_children)
    }

    /// Gets the root note for this clade, which is itself.
	pub fn root(&mut self) -> &mut Clade {
		self
	}

    /// Determines if negative (purifying) selection occured on the branch leading to this node.
	pub fn is_neg_select(&self) -> Result<bool, Error> {
		match self.tags.get(&"DNDS".to_string()) {
			Some(CladeMeta::DNDS(value)) => { Ok((*value > 0.0) && (*value <= 0.5))	},
            Some(wrong_meta) => Err(Error::MalformedMetadata("DNDS".to_string(), (*wrong_meta).clone())), 
			None => Err(Error::MissingMetadata("DNDS".to_string()))
		}
	}

	/// Helper function that updates basic metadata for a node.
    fn update_node(is_branch_length: bool, current_string: String, node: &mut Clade) -> Result<&str, io::Error> {
        if is_branch_length {
            match current_string.parse::<f64>() {
                Ok(branch_length) => {
                    (*node).branch_length = branch_length;
                    Ok("")
                },                    
                Err(_e) => Err(io::Error::new(io::ErrorKind::InvalidInput, 
                                format!("Invalid NHX/Newick Format, Non-Numeric Branch Length [{}].", current_string)))
            }
        } else {
            (*node).name = current_string;
            Ok("")
        }
    }

	/// Helper function that reads NHX metadata given field and value.
	/// 
	/// # Parameters
	/// 
	/// nhx_field: Borrowed string slice containing the field identifier.
	/// nhx_value: Mutable string object moved into the metadata information.
	/// tags:  Metadata object
   	#[allow(non_snake_case)]
	fn handle_nhx_field(nhx_field: &str, mut nhx_value: String, tags: &mut HashMap<String, CladeMeta>) -> Result<(), io::Error> {
		match nhx_field {
				"B" => match nhx_value.parse::<i64>() {
						Ok(bootstrap) => (*tags).insert(String::from("B"), CladeMeta::Bootstrap(bootstrap)),
						Err(_e) => return Err(io::Error::new(io::ErrorKind::InvalidInput, 
							format!("Invalid NHX Format, Non-Integer bootstrap Value [{}].", nhx_value)))
					},
				"S" => (*tags).insert(String::from("S"), CladeMeta::Species(nhx_value)),
				"T" => match nhx_value.parse::<i64>() {
						Ok(taxonomy) => (*tags).insert(String::from("T"), CladeMeta::Taxonomy(taxonomy)),
						Err(_e) => return Err(io::Error::new(io::ErrorKind::InvalidInput, 
							format!("Invalid NHX Format, Non-Integer Taxonomy ID [{}].", nhx_value)))
					},
				"E" => (*tags).insert(String::from("E"), CladeMeta::EC(nhx_value)),
				"D" => (*tags).insert(String::from("D"), CladeMeta::IsDuplication(&nhx_value == "Y")),
				"O" => match nhx_value.parse::<i64>() {
						Ok(external_node) => (*tags).insert(String::from("O"), CladeMeta::Orthologous(external_node)),
						Err(_e) => return Err(io::Error::new(io::ErrorKind::InvalidInput, 
							format!("Invalid NHX Format, Non-Integer External Orthologous Node [{}].", nhx_value)))
					},
				"SO" => match nhx_value.parse::<i64>() {
						Ok(external_node) => (*tags).insert(String::from("SO"), CladeMeta::SuperOrthologous(external_node)),
						Err(_e) => return Err(io::Error::new(io::ErrorKind::InvalidInput, 
							format!("Invalid NHX Format, Non-Integer External Super-Orthologous Node Length [{}].", nhx_value)))
					},
				"L" => match nhx_value.parse::<f64>() {
						Ok(likelihood) => (*tags).insert(String::from("L"), CladeMeta::LogLikelihood(likelihood)),
						Err(_e) => return Err(io::Error::new(io::ErrorKind::InvalidInput, 
							format!("Invalid NHX Format, Non-Numeric Parental Log Likelihood [{}].", nhx_value)))
					},
				"Sw" => (*tags).insert(String::from("Sw"), CladeMeta::SignificantlyWorse(&nhx_value == "Y")),
				"Co" => (*tags).insert(String::from("Co"), CladeMeta::Collapse(&nhx_value == "Y")),
				"DNDS" => match nhx_value.parse::<f64>() {
						Ok(dN_dS) => (*tags).insert(String::from("DNDS"), CladeMeta::DNDS(dN_dS)),
						Err(_e) => return Err(io::Error::new(io::ErrorKind::InvalidInput, 
							format!("Invalid NHX Format, Non-Numeric DNDS [{}].", nhx_value)))
					},
				"TYPE" => {
						match nhx_value.pop() {
							Some(clade_type) => (*tags).insert(String::from("TYPE"), CladeMeta::Type(clade_type)),
							None => (*tags).insert(String::from("TYPE"), CladeMeta::Type(' '))
						}
						
				},
				_ => if true  {
						(*tags).insert(String::from(nhx_field), CladeMeta::UserDef(nhx_value))
					} else {
						return Err(io::Error::new(io::ErrorKind::InvalidInput, 
							format!("Invalid NHX Format, Unidentified Tag [{0}] in NHX Section.", nhx_field)))
					}
		};
		Ok(())
	}

	/// Parses a newick tree from a string.
	/// 
	/// # Parameters
	/// 
	/// contents: String containing a valid newick tree.
	/// 
	/// # Return Value
	/// 
	/// Result object holding either the resulting Clade object or an I/O error.
	pub fn parse_newick_string(contents: String) -> Result<Clade, io::Error> {
		let mut child_matrix: Vec<Vec<Clade>> = vec![Vec::new()];

        let mut clades: Vec<Clade> = vec![Clade::bare()];

        let mut current_string = String::new();
        let char_list = contents.chars();
        let mut in_branch_length = false;
        for c in char_list {
            match c {
                    ';' => {
						// Update metadata for parent node.
                        match clades.last_mut() {
                            Some(last) => { Clade::update_node(in_branch_length, current_string, last)?; },
                            None => return Err(io::Error::new(io::ErrorKind::InvalidInput, 
                                                "Invalid Newick Format, Unexpected ';' (Excess ')' likely."))                            
                        }
						// Last node in the working list is the root.
                        return Ok(clades.pop().expect("This Error Should Not Occur."))
                    },
                    '(' => {
						// Add a new node to the working list, and add another tier of children.
                        clades.push(Clade::bare());
                        child_matrix.push(Vec::new());
                    },
                    ')' => {
                        match clades.pop() {
                            Some(mut new_child) => 
                                match clades.last_mut() {
                                    Some(last) => {
										// The current working node is the last child.
                                        Clade::update_node(in_branch_length, current_string, &mut new_child)?;
                                        current_string = String::from("");
                                        in_branch_length = false;
                                        (*last).add_children(child_matrix.pop().unwrap());
                                        (*last).add_child(new_child);
                                    },
                                    None => return Err(io::Error::new(io::ErrorKind::InvalidInput, 
                                                        "Invalid Newick Format, Unexpected ')' [Not Matched to '(']."))
                                }
                            None => return Err(io::Error::new(io::ErrorKind::InvalidInput, 
                                                "Invalid Newick Format, Unexpected ')'  [Not Matched to '(']."))
                        }
                    },
                    ',' =>  match clades.pop() {
                                Some(mut new_child) => {
									// We're going to start a new child, so update the metadata.
                                    Clade::update_node(in_branch_length, current_string, &mut new_child)?;
                                    current_string = String::from("");
                                    in_branch_length = false;

									// Working node is another child node, move it to the current list of children.
                                    match child_matrix.last_mut() {
                                        Some(last) => (*last).push(new_child),
                                        None => return Err(io::Error::new(io::ErrorKind::InvalidInput, 
                                                    "Invalid Newick Format, Unexpected ','."))
                                    }

									// Cue up the next child.
                                    clades.push(Clade::bare());
                                },
                                None => return Err(io::Error::new(io::ErrorKind::InvalidInput, "Invalid Newick Format, Unexpected ','."))
                            },
                    ':' => match clades.last_mut() {
                                Some(last) => {
									// We're switching from name to branch length, so assign the name and move forward.
                                    (*last).name = current_string;
                                    current_string = String::from("");
                                    in_branch_length = true
                                },
                                None => return Err(io::Error::new(io::ErrorKind::InvalidInput, 
                                                    "Invalid Newick Format, Unexpected ')' [Not Matched to '(']."))
                            },
                    ' ' => {}, 		// Spaces are ignored in newick, for formatting reasons.
					'\t' => {}, 	// So are tabs.
                    '_' => current_string.push(' '),	// Underscores are handled as spaces.
                    _ => current_string.push(c) 		// Other characters are just part of the name.
            }
        }

        Err(io::Error::new(io::ErrorKind::InvalidInput, "Invalid Newick Format, Unexpected EOF"))
	}

	/// Parses a newick file from a buffered reader. 
	/// 
	/// # Parameters
	/// 
	/// reader: Buffered Reader over an object implenting Read.
    fn parse_newick<T: Read>(mut reader: BufReader<T>) -> Result<Clade, io::Error> {
        let mut contents = String::new();
        reader.read_to_string(&mut contents)?;
        
		Clade::parse_newick_string(contents)
	}
	
	/// Export leaf values
	pub fn export_leaves(&self) -> HashMap<String, HashMap<String, CladeMeta>> {
		let mut tip_values: HashMap<String, HashMap<String, CladeMeta>> = HashMap::new();

		if self.children.is_empty() {
			tip_values.insert(self.name.clone(), self.tags.clone());
		}
		for child in self.children.iter() {
			tip_values.extend(child.export_leaves());
		}
		
		tip_values
	}

	/// Export leaf values to CSV
	pub fn write_leaves_csv(&self, file_name: &str) -> Result<(), io::Error> {
		let mut leaf_values: HashMap<String, HashMap<String, CladeMeta>> = self.export_leaves();
		let mut keys: IndexSet<String> = IndexSet::new();
		let mut out_string: String = String::from("Name");

		for (_leaf_name, leaf_data) in leaf_values.iter() {
			for (tag, _tag_value) in leaf_data.iter() {
				keys.insert(tag.clone());
			}
		}

		for key in keys.iter() {
			out_string.push(',');
			out_string.push_str(key);
		}

		out_string.push('\n');

		for (leaf_name, leaf_data) in leaf_values.iter_mut() {
			out_string.push_str(leaf_name);
			for key in keys.iter() {
				out_string.push_str(&format!(",{}", leaf_data.entry(key.clone()).or_default()));
			}
			out_string.push('\n');
		}

		match File::create(file_name) {
            Ok(mut handle) => {
                handle.write_all(out_string.as_bytes())?;
				Ok(())
            },
            Err(e) => Err(e)
        }
	}

	/// Parses an NHX file from a buffered reader. 
    fn parse_nhx<T: Read>(mut reader: BufReader<T>) -> Result<Clade, io::Error> {
        let mut contents = String::new();
        reader.read_to_string(&mut contents)?;
        let mut child_matrix: Vec<Vec<Clade>> = vec![Vec::new()];

        let mut clades: Vec<Clade> = vec![Clade::bare()];
        //let mut current_children: Vec<Clade> = Vec::new();

        let mut current_string = String::new();
        let char_list = contents.chars();

		// State flags what kind of data parser is going through.
        let mut in_branch_length = false;		// We're in the branch length area.
		let mut check_nhx_section = false;		// We're reading through (what we think is) the [&&NHX section but haven't reached the end.
		let mut nhx_section = false;			// We're in an NHX Section.
		let mut newick_meta_ready = true;		// We're ready for a newick metadata section; not true if we've just exited an nhx section.
		
		let mut nhx_prompt = String::new();
		let mut nhx_field = String::new();


        for c in char_list {
			if nhx_section {
				match c {
					'=' => {
						nhx_field = current_string;
						current_string = String::new();
					},
					':' => {
						if &nhx_field != "" {
							if let Some(node) = clades.last_mut() {
								Clade::handle_nhx_field(&nhx_field, current_string, &mut node.tags)?;
								nhx_field = String::new();
								current_string = String::new();
							}
						}					

					},
					']' => {
						if &nhx_field != "" {
							if let Some(node) = clades.last_mut() {
								Clade::handle_nhx_field(&nhx_field, current_string, &mut node.tags)?;
								nhx_field = String::new();
								current_string = String::new();
								newick_meta_ready = false;
							}
						}
						nhx_section = false
					},
                    ' ' | '\t' | '\n' | '\r' | '\x0C' => {}, 	// Spatial characters (\x0C is line feed) are ignored
					_ => current_string.push(c) 
				}
			} else {
				match c {
					'[' => {
						check_nhx_section = true;
						nhx_prompt.push(c);
					},
					'&' | 'N' | 'H' => 
						if check_nhx_section {
							nhx_prompt.push(c);
							if "[&&NHX".find(&nhx_prompt) == None {
								check_nhx_section = false;
								current_string.push_str(&nhx_prompt);
								nhx_prompt = String::new();
							}
						} else {
							current_string.push(c);
						}
					'X' => 
						if check_nhx_section {
							nhx_prompt.push(c);
							if nhx_prompt == "[&&NHX" {		// Is this comparison good?
								check_nhx_section = false;
								nhx_section = true;
								nhx_prompt = String::new();
								match clades.last_mut() {
									Some(last) => { Clade::update_node(in_branch_length, current_string, last)?; },
									None => return Err(io::Error::new(io::ErrorKind::InvalidInput, 
														"Invalid NHX Format, Unexpected [&&NHX Section"))                            
								}
								current_string = String::from("");
								in_branch_length = false;
							} else {
								check_nhx_section = false;
								current_string.push_str(&nhx_prompt);
								nhx_prompt = String::new();								
							}
						} else {
							current_string.push(c);
						}
                    ';' => {
						// Update metadata for parent node.
                        match clades.last_mut() {
                            Some(last) => { Clade::update_node(in_branch_length, current_string, last)?; },
                            None => return Err(io::Error::new(io::ErrorKind::InvalidInput, 
                                                "Invalid NHX Format, Unexpected ';' (Excess ')' likely."))                            
                        }
						// Last node in the working list is the root.
                        return Ok(clades.pop().expect("This Error Should Not Occur."))
                    },
                    '(' => {
						// Add a new node to the working list, and add another tier of children.
                        clades.push(Clade::bare());
                        child_matrix.push(Vec::new());
                    },
                    ')' => {
                        match clades.pop() {
                            Some(mut new_child) => 
                                match clades.last_mut() {
                                    Some(last) => {
										// The current working node is the last child.
										if newick_meta_ready {
											Clade::update_node(in_branch_length, current_string, &mut new_child)?;
											current_string = String::from("");
											in_branch_length = false;
										}
                                        (*last).add_children(child_matrix.pop().unwrap());
                                        (*last).add_child(new_child);
										newick_meta_ready = true
                                    },
                                    None => return Err(io::Error::new(io::ErrorKind::InvalidInput, 
                                                        "Invalid NHX Format, Unexpected ')' [Not Matched to '(']."))
                                }
                            None => return Err(io::Error::new(io::ErrorKind::InvalidInput, 
                                                "Invalid NHX Format, Unexpected ')'  [Not Matched to '(']."))
                        }
                    },
                    ',' =>  match clades.pop() {
                                Some(mut new_child) => {
									// We're going to start a new child, so update the metadata if we're not just out of nhx.
									if newick_meta_ready {
										Clade::update_node(in_branch_length, current_string, &mut new_child)?;
										current_string = String::from("");
										in_branch_length = false;
									}

									// Working node is another child node, move it to the current list of children.
                                    match child_matrix.last_mut() {
                                        Some(last) => (*last).push(new_child),
                                        None => return Err(io::Error::new(io::ErrorKind::InvalidInput, 
                                                    "Invalid NHX Format, Unexpected ','."))
                                    }

									// Cue up the next child.
                                    clades.push(Clade::bare());
									newick_meta_ready = true
                                },
                                None => return Err(io::Error::new(io::ErrorKind::InvalidInput, "Invalid NHX Format, Unexpected ','."))
                            },
                    ':' => match clades.last_mut() {
                                Some(last) => {
									// We're switching from name to branch length, so assign the name and move forward.
                                    (*last).name = current_string;
                                    current_string = String::from("");
                                    in_branch_length = true
                                },
                                None => return Err(io::Error::new(io::ErrorKind::InvalidInput, 
                                                    "Invalid NHX Format, Unexpected ')' [Not Matched to '(']."))
                            },
                    ' ' | '\t' | '\n' | '\r' | '\x0C' => {}, 	// Spatial characters (\x0C is line feed) are ignored
                    '_' => current_string.push(' '),			// Underscores are handled as spaces.
                    _ => current_string.push(c) 				// Other characters are just part of the name / branch length.
            	}
			}
        }

        Err(io::Error::new(io::ErrorKind::InvalidInput, "Invalid NHX Format, Unexpected EOF"))
    }

	/// Parses a nexus file from a buffered reader.
    fn parse_nexus<T: Read>(mut reader: BufReader<T>) -> Result<Clade, io::Error> {

        let mut contents = String::new();
        reader.read_to_string(&mut contents)?;

        Err(io::Error::new(io::ErrorKind::Other, "Not Yet Implemented"))
    }

	/// Reads a file into a Clade object.
    pub fn read(file_name: &str, file_type: &str) -> Result<Clade, io::Error> {

        let parse_func = match file_type {
            "newick" => Clade::parse_newick,
            "nexus" => Clade::parse_nexus,
            "nhx" => Clade::parse_nhx,
			"extension" => match file_name.rsplit('.').next() {
				Some(extension) => match extension {
					"newick" => Clade::parse_newick,
					"nexus" => Clade::parse_nexus,
					"nhx" => Clade::parse_nhx,
					_ => return Err(io::Error::new(io::ErrorKind::InvalidInput, "Unknown Extension"))

				},
				None => return Err(io::Error::new(io::ErrorKind::InvalidInput, "Invalid Filename"))
			}
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
}

impl Clade {

	/// Adds new annotations to the named nodes in the tree.
	/// Modifies tree in place.
	/// 
	/// # Parameters
	/// 
	/// annotations:  Hashmap of Names to match and the annotations to merge into it.
	/// Annotations already in the tree will be overwritten.
	/// 
	pub fn annotate(&mut self, annotations: &HashMap<String, HashMap<String, CladeMeta>>) {
		for name in annotations.keys() {
			if self.name.find(&name.replace('_', " ")).is_some() {
				self.tags.extend(annotations[name].clone());
			};
		}

		for child in self.children.iter_mut() {
			child.annotate(&annotations)
		}
	}

	pub fn annotate_from_csv(&mut self, file_name: &str, separator: &str) -> Result<(), io::Error> {
		let mut annotations: HashMap<String, HashMap<String, CladeMeta>> = HashMap::new();

		match File::open(file_name) {
			Ok(file) => {
				let buf_reader = BufReader::new(file);
				for line in buf_reader.lines() {
					let annotation_list: Vec<String> = match line {
						Ok(names) => names.split(separator).map(|s| s.to_string()).collect(),
						Err(e) => return Err(e)
					};
					let mut annotation_hash: HashMap<String, CladeMeta> = HashMap::new();
					for i in 1..annotation_list.len() {
						let nhx_entry: Vec<String> = annotation_list[i].split('=').map(|s| s.to_string()).collect();
						Clade::handle_nhx_field(&nhx_entry[0], nhx_entry[1].clone(), &mut annotation_hash)?;
					}
					annotations.insert(annotation_list[0].clone(), annotation_hash);
				}
			},
			Err(e) => return Err(e)
		}

		self.annotate(&annotations);
		Ok(())
	}


	/// Finds the most recent common ancestor for a list of species on a clade, returns its node.
	/// 
	/// # Parameters
	/// 
	/// name_list:
	/// 	List of node names to find the common ancestor of.  Values consumed by the function
	/// 	if found.
	/// 
	/// # Return Value
	/// 
	/// 	Option value holding reference to common ancestor of all given names.
	/// 	If value is "None", no common ancestor was found, and the name_list contains
	/// 	all the names that were not located. 
	pub fn find_common_ancestor(&self, name_list: &mut Vec<String>) -> Option<&Clade> {
		let mut trimmable_list: Vec<String> = name_list.clone();
		let mut partial_flag = false;
		for child in self.children.iter() {
			if let Some(_result) = child.check_for_names(&mut trimmable_list) {
				if !partial_flag {
					// This child has all the names we want to find, so it contains the common ancestor.
					return child.find_common_ancestor(name_list);
				} else {
					// We've found everything and are the common ancestor.
					return Some(self);
				}
			} else if trimmable_list.len() != name_list.len() {
				// This child has some of the names we want to find.
				// If there is a common ancestor, we're it.
				// Set flag so we know we're partway done.
				partial_flag = true;
			}
		}
		
		// List of names to prune.
		let mut to_prune: Vec<usize> = Vec::new();

		// Find all matches - searches for partial matches (handles subspecies)
		for i in 0..name_list.len() {
			if let Some(j) = self.name.find(&name_list[i].replace('_', " ")) {
				to_prune.push(i);
				println!("{}:{}:{}", j, name_list[i], self.name);
			}
		}

		to_prune.reverse();
		for pruning in to_prune {
			name_list.remove(pruning);
		}
		
		if name_list.is_empty() {
			// If we've found all the names, so this is the common ancestor.
			Some(self)
		} else {
			// Otherwise there is no common ancestor.
			None
		}
	}

	fn check_for_names(&self, name_list: &mut Vec<String>) -> Option<&Clade> {
		for child in self.children.iter() {
			if let Some(result) = child.check_for_names(name_list) {
				// If an ancestor has been found, keep passing it down.
				// The function call will trim the name list where appropriate.
				return Some(result)
			}
		}

		// List of names to prune.
		let mut to_prune: Vec<usize> = Vec::new();

		// Find all matches - searches for partial matches (handles subspecies)
		for i in 0..name_list.len() {
			if let Some(j) = self.name.find(&name_list[i].replace('_', " ")) {
				to_prune.push(i);
				println!("{}:{}:{}", j, name_list[i], self.name);
			}
		}

		to_prune.reverse();
		for pruning in to_prune {
			name_list.remove(pruning);
		}
		
		if name_list.is_empty() {
			// If we've found all the names, so return this node.
			Some(self)
		} else {
			// Otherwise we haven't yet.
			None
		}
	}


	/// Trims unwanted nodes from a tree.
	/// Trims tree in place.  Working off of Hunter Cameron's prune_phylo_tree.py, as
	/// BioPython can't handle some valid tree files.  Some modifications for my own uses.
	/// 
	/// # Parameters
	/// 
	/// &self
	/// nodes: Nodes Hashmap (e.g. generated by parse_name_list)
	/// 
	pub fn trim(&mut self, nodes: &HashMap<String, String>) {
		// store a list of nodes that need to be pruned to avoid pruning while traversing
		// Because HashMaps don't implement equality we can't do this the clean way
		// with &Clade & remove_item
		let mut to_prune: Vec<usize> = Vec::new();

		for (child_index, child) in self.children.iter_mut().enumerate() {
			// Trim child
			child.trim(nodes);

			if child.children.is_empty() {
				// Do it this way so we can do parital matches on the name list.
				// TODO: make name list regexes
				let mut found = false;
				for name in nodes.keys() {
					if child.name.find(&name.replace('_', " ")).is_some() {
						found = true;
						// If the key and value of the hashmap are the same, don't rename.
						if nodes[name] != *name {
							child.name = nodes[name].clone();
						};
						break;					
					}
				}
				if !found { to_prune.push(child_index) }
			}
		}

		to_prune.reverse();

		for clade_pos in to_prune {
			self.children.remove(clade_pos);
		}
	}

/*

        # after the pruning, if there is only one child, there is no point in having the internal node
        # must make a copy because root.clades is modified with each call to collapse
        for child in root.clades[:]:
            if len(child.clades) == 1:
                root.collapse(child)


*/


	/// Serializes a clade object into a newick-formatted string.
    fn serialize_newick(&self) -> String {
        if !self.children.is_empty() {
            let mut child_data = String::from("");
            for child in &self.children {
                child_data.push_str(&child.serialize_newick());
                child_data.push(',');
            }
            child_data.pop();
            format!("({}){}{}", child_data, self.name, 
                if self.branch_length != 0.0 { format!(":{}", self.branch_length) } else { String::from("") })
        } else {
            format!("{}{}", self.name, if self.branch_length != 0.0 { format!(":{}", self.branch_length) } else { String::from("") })
        }
    }

    /// Serializes a clade object into an nhx-formatted string.
	fn serialize_nhx(&self) -> String {
		let mut node = format!("{}{}", self.name, if self.branch_length != 0.0 { format!(":{}", self.branch_length) } else { String::from("") });
		
		if !self.tags.is_empty() {
			let mut tag_data = String::from("");
			for (tag, value) in &self.tags {
				tag_data.push_str(&format!(":{}={}", tag, value))
			}
			node = format!("{}[&&NHX{}]", node, tag_data)
		}

        if !self.children.is_empty() {
            let mut child_data = String::from("");
            for child in &self.children {
                child_data.push_str(&child.serialize_nhx());
                child_data.push(',');
            }
            child_data.pop();
            node = format!("({}){}", child_data, node)
        }
		node
	}

	/// Writes a clade object to a file.
    pub fn write(&self, file_name: &str, file_type: &str) -> Result<(), io::Error> {

        match File::create(file_name) {
            Ok(mut handle) => {
                let write_string = match file_type {
                    "newick" => format!("{};", self.serialize_newick()),
					"nhx" => format!("{};", self.serialize_nhx()),
					"extension" => match file_name.rsplit('.').next() {
						Some(extension) => match extension {
							"newick" => format!("{};", self.serialize_newick()),
							"nhx" => format!("{};", self.serialize_nhx()),
							_ => return Err(io::Error::new(io::ErrorKind::InvalidInput, "Unknown Extension"))
						},
						None => return Err(io::Error::new(io::ErrorKind::InvalidInput, "No Extension to Autodetect"))
					}
                    _ => return Err(io::Error::new(io::ErrorKind::InvalidInput, "Invalid Format"))
                };
                
                handle.write_all(write_string.as_bytes())?;
            },
            Err(e) => return Err(e)
        }
        Ok(())
    }
}
