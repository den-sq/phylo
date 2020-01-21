//! Simple statistical caluclations for mean / median / mode used elsewhere in the module.
extern crate num;
use std::ops::{Add, Div};
use std::default::Default;
use self::num::FromPrimitive;
use std::fmt;
use std::error;
use std::io;
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufReader, BufRead, Write};
use serde::{Serialize};

#[derive(Debug, Clone)]
/// Error for attempting to calculate a statistical value for an empty slice.
pub struct EmptySliceError;

// Generation of an error is completely separate from how it is displayed.
impl fmt::Display for EmptySliceError {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "No Values in Slice")
    }
}

// This is important for other errors to wrap this one.
impl error::Error for EmptySliceError {
    fn source(&self) -> Option<&(dyn error::Error + 'static)> {
        // Generic error, underlying cause isn't tracked.
        None
    }
}


/// Parses list of names with optional rename targets
/// List is from a file with a name and an (optional) rename target on separate lines, separated by tabs.
/// Inspired by Hunter Cameron's prune_phylo_tree.py, but BioPython can't handle some valid tree files.
/// 
/// # Parameters
/// 
/// filename: String containing filename of names.
/// 
/// # Return
/// 
/// Result object containing renaming hashmap or IO error.
pub fn parse_name_list(filename: &str) -> Result<HashMap<String, String>, io::Error> {
	let mut name_list: HashMap<String, String> = HashMap::new();
	match File::open(filename) {
		Ok(file) => {
			let buf_reader = BufReader::new(file);
			for line in buf_reader.lines() {
				let rename: Vec<String> = match line {
					Ok(names) => names.split('\t').map(|s| s.to_string()).collect(),
					Err(e) => return Err(e)
				};
				if rename.len() == 1 {
					name_list.insert(rename[0].clone(), rename[0].clone());
				} else if rename.len() == 2 {
					name_list.insert(rename[0].clone(), rename[1].clone());
				} else {
					// TODO add error
				}
			}
		},
		Err(e) => panic!("Could Not Open Name List, Exiting: {}", e)
	}

	Ok(name_list)
}

/// Gets the mean value given an iterator.
/// 
/// # Panic
/// 
/// Will panic if no values in iterator.
pub fn arithmatic_mean<I, V>(vals: I) -> V
	where I: Iterator, 
		V: Add<I::Item, Output = V> + Div<V, Output = V> + Default + FromPrimitive {
	let mut total: V = Default::default();
	let mut count = 0;
	for item in vals {
		total = total + item;
		count += 1;
	}
	total / FromPrimitive::from_usize(count).unwrap()
}

/// Gets the median value in a slice of values.
/// 
/// Slice can be unsorted but will be sorted in place.
/// 
/// # Errors
/// 
/// Slice must not be empty.
pub fn arithmatic_median<V>(vals: &mut [V]) -> Result<V, EmptySliceError>
	where V: Add<V, Output = V> + fmt::Debug + Div<V, Output = V> + FromPrimitive + Copy + std::cmp::PartialOrd {
	if vals.is_empty() {
		Err(EmptySliceError {})
	} else if vals.len() == 1 {
		Ok(vals[0])
	} else {
		//println!("{:?}", vals);
		vals.sort_by(|a, b| a.partial_cmp(b).unwrap());
		// Bitwise And - Is Even
		if vals.len() & 1 == 0 {
			let x = vals.len() / 2;
			Ok((vals[x] + vals[x - 1]) / FromPrimitive::from_usize(2).unwrap())
		} else {
			Ok(vals[(vals.len() - 1) / 2])
		}
	}
	
}

/// Writes serializable data to a datastore.
pub fn write_serializable_data<T:Serialize> (data: &T, file_path: &str) {
	let file_name = file_path.split('/').last().unwrap().to_string();
	match File::create(file_path) {
		Ok(mut datastore) =>
			match serde_json::to_string(data) { 
				Ok(data_dump) => if let Err(e) = datastore.write_all(data_dump.as_bytes()) { 
					eprintln!("Failed to Write Records ({}):{}", file_name, e);
				},
				Err(e) => eprintln!("Failed to Write Records ({}):{}", file_name, e)
			},
		Err(e) => eprintln!("Failed to Write Records ({}):{}", file_name, e)
	};
}

/// Gets the maximum value in a slice of values.
/// 
/// # Errors
/// 
/// Slice must not be empty.
pub fn arithmatic_max<V>(vals: &[V]) -> Result<V, EmptySliceError>
	where V: FromPrimitive + Copy + std::cmp::PartialOrd {
	if vals.is_empty() {
		return Err(EmptySliceError {});
	}
	let mut max: V = vals[0];
	for item in vals {
		if *item > max {
			max = *item
		}
	}
	Ok(max)
}

#[derive(Copy, Clone, Debug, Default)]
/// Statistical Range Values.
pub struct Ranges {
    pub mean: f64,
    pub median: f64,
    pub max: f64,
}

impl Ranges {
	/// Generates a new, dummy Range with NAN values.
	pub fn new() -> Ranges {
		Ranges {
			mean: std::f64::NAN,
			median: std::f64::NAN,
			max: std::f64::NAN
		}
	}
}