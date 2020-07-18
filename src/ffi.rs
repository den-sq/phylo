//! # Phylo FFI
//! 
//! C Foreign Function Interface for use of the phylo phylogenetic data analysis library.
use crate::alignment::Alignment;
use crate::seq::Seq;
use std::process::Command;
use std::process::Stdio;


use std::cell::RefCell;
use std::error::Error;
use std::io::Write;
use std::slice;
use log::{warn, error};
use std::os::raw::{c_int, c_char};
use std::ffi::{CStr};
use std::ptr;
use crate::clade::Clade;
use crate::paml::PAML;
use crate::phenotype::{ProcessedData, PhylogeneticSet, get_matching_entries, get_divergent_entries};
use crate::prot::{ProteinStructure, ProteinData};

use serde_json::{Value};

use std::collections::HashMap;
use std::fs::File;
use std::io::BufReader;

use rand::prelude::*;
use rand::distributions::Uniform;

use std::io::stdout;
use csv;
pub use crate::utils::Ranges;

type BootstrapResult = HashMap<String, HashMap<String, u64>>;

// Error Handling from Michael F. Bryan's work.
thread_local!{
    static LAST_ERROR: RefCell<Option<Box<dyn Error>>> = RefCell::new(None);
}

/// Update the most recent error, clearing whatever may have been there before.
pub fn update_last_error<E: Error + 'static>(err: E) {
    error!("Setting LAST_ERROR: {}", err);

    {
        // Print a pseudo-backtrace for this error, following back each error's
        // cause until we reach the root error.
        let mut cause = err.source();
        while let Some(parent_err) = cause {
            warn!("Caused by: {}", parent_err);
            cause = parent_err.source();
        }
    }

    LAST_ERROR.with(|prev| {
        *prev.borrow_mut() = Some(Box::new(err));
    });
}

/// Retrieve the most recent error, clearing it in the process.
pub fn take_last_error() -> Option<Box<dyn Error>> {
    LAST_ERROR.with(|prev| prev.borrow_mut().take())
}

/// Calculate the number of bytes in the last error's error message **not**
/// including any trailing `null` characters.
#[no_mangle]
pub extern "C" fn last_error_length() -> c_int {
    LAST_ERROR.with(|prev| match *prev.borrow() {
        Some(ref err) => err.to_string().len() as c_int + 1,
        None => 0,
    })
}

/// Write the most recent error message into a caller-provided buffer as a UTF-8
/// string, returning the number of bytes written.
///
/// # Note
///
/// This writes a **UTF-8** string into the buffer. Windows users may need to
/// convert it to a UTF-16 "unicode" afterwards.
///
/// If there are no recent errors then this returns `0` (because we wrote 0
/// bytes). `-1` is returned if there are any errors, for example when passed a
/// null pointer or a buffer of insufficient size.
/// 
///  # Safety
/// 
/// buffer is a raw pointer; checked for null before error message is copied
/// into it with copy_nonoverlapping.
#[no_mangle]
pub unsafe extern "C" fn last_error_message(buffer: *mut c_char, length: c_int) -> c_int {
    if buffer.is_null() {
        warn!("Null pointer passed into last_error_message() as the buffer");
        return -1;
    }

    let last_error = match take_last_error() {
        Some(err) => err,
        None => return 0,
    };

    let error_message = last_error.to_string();

    let buffer = slice::from_raw_parts_mut(buffer as *mut u8, length as usize);

    if error_message.len() >= buffer.len() {
        warn!("Buffer provided for writing the last error message is too small.");
        warn!(
            "Expected at least {} bytes but got {}",
            error_message.len() + 1,
            buffer.len()
        );
        return -1;
    }

    ptr::copy_nonoverlapping(
        error_message.as_ptr(),
        buffer.as_mut_ptr(),
        error_message.len(),
    );

    // Add a trailing null so people using the string as a `char *` don't
    // accidentally read into garbage.
    buffer[error_message.len()] = 0;

    error_message.len() as c_int
}



/// FFI Function that reads a tree file of caller provided type at caller provided location, 
/// and returns the tree found in the file at that location, via a pointer to the root.
/// 
/// # Arguments
/// 
/// * 'file_name' - C-Style String holding the full file location.
/// * 'file_type' - C-Style String denoting the file type ("newick", "nexus", ""nhx", "extension")
/// 
/// # Safety
/// 
/// Parameters are all C-Style strings; can
/// 1: Fail Null Assertion (Causes Panic)
/// 2: Is Unsafe if Parameter does not point to a null-terminated string.
#[no_mangle]
pub unsafe extern "C" fn read_tree(file_name: *const c_char, file_type: *const c_char) -> *mut Clade {

	let inner_filename = {
		assert!(!file_name.is_null());
		CStr::from_ptr(file_name)
	};

	let inner_filetype = {
		assert!(!file_type.is_null());
		CStr::from_ptr(file_type)
	};

    let safe_filename = match inner_filename.to_str() {
		Ok(file_name) => String::from(file_name),
		Err(e) =>  {
			update_last_error(e);
			return Box::into_raw(Box::new(Clade::bare()))
		}
    };

	let safe_filetype = match inner_filetype.to_str() {
		Ok(file_type) => String::from(file_type),
		Err(e) =>  {
			update_last_error(e);
			return Box::into_raw(Box::new(Clade::bare()))
		}
	};

	match Clade::read(&safe_filename, &safe_filetype) {
		Ok(tree) => { Box::into_raw(Box::new(tree)) },
		Err(e) => {
			update_last_error(e);
			Box::into_raw(Box::new(Clade::bare()))
		}
	}
}

/// bootstraps a dataset based on passed in totals.
/// 
/// # Safety
/// 
/// Parameters are all C-Style strings; can
/// 1: Fail Null Assertion (Causes Panic)
/// 2: Is Unsafe if Parameter does not point to a null-terminated string.
#[no_mangle]
pub unsafe extern "C" fn bootstrap(path: *mut c_char, counts: *mut c_char, related_counts: *mut c_char, data: *mut c_char, matrix: *mut c_char, lookup: *mut c_char) {
    let temp_path = {
        assert!(!path.is_null());
        CStr::from_ptr(path)
    };

    let temp_counts = {
        assert!(!counts.is_null());
        CStr::from_ptr(counts)
    };

	let temp_related_counts = {
        assert!(!related_counts.is_null());
        CStr::from_ptr(related_counts)
    };

    let temp_data = {
        assert!(!data.is_null());
        CStr::from_ptr(data)
    };

    let temp_matrix = {
        assert!(!matrix.is_null());
        CStr::from_ptr(matrix)
    };

    let temp_lookup = {
        assert!(!lookup.is_null());
        CStr::from_ptr(lookup)
    };

    bootstrap_local(temp_path.to_str().unwrap(), temp_counts.to_str().unwrap(), temp_related_counts.to_str().unwrap(),
                    temp_data.to_str().unwrap(), temp_matrix.to_str().unwrap(), temp_lookup.to_str().unwrap())
}

fn bootstrap_local(path: &str, counts: &str, related_counts: &str, data: &str, matrix: &str, lookup: &str) {
    let totals: Value;
	let related_totals: Value;
    let calc_data: HashMap<String, ProcessedData>;
    let lookup_hash: Value;
    let mut matrix_hash: HashMap<String, Vec<String>> = HashMap::new();

    match File::open(format!("{}{}", path, counts)) {
        Ok(file) => {
            let buf_reader = BufReader::new(file);
            totals = match serde_json::from_reader(buf_reader) {
                Ok(data) => data,
                Err(e) => panic!("Could Not Parse Stored Data [totals], Exiting: {}", e)
            };
        },
        Err(e) => panic!("Could Not Open Stored Data [totals], Exiting: {}", e)
    };

	match File::open(format!("{}{}", path, related_counts)) {
        Ok(file) => {
            let buf_reader = BufReader::new(file);
            related_totals = match serde_json::from_reader(buf_reader) {
                Ok(data) => data,
                Err(e) => panic!("Could Not Parse Stored Data [totals], Exiting: {}", e)
            };
        },
        Err(e) => panic!("Could Not Open Stored Data [totals], Exiting: {}", e)
    };

    match File::open(format!("{}{}", path, data)) {
        Ok(file) => {
            let buf_reader = BufReader::new(file);
            calc_data = match serde_json::from_reader(buf_reader) {
                Ok(data) => data,
                Err(e) => panic!("Could Not Parse Stored Data [calc_data], Exiting: {}", e)
            };
        },
        Err(e) => panic!("Could Not Open Stored Data [calc_data], Exiting: {}", e)
    }

    match File::open(matrix) {
        Ok(file) => {
            let buf_reader = BufReader::new(file);
            let mut csv_reader = csv::Reader::from_reader(buf_reader);
            for result in csv_reader.records() {
                let record: Vec<String> = result.unwrap().iter().map(|s| s.to_owned()).collect();

                matrix_hash.insert(record[0].clone(), record);
            }
        },
        Err(e) => panic!("Could Not Open Stored Data [matrix], Exiting: {}", e)
    };

    match File::open(format!("{}{}", path, lookup)) {
        Ok(file) => {
            let buf_reader = BufReader::new(file);
            lookup_hash = match serde_json::from_reader(buf_reader) {
                Ok(data) => data,
                Err(e) => panic!("Could Not Parse Stored Data [lookup], Exiting: {}", e)
            };
        },
        Err(e) => panic!("Could Not Open Stored Data [lookup], Exiting: {}", e)
    };


    let tabulation = bootstrap_helper(&totals, &related_totals, &calc_data, &matrix_hash, &lookup_hash);

    match File::create("bootstrap.json") {
        Ok(mut datastore) => 
            match serde_json::to_string(&tabulation) { 
                Ok(data_dump) => match datastore.write_all(data_dump.as_bytes()) 
                    { Ok(()) => (), Err(e) => eprintln!("Could Write bootstrap Data {}", e) }, 
                Err(e) => eprintln!("Could Write bootstrap Data {}", e) 
            }, 
        Err(e) => eprintln!("Could Not Create File For bootstrap  Data: {}", e)
    };
}

///	FFI Function to Run Pathway and Gene Enrichment Analysis on phenomic and genomic data.
/// 
/// # Arguments
/// 
/// * 'path' - C-Style String holding the path to the directory for 'counts' and 'data'
/// * 'counts' - C-Style String Filename to save the total count across all species of selected genes and pathways, and genes and pathways.
/// * 'data' - C-Style String Filename to load the stored data to analyze.  Format is a .json HashMap/Dictionary (String Key, ProcessedData Value)
/// * 'matrix' - C-Style String full path to the .csv file with phenotype information.
/// 
/// # Note
/// 
/// Parameters are checked for null via assert!
/// If it encounters an error, runs update_last_error(), prints to stderr, and exits.
/// Result data is written to "analysis.json" in the "path" directory.
/// 
/// # Safety
/// 
/// Parameters are all C-Style strings; can
/// 1: Fail Null Assertion (Causes Panic)
/// 2: Is Unsafe if Parameter does not point to a null-terminated string.
#[no_mangle]
pub unsafe extern "C" fn basic_analysis(path: *mut c_char, counts: *mut c_char, data: *mut c_char, matrix: *mut c_char) {
     let temp_path = {
        assert!(!path.is_null());
        CStr::from_ptr(path)
    };

    let temp_counts = {
        assert!(!counts.is_null());
        CStr::from_ptr(counts)
    };

    let temp_data = {
        assert!(!data.is_null());
        CStr::from_ptr(data)
    };

    let temp_matrix = {
        assert!(!matrix.is_null());
        CStr::from_ptr(matrix)
    };

    let calc_data: HashMap<String, ProcessedData>;
    let mut matrix_hash: HashMap<String, Vec<String>> = HashMap::new();

    match File::open(format!("{}{}", temp_path.to_str().unwrap(), temp_data.to_str().unwrap())) {
        Ok(file) => {
            let buf_reader = BufReader::new(file);
            calc_data = match serde_json::from_reader(buf_reader) {
                Ok(data) => data,
                Err(e) => {
					eprintln!("Could Not Parse Stored Data [data], Exiting: {}", e);
					update_last_error(e);
					return
				}
            };
        },
        Err(e) => {
			eprintln!("Could Not Open Stored Data [data], Exiting: {}", e);
			update_last_error(e);
			return
		}
    }

    match File::open(&String::from(temp_matrix.to_str().unwrap())) {
        Ok(file) => {
            let buf_reader = BufReader::new(file);
            let mut csv_reader = csv::Reader::from_reader(buf_reader);
            for result in csv_reader.records() {
                let record: Vec<String> = result.unwrap().iter().map(|s| s.to_owned()).collect();

                matrix_hash.insert(record[0].clone(), record);
            }
        },
        Err(e) => {
			eprintln!("Could Not Open Stored Data [matrix], Exiting: {}", e);
			update_last_error(e);
			return
		}
    };

    let mut totals: HashMap<String, u64>  = 
    [   
        (String::from("gene_pos"), 0),
        (String::from("gene_count"), 0),
        (String::from("pathway_pos"), 0),
        (String::from("pathway_count"), 0),
    ].iter().cloned().collect();
    let analysis = enrichment_analysis(&mut totals, &calc_data, &matrix_hash);

    match File::create(format!("{}analysis.json", temp_path.to_str().unwrap())) {
         Ok(mut datastore) => 
            match serde_json::to_string(&analysis) { 
                Ok(data_dump) => match datastore.write_all(data_dump.as_bytes()) { 
					Ok(()) => (), 
					Err(e) => {
						eprintln!("Could Not Write Serialized Analytics Data {}", e);
						update_last_error(e);
						return
					}
				},	
                Err(e) => {
					eprintln!("Could Not Serialize Analytics Data {}", e);
					update_last_error(e);
					return
				}
            }, 
        Err(e) => eprintln!("Could Not Create File For Analysis Data: {}", e)       
    }

    match File::create(format!("{}{}", temp_path.to_str().unwrap(), temp_counts.to_str().unwrap())) {
        Ok(mut datastore) => 
            match serde_json::to_string(&totals) { 
                Ok(data_dump) => match datastore.write_all(data_dump.as_bytes()) { 
					Ok(()) => (), 
					Err(e) => {
						eprintln!("Could Not Write Serialized Analytics Data {}", e); 
						update_last_error(e);
					}, 
 	        	},
				Err(e) => {
					eprintln!("Could Not Serialize Analytics Data {}", e);
					update_last_error(e);
				}
			},

        Err(e) => {
			eprintln!("Could Not Create File For Totals: {}", e);
			update_last_error(e);
		}
    }
}

/// Runs genetic comparison analyses 
/// 
/// # Safety
/// 
/// Parameters are all C-Style strings; can
/// 1: Fail Null Assertion (Causes Panic)
/// 2: Is Unsafe if Parameter does not point to a null-terminated string.
#[no_mangle]
pub unsafe extern "C" fn related_basic_analysis(path: *mut c_char, counts: *mut c_char, data: *mut c_char, matrix: *mut c_char, lookup: *mut c_char) {
     let temp_path = {
        assert!(!path.is_null());
        CStr::from_ptr(path)
    };

    let temp_counts = {
        assert!(!counts.is_null());
        CStr::from_ptr(counts)
    };

    let temp_data = {
        assert!(!data.is_null());
        CStr::from_ptr(data)
    };

    let temp_matrix = {
        assert!(!matrix.is_null());
        CStr::from_ptr(matrix)
    };

    let temp_lookup = {
        assert!(!lookup.is_null());
        CStr::from_ptr(lookup)
    };

    let calc_data: HashMap<String, ProcessedData>;
    let lookup_hash: Value;
    let mut matrix_hash: HashMap<String, Vec<String>> = HashMap::new();

    match File::open(format!("{}{}", temp_path.to_str().unwrap(), temp_data.to_str().unwrap())) {
        Ok(file) => {
            let buf_reader = BufReader::new(file);
            calc_data = match serde_json::from_reader(buf_reader) {
                Ok(data) => data,
                Err(e) => panic!("Could Not Parse Stored Data, Exiting: {}", e)
            };
        },
        Err(e) => panic!("Could Not Open Stored Data, Exiting: {}", e)
    }

    match File::open(format!("{}{}", temp_path.to_str().unwrap(), temp_lookup.to_str().unwrap())) {
        Ok(file) => {
            let buf_reader = BufReader::new(file);
            lookup_hash = match serde_json::from_reader(buf_reader) {
                Ok(data) => data,
                Err(e) => panic!("Could Not Parse Stored Data, Exiting: {}", e)
            };
        },
        Err(e) => panic!("Could Not Open Stored Data, Exiting: {}", e)
    }

    match File::open(&String::from(temp_matrix.to_str().unwrap())) {
        Ok(file) => {
            let buf_reader = BufReader::new(file);
            let mut csv_reader = csv::Reader::from_reader(buf_reader);
            for result in csv_reader.records() {
                let record: Vec<String> = result.unwrap().iter().map(|s| s.to_owned()).collect();

                matrix_hash.insert(record[0].clone(), record);
            }
        },
        Err(e) => panic!("Could Not Open Stored Data, Exiting: {}", e)
    };

    let mut totals: HashMap<String, u64> =
    [   
        (String::from("related_gene_pos"), 0),
        (String::from("related_gene_count"), 0),
        (String::from("related_pathway_pos"), 0),
        (String::from("related_pathway_count"), 0)
    ].iter().cloned().collect();
    let analysis = related_enrichment_analysis(&mut totals, &calc_data, &matrix_hash, &lookup_hash);

    match File::create(format!("{}related_analysis.json", temp_path.to_str().unwrap())) {
         Ok(mut datastore) =>
            match serde_json::to_string(&analysis) {
                Ok(data_dump) => match datastore.write_all(data_dump.as_bytes())
                    { Ok(()) => (), Err(e) => eprintln!("Could Write Serialized Analytics Data {}", e) },
                Err(e) => eprintln!("Could Not Serialize Analytics Data {}", e)
            },
        Err(e) => eprintln!("Could Not Create File For Analysis Data: {}", e)
    }

    match File::create(format!("{}{}", temp_path.to_str().unwrap(), temp_counts.to_str().unwrap())) {
        Ok(mut datastore) =>
            match serde_json::to_string(&totals) {
                Ok(data_dump) => match datastore.write_all(data_dump.as_bytes())
                    { Ok(()) => (), Err(e) => eprintln!("Could Write Serialized Analytics Data {}", e) }
                Err(e) => eprintln!("Could Not Serialize Analytics Data {}", e)
            },
        Err(e) => eprintln!("Could Not Create File For Totals: {}", e)
    }
}

#[allow(non_snake_case)]
fn related_enrichment_analysis(totals: &mut HashMap<String, u64>, calc: &HashMap<String, ProcessedData>, matrix: &HashMap<String, Vec<String>>, lookup: &Value) -> HashMap<String, HashMap<String, HashMap<String, u64>>> {

    let mut results: HashMap<String, HashMap<String, HashMap<String, u64>>> =
    [
        (String::from("related_pathway_positive"), HashMap::new()),
        (String::from("related_gene_positive"), HashMap::new()),
    ].iter().cloned().collect();

    let mut counts: HashMap<String, HashMap<String, HashMap<String, u64>>> =
    [
        (String::from("related_pathway_count"), HashMap::new()),
        (String::from("related_gene_count"), HashMap::new()),
    ].iter().cloned().collect();

    let dummy: Vec<Value> = Vec::new(); // Dummy for references later.

    for (species, data) in calc.iter() {
        if matrix.contains_key(species) {
            for (_gene, _pathways) in data.pathways.iter() {
                for phenotype in matrix[species].iter() {
                    let edit_phenotype = if phenotype == "" { "empty" } else { phenotype };

                    let pos = results.get_mut("related_pathway_positive").unwrap().entry(String::from(edit_phenotype)).or_insert_with(HashMap::new);
                    let count = counts.get_mut("related_pathway_count").unwrap().entry(String::from(edit_phenotype)).or_insert_with(HashMap::new);

                    for (seq_id, dNdS) in &data.related_selection {
                        let pathway_list = match lookup["pathways_from_seq"][seq_id].as_array() {
                            Some(v) => v,
                            None => &dummy
                        };
						// Pathways
                        for pathway in pathway_list.iter().map(|x| String::from(x.as_str().unwrap())) {
                            if (*dNdS) >= 1.000 {
                                pos.entry(pathway.clone()).and_modify(|x| { *x += 1} ).or_insert(1);
                                totals.entry(String::from("related_pathway_pos")).and_modify(|x| { *x += 1});
                            }
                            count.entry(pathway).and_modify(|x| { *x += 1} ).or_insert(1);
                            totals.entry(String::from("related_pathway_count")).and_modify(|x| { *x += 1});    
						}
					}
                }
            }

			for phenotype in matrix[species].iter() {
				let edit_phenotype = if phenotype == "" { "empty" } else { phenotype };
				let pos = results.get_mut("related_gene_positive").unwrap()
								.entry(String::from(edit_phenotype)).or_insert_with(HashMap::new);
				let count = counts.get_mut("related_gene_count").unwrap()
								.entry(String::from(edit_phenotype)).or_insert_with(HashMap::new);

				for (seq_id, dNdS) in &data.related_selection {
					
					if lookup["gene_from_seq"][seq_id].is_string() {
						let gene = String::from(lookup["gene_from_seq"][seq_id].as_str().unwrap());
						if (*dNdS) >= 1.00 { pos.entry(gene.clone()).and_modify(|x| { *x += 1 }).or_insert(1); }
						count.entry(gene).and_modify(|x| { *x += 1 }).or_insert(1);
					} else {
						println!("Not String: gene_from_seq-{}: {}", seq_id, lookup["gene_from_seq"][seq_id]);
					}
				}
			}
        }
    }    

    results.extend(counts);
    results
}


///	Runs Pathway and Gene Enrichment Analysis on phenomic and genomic data.
/// 
/// # Arguments
/// 
/// * 'totals' - C-Style String holding the path to the directory for 'counts' and 'data'
/// * 'calc' - C-Style String Filename to save the total count across all species of selected genes and pathways, and genes and pathways.
/// * 'matrix' - C-Style String full path to the .csv file with phenotype information.
/// 
/// # Note
/// 
/// 
/// If it encounters an error, runs update_last_error(), prints to stderr, and exits.
/// Result data is written to "analysis.json" in the "path" directory.
/// 
#[allow(non_snake_case)]
fn enrichment_analysis(totals: &mut HashMap<String, u64>, calc: &HashMap<String, ProcessedData>, matrix: &HashMap<String, Vec<String>>) -> HashMap<String, HashMap<String, HashMap<String, u64>>> {

    let mut results: HashMap<String, HashMap<String, HashMap<String, u64>>> =
    [
        (String::from("pathway_positive"), HashMap::new()),
        (String::from("gene_positive"), HashMap::new()),
    ].iter().cloned().collect();

    let mut counts: HashMap<String, HashMap<String, HashMap<String, u64>>> =
    [
        (String::from("pathway_count"), HashMap::new()),
        (String::from("gene_count"), HashMap::new()),
    ].iter().cloned().collect();


    for (species, data) in calc.iter() {
        if matrix.contains_key(species) {
			println!("Pathway Count: {}; Selection Count: {}", data.pathways.len(), data.selection.len());

			// Pathway enrinchment is determined by the individual gene enrichments that constitute the pathway.
            for (gene, pathways) in data.pathways.iter() {

                if data.selection.contains_key(gene) {
                    for phenotype in matrix[species].iter() {
                        let edit_phenotype = if phenotype == "" { "empty" } else { phenotype };
                        let pos = results.get_mut("pathway_positive").unwrap().entry(String::from(edit_phenotype)).or_insert_with(HashMap::new);
                        let count = counts.get_mut("pathway_count").unwrap().entry(String::from(edit_phenotype)).or_insert_with(HashMap::new);
                        for dNdS in data.selection[gene].iter() {

							//println!("Inner Pathway Count: {}", pathways.len());
                            for pathway in pathways.iter() {
                                if (*dNdS) >= 1.000 {
									// When a gene sequence for a species is postiviely selected, it is marked
									// as enriching the pathway it is on for each phenotype it species has.
                                    pos.entry(pathway.clone()).and_modify(|x| { *x += 1} ).or_insert(1);
                                    totals.entry(String::from("pathway_pos")).and_modify(|x| { *x += 1});
                                }
                                count.entry(pathway.clone()).and_modify(|x| { *x += 1} ).or_insert(1);
                                totals.entry(String::from("pathway_count")).and_modify(|x| { *x += 1});
                            }
                        }
                    }
                }                
            }

			for (gene, dNdS_set) in data.selection.iter() {
				for phenotype in matrix[species].iter() {
					let edit_phenotype = if phenotype == "" { "empty" } else { phenotype };
					let pos = results.get_mut("gene_positive").unwrap()
									.entry(String::from(edit_phenotype)).or_insert_with(HashMap::new)
									.entry(gene.clone()).or_insert(0);
					let count = counts.get_mut("gene_count").unwrap()
									.entry(String::from(edit_phenotype)).or_insert_with(HashMap::new)
									.entry(gene.clone()).or_insert(0);
					for dNdS in dNdS_set {
						// When a gene sequence for a species is postiviely selected, it is marked 
						// as enriched for each phenotype associated with that species.
						if (*dNdS) >= 1.00 { (*pos) += 1; totals.entry(String::from("gene_pos")).and_modify(|x| { *x += 1}); }
						(*count) += 1;
						totals.entry(String::from("gene_count")).and_modify(|x| { *x += 1 });
					}
				}
			}
			stdout().flush().unwrap();
        }
    }

    results.extend(counts);
    results
}

#[allow(non_snake_case)]
fn bootstrap_helper(totals: &Value, related_totals: &Value, calc: &HashMap<String, ProcessedData>, matrix: &HashMap<String, Vec<String>>, lookup: &Value) -> HashMap<String, Vec<BootstrapResult>> {

    let mut bootstraps: HashMap<String, Vec<BootstrapResult>> = 
    [   (String::from("pathway"), Vec::new()),
        (String::from("related_pathway"), Vec::new()),
        (String::from("gene"), Vec::new()),
        (String::from("related_gene"), Vec::new())
    ].iter().cloned().collect();
    let dummy: Vec<Value> = Vec::new(); // Dummy for references later.
    let bootstrap_runs = 1000;
	let mut sampler = rand::thread_rng();
    for _i in 0..bootstrap_runs {
        for (_key, value) in bootstraps.iter_mut() {
            value.push(HashMap::new());
        }
    }
	for (_key, value) in bootstraps.iter_mut() {
		println!("{}", value.len())
	}
    // 
    let p_sample = Uniform::from(0..totals["pathway_count"].as_u64().unwrap());
    let p_threshold = totals["pathway_pos"].as_u64().unwrap();
    
    let pr_sample = Uniform::from(0..related_totals["related_pathway_count"].as_u64().unwrap());
    let pr_threshold = related_totals["related_pathway_pos"].as_u64().unwrap();

    let g_sample = Uniform::from(0..totals["gene_count"].as_u64().unwrap());
    let g_threshold = totals["gene_pos"].as_u64().unwrap();
    
    let gr_sample = Uniform::from(0..related_totals["related_gene_count"].as_u64().unwrap());
    let gr_threshold = related_totals["related_gene_pos"].as_u64().unwrap();

	for (species, data) in calc.iter() {
		if matrix.contains_key(species) {
			for (gene, pathways) in data.pathways.iter() {
				if data.selection.contains_key(gene) {
                    for phenotype in matrix[species].iter() {
                        let edit_phenotype = if phenotype == "" { "empty" } else { phenotype };
                        let counter = bootstraps.get_mut("pathway").unwrap();
                        for _dNdS in data.selection[gene].iter() {
                            for pathway in pathways.iter() {
                                for i in 0..bootstrap_runs {
                                    let t = counter.get_mut(i).unwrap()
                                        .entry(String::from(edit_phenotype)).or_insert_with(HashMap::new)
                                        .entry(pathway.clone()).or_insert(0);
                                    let res_num = p_sample.sample(&mut sampler);
                                    if res_num <= p_threshold {
                                        (*t) += 1;
                                    }
                                }
                            }
                        }
                    }
                }
                for phenotype in matrix[species].iter() {
                    let edit_phenotype = if phenotype == "" { "empty" } else { phenotype };
                    let counter = bootstraps.get_mut("related_pathway").unwrap();
                    for (seq_id, _dNdS) in data.related_selection.iter() {
                        // Pathways
                        let pathway_list = match lookup["pathways_from_seq"][seq_id].as_array() {
                            Some(v) => v,
                            None => &dummy
                        };
                        for pathway in pathway_list {
                            for i in 0..bootstrap_runs {
                                let t = counter.get_mut(i).unwrap()
                                    .entry(String::from(edit_phenotype)).or_insert_with(HashMap::new)
                                    .entry(String::from(pathway.as_str().unwrap())).or_insert(0);
                                let res_num = pr_sample.sample(&mut sampler);
                                if res_num <= pr_threshold {
                                    (*t) += 1;
                                }
                            }
                        }
                    }
                }
			}

			for (gene, dNdS_set) in data.selection.iter() {
				for phenotype in matrix[species].iter() {
					let edit_phenotype = if phenotype == "" { "empty" } else { phenotype };
					let counter = bootstraps.get_mut("gene").unwrap();
					for _val in dNdS_set {
						for i in 0..bootstrap_runs {
							let t = counter.get_mut(i).unwrap()
										.entry(String::from(edit_phenotype)).or_insert_with(HashMap::new)
										.entry(gene.clone()).or_insert(0);
							let res_num = g_sample.sample(&mut sampler);
							if res_num <= g_threshold {
								(*t) += 1;
							}
							if (*t) > 100 {
								println!("{}:{}:{}:{}", String::from(edit_phenotype), gene.clone(), i, (*t));                   
							}
						}
					}
				}
			}
			
			for phenotype in matrix[species].iter() {
				let edit_phenotype = if phenotype == "" { "empty" } else { phenotype };
				let g_counter = bootstraps.get_mut("related_gene").unwrap();
				for (seq_id, _dNdS) in data.related_selection.iter() {
                    for i in 0..bootstrap_runs {
                        // Genes
                        let t = g_counter.get_mut(i).unwrap()
                                    .entry(String::from(edit_phenotype)).or_insert_with(HashMap::new)
                                    .entry(seq_id.clone()).or_insert(0);
                        let res_num = gr_sample.sample(&mut sampler);
                        if res_num <= gr_threshold {
                            (*t) += 1;
                        }
                    }
				}
			}
		} else {
			println!("{} Not Found in Matrix", species);
		}
	}
    bootstraps
}

fn run_mafft (short: &str, long: &str) -> Result<Alignment<Seq>, std::io::Error> {
	let mafft_run = Command::new("mafft-profile")
					.args(&[short, long])
					.stdout(Stdio::piped())
					.output()
					.expect("Failed to Execute Process");
	Alignment::<Seq>::read_bytes(mafft_run.stdout, "fasta")
}

fn ptr_to_string(raw: *mut c_char) -> Result<String, std::str::Utf8Error> {
	match unsafe {
        assert!(!raw.is_null());	
		CStr::from_ptr(raw)
	}.to_str() {
		Ok(rust_str) => Ok(rust_str.to_string()),
		Err(e) => Err(e)
	}
}

#[no_mangle]
pub extern "C" fn get_protein_identity(alignment_path: *mut c_char, protein_path: *mut c_char) -> *const Ranges {
	match File::open(&ptr_to_string(protein_path).unwrap()) {
		Ok(file) => {
			let buf_reader = BufReader::new(file);
			
			let prot: ProteinStructure = match serde_json::from_reader(buf_reader) {
				Ok(protein) => protein,
				Err(e) => {
					eprintln!("{}", e);
					return Box::into_raw(Box::new(Ranges::new()));
				}
			};
			let second_path = &ptr_to_string(alignment_path).unwrap();
			prot.get_sequence().write("pdb_a.fasta", "fasta").unwrap();
			println!("{}", second_path);
			let new_align = run_mafft("pdb_a.fasta", second_path).unwrap();
			new_align.write("new_align.fasta", "fasta").unwrap();
			println!("Updated Alignment Length: {}", new_align.members.len());
			if new_align.members.len() == 1 { 
				eprintln!("One Length Alignment Found at {}", second_path);
				Box::into_raw(Box::new(Ranges::new()))
			} else {
				Box::into_raw(Box::new(new_align.get_root_identity()))
			}
		},
		Err(e) => {
			eprintln!("{:?}", e);
			Box::into_raw(Box::new(Ranges::new()))
		}
	}
}

#[no_mangle]
pub extern "C" fn get_identity(alignment_path: *mut c_char, alignment_format: *mut c_char) -> *const Ranges {
	let path = &ptr_to_string(alignment_path).unwrap();
	let format = &ptr_to_string(alignment_format).unwrap();
	match Alignment::<Seq>::read(path, format) {
		Ok(align) => Box::into_raw(Box::new(align.get_identity())),
		Err(e) => {
			eprintln!("Error:{}", e);
			Box::into_raw(Box::new(Ranges::new()))
		}
	}
}

/// Drops a ranges object.
/// 
/// # Safety
/// 
/// Parameters is raw pointer 
/// 1: Fail Null Assertion (Causes Panic)
/// 2: Is unpredicatable if pointer is not to a ranges object.
pub unsafe extern "C" fn drop_ranges(raw: *mut Ranges) {
    assert!(!raw.is_null());
    Box::from_raw(raw);
}

/// Processes Data for Convergent Evo Analysis
/// 
/// # Safety
/// 
/// Parameters are all C-Style strings; can
/// 1: Fail Null Assertion (Causes Panic)
/// 2: Is Unsafe if Parameter does not point to a null-terminated string.

#[no_mangle]
pub unsafe extern "C" fn process(db_path: *mut c_char, db_data: *mut c_char) {
	// Get the filenames 
	let temp_str = {
        assert!(!db_data.is_null());	
		CStr::from_ptr(db_data)
	};
    let temp_path = {
        assert!(!db_path.is_null());
        CStr::from_ptr(db_path)
    };

    let mut read_filename = String::from(temp_path.to_str().unwrap());
    let mut write_filename = String::from(temp_path.to_str().unwrap());
    read_filename.push_str(temp_str.to_str().unwrap());
    write_filename.push_str("calc_data.json");

    process_local(read_filename, write_filename)
}

fn process_local(read_filename: String, write_filename: String) {
	// Parse Stored Data (formatted as a Phylogenetic Set)
    let v: PhylogeneticSet;
    match File::open(&read_filename) {
        Ok(file) => {
            let buf_reader = BufReader::new(file);
            v = match serde_json::from_reader(buf_reader) {
                Ok(data) => data,
                Err(e) => panic!("Could Not Parse Stored Data, Exiting: {}", e)
            };
        },
        Err(e) => panic!("Could Not Open Stored Data, Exiting: {}", e)
    };

    println!("Opened File.");
	let mut r: HashMap<String, ProcessedData> = HashMap::new();
	for (species, data) in v.species {
		// Copy over the pathway and genetic information.
		r.insert(species.clone(), ProcessedData::from(data.pathways, data.genes));

		// Load up the PAML file data to get selection information and related and ancestral sequences.
		for (gene_id, files) in data.files.iter() {
			match PAML::read(&files.paml) {
				Ok(mut paml_data) => {
					if paml_data.append(&str::replace(&files.paml, ".subTree_tree_1.paml_rooted", "_1.RST")).is_some() {
                        paml_data.regenerate_records();
                        
                        // Looking at the Selection Records.
                        for record in paml_data.records {
                            if let Some(val) = r.get_mut(&species) {
                                if record.name == format!("T{}", gene_id) {
                                    // Record contains a ID linkable to TAED Data.
                                    if !val.selection.contains_key(&files.gene_name) { 
                                        val.selection.insert(files.gene_name.clone(), Vec::new());
                                    }
                                    val.selection.get_mut(&files.gene_name).unwrap().push(record.dNdS());
                                    val.ancestor_pairs.insert(gene_id.clone(), (record.id, record.ancestor));			
                                } else if record.name.len() > 3 {
                                    val.related_selection.insert(record.name[1..].to_string(), record.dNdS());
                                    val.ancestor_pairs.insert(gene_id.clone(), (record.id, record.ancestor));
                                } else if record.name.parse::<u32>().is_ok() {
                                    val.ancestor_pairs.insert(gene_id.clone(), (record.id, record.ancestor));
                                }

                            }
                        }
    				}
                },
				Err(e) => {
                    println!("Error {} Reading PAML Data {}", e, str::replace(&files.paml, ".subTree_tree_1.paml_rooted", "_1.RST"));
				}
			}

			// Get phylogenetic related species from associated tree file.
			match Clade::read(&files.tree, "newick") {
				Ok(tree) => {
                    if let Some(val) = r.get_mut(&species) {
                        //if !val.trees.contains_key(&files.gene_name) { val.trees.insert(files.gene_name.clone(), Vec::new()); }
                        val.related_arctic_species.append(&mut get_matching_entries(&tree, &v.species_list));
                        val.related_nonarctic_species.append(&mut get_divergent_entries(&tree, &v.species_list));
                        //val.trees.get_mut(&files.gene_name).unwrap().push(tree);
                    }					
				},
				Err(e) => {
                    println!("Error {} Reading Newick Data {}", e, &files.tree);
				}
			}
		}	
	}
	
	// Writes analytics data to file.
    match File::create(&write_filename) {
        Ok(mut datastore) => 
            match serde_json::to_string(&r) { 
                Ok(data_dump) => match datastore.write_all(data_dump.as_bytes()) 
                    { Ok(()) => (), Err(e) => eprintln!("Could Write Serialized Analytics Data {}", e) }, 
                Err(e) => eprintln!("Could Not Serialize Analytics Data {}", e) 
            }, 
        Err(e) => eprintln!("Could Create File For Pre-Calculated Data: {}", e)
    };
}