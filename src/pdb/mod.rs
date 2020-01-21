///! Tools for fetch, parsing, and converting PDB Data.
use std::{io, fs::File, io::{BufRead, BufReader, SeekFrom}};
use std::io::prelude::*;
use std::collections::HashMap;
use std::rc::Rc;
use std::cell::RefCell;
use std::str;
use arrayref::array_ref;

pub type PDBID = [u8; 4];
pub struct Pos(pub usize, pub char);
pub type CloseContact = (Rc<RefCell<Atom>>, Rc<RefCell<Atom>>);

impl Pos {
	pub fn index (&self) -> usize {
		// Sequence # starts at 1; insertion codes start at 'A'
		// To flatten into an index, get insertion code, plus number of chars past A, minus 1
		self.0 + if self.1 == ' ' { 0 } else { (self.1 as u8 - b'A') as usize } - 1
	}
}

/// Journal information that can be stored in PDB
pub struct Journal {
	pub authors: Vec<String>,
	pub title: Vec<String>,
	pub editors: Vec<String>,
	pub reference: String,
	pub publisher: String,
	pub refn: Vec<String>,
	pub pubmed_id: usize,
	pub doi: String
}

/// Database info for a protein record.
pub struct DBInfo {
	pub id: PDBID,
	pub chain: char,
	pub db_name: String,
	pub db_accession: String,
	pub db_id: String,
	pub start: Pos,
	pub end: Pos,
	pub db_start: Pos,
	pub db_end: Pos
}


pub struct Compound { }

pub struct Source { }

pub struct Molecule { 
	pub compound: Compound,
	pub source: Source
}

/// All Atom Data Handled by PDB.
pub struct Atom {
	pub serial: usize,
	pub name: [u8; 4],
	pub altloc: char,
	pub coord: Coordinates,
	pub occupancy: f64,
	pub temp: f64,
	pub element: [u8; 2],
	pub charge: [u8; 2],
	pub connections: Vec<Rc<RefCell<Atom>>>,
	pub links: Vec<Rc<Connectivity>>,
    pub anisotropic: [u16; 6]
}

/// XYZ Coordinates used to locate Atoms
pub struct Coordinates {
	pub x: f64,
	pub y: f64,
	pub z: f64
}

/// 
pub struct Strand {
	pub sense: bool, // True for parallel, False for antiparallel
	pub bonded: Rc<RefCell<Atom>>,
	pub remote: Rc<RefCell<Atom>>
}

pub enum SecondaryStructure {
	Helix(u32, String, Pos, Pos),
	Sheet(String, Vec<(Pos, Pos, Option<Strand>)>),
	SSBond((Rc<RefCell<Residue>>, usize), (Rc<RefCell<Residue>>, usize), f64)
}

pub enum Connectivity {
	CisPep(Rc<RefCell<Residue>>, Rc<RefCell<Residue>>, f64),
	Link((Rc<RefCell<Atom>>, usize), (Rc<RefCell<Atom>>, usize), f64),
	SSBond((Rc<RefCell<Residue>>, usize), (Rc<RefCell<Residue>>, usize), f64)
}

pub struct Residue { 
	pub name: [u8; 3],
	pub atoms: Vec<Rc<RefCell<Atom>>>,
	pub chain: char,
	pub seq: usize,
	pub insertion_code: char,
	pub connectivity: Vec<Rc<Connectivity>>,
	pub sites: Vec<Rc<Site>>,
	pub ss: Option<SecondaryStructure>
}

pub struct Model {
	pub kind: String,
	pub residues: Vec<Rc<RefCell<Residue>>>,
}

impl Model {
    pub fn bare() -> Model {
        Model {
            kind: String::new(),
            residues: Vec::new()
        }
    }
}


pub struct Site {
	pub name: String,
	pub residues: Vec<Rc<RefCell<Residue>>>,
	pub detail: Rc<Remark>
}


pub struct Group {
	pub name: String,
	pub chain: Vec<char>,
	pub components: Vec<String>,
	pub description: String
}


pub struct CoordinationAngleSet(Rc<RefCell<Atom>>, Rc<RefCell<Atom>>, f64);


pub enum Remark {
    Refinement(String),
	Biomolecule(String),
	Transformations(String),
	SpecialPosition(Vec<Rc<RefCell<Atom>>>),
	Compound(String, Vec<Group>),
	SourceDetails(String),
	MissingResidues(Vec<Residue>),
	MissingAtoms(Vec<Atom>),
	ZeroOccupancyResidues(Vec<Rc<RefCell<Residue>>>),
	ZeroOccupancyAtoms(Vec<Rc<RefCell<Atom>>>),
	CloseContacts(Vec<CloseContact>, String),
	CloseContactsSameUnit(Vec<CloseContact>, String),
	/* For Later Implementation
	UnusualCovalentBondLengths(HashMap<Rc<RefCell<Residue>>, (usize, usize)>, String),
	UnusualCovalentBondAngles(HashMap<Rc<RefCell<Residue>>, (usize, usize)>, String),
	UnusualTorsionAngles(HashMap<Rc<RefCell<Residue>>, (f64, f64)>, String),
	UnusualConformation(Vec<(Rc<RefCell<Residue>>, Rc<RefCell<Residue>>, f64)>, String),
	UnusualPlanarGroups(HashMap<Rc<RefCell<Residue>>, (f64, String)>, String),
	UnusualMainChainPlanarity(HashMap<Rc<RefCell<Residue>>, f64>, String),
	UnusualChiralCenters(HashMap<Rc<RefCell<Residue>>, (f64, char, char, String)>, String),
	*/
	DistantSolventAtoms(Vec<Rc<RefCell<Atom>>>),
	Heterogens(String),
	NonPolymerMissingAtoms(Vec<Residue>),
	NonPolymerZeroOccupancy(Vec<Residue>),
	//MetalCoordination(HashMap<Rc<RefCell<Residue>>, CoordinationAngleSet>),
	Inhibitor(String),
	HelixDetails(String),
	SheetDetails(String),
	SiteDetails(String),
	RelatedEntries(String, HashMap<PDBID, (String, String)>),
	Sequence(String),

    Error(String)
}
impl Remark {
	pub fn code(&self) -> u32 {
		match self {
            Remark::Refinement(_) => 3,
			Remark::Biomolecule(_) => 300,
			Remark::Transformations(_) => 350,
			Remark::SpecialPosition(_) => 375,
			Remark::Compound(_,_) => 400,
			Remark::SourceDetails(_) => 450,
			Remark::MissingResidues(_) => 465,
			Remark::MissingAtoms(_) => 470,
			Remark::ZeroOccupancyResidues(_) => 475,
			Remark::ZeroOccupancyAtoms(_) => 480,

			Remark::CloseContacts(_, _) | Remark::CloseContactsSameUnit(_, _) | 
/*			Remark::UnusualCovalentBondLengths(_, _) | 
			Remark::UnusualCovalentBondAngles(_, _) | Remark::UnusualTorsionAngles(_, _) | 
			Remark::UnusualConformation(_, _) | Remark::UnusualPlanarGroups(_, _) | 
			Remark::UnusualMainChainPlanarity(_, _) | 
			Remark::UnusualChiralCenters(_, _) => 500,
*/
			Remark::DistantSolventAtoms(_) => 525,
			Remark::Heterogens(_) => 600,
			Remark::NonPolymerMissingAtoms(_) => 610,
			Remark::NonPolymerZeroOccupancy(_) => 615,
//			Remark::MetalCoordination(_) => 620,
			Remark::Inhibitor(_) => 630,
			Remark::HelixDetails(_) => 650,
			Remark::SheetDetails(_) => 700,
			Remark::SiteDetails(_) => 800,
			Remark::RelatedEntries(_,_) => 900,
			Remark::Sequence(_) => 999,

            Remark::Error(_) => 0
		}
	}

    pub fn parse_remark<U: Read + Seek>(reader: &mut BufReader<U>, remark_num: u32) -> Result<Remark, io::Error> {
        let mut line = String::new();
        reader.read_line(&mut line)?;
        let mut old_remark_num: u32 = remark_num;
		#[allow(unused_mut)]
        let mut remark = Remark::Error("No Valid Remark Section".to_string());
		
        while &line[0..6] == "REMARK" && (old_remark_num == remark_num)  {
            old_remark_num = match &line[7..10].trim_start().parse::<u32>() {
                Ok(remark_num) => *remark_num,
                Err(_e) => { 
					/* Something*/ 
					0
				}
            };
			line = String::new();
            reader.read_line(&mut line)?;
			print!(":{}", line);
        }
        reader.seek(SeekFrom::Current(-(line.len() as i64)))?;

        Ok(remark)
    }
}


pub struct PDB {
	// PDB Metadata Information
	pub id: PDBID,
	pub date: String,
	pub title: String,
	pub authors: Vec<String>,
	pub classification: String,
	pub caveat: String,
	pub keywords: Vec<String>,
	pub experimental_data: Vec<String>,
	pub related: Vec<PDBID>,
	pub model_count: usize,
	pub model_type: String,
	pub revision_history: String,
	pub journal: String,

    pub atoms: Vec<Rc<RefCell<Atom>>>,
	pub residues: Vec<Rc<RefCell<Residue>>>,
	pub structures: Vec<Rc<SecondaryStructure>>,
	pub connectivity: Vec<Rc<Connectivity>>,

	// Structural Models
	pub models: Vec<Model>,
	pub sites: Vec<Rc<Site>>,

	// Crosslinks to Database Info

	// Database
	pub db_ref: Vec<DBInfo>,

	// Remarks - Comments & Autogenerated Info
	pub remarks: Vec<Remark>,

	// PDB Molecules
	pub molecules: Vec<Molecule>
}

impl PDB {
	const EMPTYID: [u8; 4] = [b' ', b' ', b' ', b' '];

	pub fn bare() -> PDB {
		PDB {
			id: PDB::EMPTYID,
			date: String::new(),
			title: String::new(),
			authors: Vec::new(),
			classification: String::new(),
			caveat: String::new(),
			keywords: Vec::new(),
			experimental_data: Vec::new(),
			related: Vec::new(),
			revision_history: String::new(),
			models: Vec::new(),
			sites: Vec::new(),
			remarks: Vec::new(),
			molecules: Vec::new(),

			db_ref: Vec::new(),

            atoms: Vec::new(),
			residues: Vec::new(),
			structures: Vec::new(),
			connectivity: Vec::new(),

            model_count: 0,
            model_type: String::new(),

			journal: String::new(),
		}
	}

	pub fn find_atom_in_residue (&self, res_seq: usize, atom_name: [u8; 4]) -> Option<Rc<RefCell<Atom>>> {
		for i in 0..self.residues[res_seq - 1].borrow().atoms.len() {
			if self.residues[res_seq - 1].borrow().atoms[i].borrow().name == atom_name {
				return Some(Rc::clone(&self.residues[res_seq - 1].borrow().atoms[i]));
			}
		}
		None
	}

	pub fn read(file_name: &str) -> Result<PDB, io::Error> {
		let mut data: PDB = PDB::bare();
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
		fn atov<T: str::FromStr> (a: &[u8]) -> T 
			where <T as std::str::FromStr>::Err: std::fmt::Debug { 
			str::from_utf8(a).unwrap().trim_start().parse().unwrap() 
		};
		
        // Push in a default model for PDB files with only 1 model (where MODEL sections are not required).
        data.models.push(Model {
            kind: String::new(),
            residues: Vec::new()
        });

		match File::open(file_name) {
			Ok(file) => {
				let mut buf_reader = BufReader::new(file);
				let mut file_data: Vec<u8> = Vec::new();
                let mut cur_model: usize;
				let mut sheet_buffer: Vec<Vec<u8>> = Vec::new();
				let mut link_buffer: Vec<Vec<u8>> = Vec::new();
				let mut line_size: usize = buf_reader.read_until(b'\n', &mut file_data)?;

				while line_size > 4 {
					let line = file_data[(file_data.len() - line_size)..].to_vec();
					
					print!("?{}:{}", line_size, str::from_utf8(&line).unwrap());
					match str::from_utf8(&line[0..6]).unwrap() {
						// Skipped Title Sections: TODO
						"OBSLTE" | "COMPND" | "SOURCE" | "REVDAT" | "SPRSDE" | "JRNL  " => {},

						// Title Sections
						"HEADER" => {
							data.classification = str::from_utf8(&line[10..50]).unwrap().trim().to_string();
							// data.date = ???chrono?? line[50..59] TODO
							data.id[0..4].clone_from_slice(&line[62..66]) //  for i in 0..4 {	data.id[i] = line[62 + i]; }
						},
						"TITLE " => data.title.push_str(str::from_utf8(&line[10..line_size]).unwrap().trim()),
						"SPLIT " => {
							let mut start: usize = 11;
							while line[start..start + 4] != PDB::EMPTYID {
								data.related.push(array_ref![line, start, 4].clone());
								start += 5;
								if start + 4 >= line_size { break; }
							}
						},
						"CAVEAT" =>	data.caveat.push_str(str::from_utf8(&line[19..line_size]).unwrap().trim()),
						"KEYWDS" => parse_separator_list(&mut data.keywords, ',', str::from_utf8(&line[10..line_size]).unwrap()),
						"EXPDTA" => parse_separator_list(&mut data.experimental_data, ';', str::from_utf8(&line[10..line_size]).unwrap()),
                        "NUMMDL" => data.model_count = atov(&line[10..14]),
                        "MDLTYP" => data.model_type.push_str(str::from_utf8(&line[10..line_size]).unwrap()),
                        "AUTHOR" => parse_separator_list(&mut data.authors, ',', str::from_utf8(&line[10..line_size]).unwrap()), // Could ignore breaks after a comma for faster read I think

						// REMARK Section: TODO
						"REMARK" => {
                            buf_reader.seek(SeekFrom::Current(-(line_size as i64)))?;
                            match Remark::parse_remark(&mut buf_reader, str::from_utf8(&line[7..10]).unwrap().trim_start().parse::<u32>().unwrap()) {
                                Ok(new_remark) => data.remarks.push(new_remark),
                                Err(e) => { return Err(io::Error::new(io::ErrorKind::InvalidData, format!("Invalid REMARK Section {}", e))); }
                            }
                        },

						// Skipped Primary Structure Sections
						"SEQADV" | "MODRES" => {}, // TODO

						// Primary Structure Sections
						"DBREF " =>
                            data.db_ref.push(DBInfo {
								id: *array_ref!(line, 7, 4),
								chain: line[13] as char,
								db_name: str::from_utf8(&line[26..32]).unwrap().to_string(),
								db_accession: str::from_utf8(&line[33..41]).unwrap().to_string(),
								db_id: str::from_utf8(&line[42..54]).unwrap().to_string(),
								start: Pos(atov(&line[14..18]), line[18] as char),
								end: Pos(atov(&line[20..24]), line[24] as char),
								db_start: Pos(atov(&line[55..60]), line[60] as char),
								db_end: Pos(atov(&line[62..67]), line[67] as char),
							}),
                        "DBREF1" =>
                            data.db_ref.push(DBInfo {
								id: *array_ref!(line, 7, 4),
								chain: line[13] as char,
								db_name: str::from_utf8(&line[26..32]).unwrap().to_string(),
								db_accession: String::new(),
								db_id: str::from_utf8(&line[47..67]).unwrap().to_string(),
								start: Pos(atov(&line[14..18]), line[18] as char),
								end: Pos(atov(&line[20..24]), line[24] as char),
								db_start: Pos(0, ' '),
								db_end: Pos(0, ' '),
							}),
                        "DBREF2" => {
							
                            match data.db_ref.last_mut() {
								Some(db_rec) => {
									db_rec.db_accession = str::from_utf8(&line[18..40]).unwrap().to_string();
									db_rec.db_start = Pos(atov(&line[55..60]), line[60] as char);
									db_rec.db_end = Pos(atov(&line[62..67]), line[67] as char);
								},
								None => { 
									return Err(io::Error::new(io::ErrorKind::InvalidData, "Missing DBREF2 Section"));
								}
							}
                        },
						"SEQRES" => {
							let mut start: usize = 19;
							while line[start..start + 3] != [b' ', b' ', b' '] {
								let res = Rc::new(RefCell::new(Residue {
									name: *array_ref![line, start, 3],
									atoms: Vec::new(),
									chain: line[11] as char,
									seq: data.residues.len() + 1,
									insertion_code: ' ',
									ss: None,
									connectivity: Vec::new(),
									sites: Vec::new(),
								}));
								data.residues.push(res);
								start += 4;
								if start + 3 >= line_size { break; }
							}
						},

						// Skipped Heterogen Sections TODO
						"HET   " | "HETNAM" | "HETSYN" | "FORMUL" => {},

						// Secondary Structure Sections.
                        "HELIX " =>
							data.structures.push(Rc::new(SecondaryStructure::Helix(
								atov(&line[7..10]),
								str::from_utf8(&line[40..70]).unwrap().to_string(),
								Pos(atov(&line[21..25]), line[25] as char),
								Pos(atov(&line[33..37]), line[37] as char),
							))),
                        "SHEET " => sheet_buffer.push(line), // We buffer it so we can link up with atoms later.
                        "SSBOND" => {
							let first = &data.residues[Pos(atov(&line[17..21]), line[21] as char).index()];
							let second = &data.residues[Pos(atov(&line[31..35]), line[35] as char).index()];

							let connection = Rc::new(Connectivity::SSBond(
								(Rc::clone(first), atov(&line[59..65])),
								(Rc::clone(second), atov(&line[66..72])),
								atov(&line[73..78])
							));
							first.borrow_mut().connectivity.push(Rc::clone(&connection));
							second.borrow_mut().connectivity.push(Rc::clone(&connection));
							data.connectivity.push(connection);
						},
						"LINK  " => link_buffer.push(line), // Another Buffer for Atom Connectivity
						"CISPEP" =>  {
							let first = &data.residues[Pos(atov(&line[17..21]), line[21] as char).index()];
							let second = &data.residues[Pos(atov(&line[31..35]), line[35] as char).index()];

							let connection = Rc::new(Connectivity::CisPep(
								Rc::clone(first),
								Rc::clone(second),
								atov(&line[53..59])
							));
							first.borrow_mut().connectivity.push(Rc::clone(&connection));
							second.borrow_mut().connectivity.push(Rc::clone(&connection));
							data.connectivity.push(connection);
						},

						"SITE  " => {
							let mut temp: Rc<Site>;
							let site = if atov::<u32>(&line[7..10]) > 1 {
								data.sites.last_mut().unwrap()
							} else {
								temp = Rc::new(Site {
									name: str::from_utf8(&line[11..14]).unwrap().to_string(),
									detail: Rc::new(Remark::SiteDetails("Empty".to_string())),
									residues: Vec::new(),
								});
								&mut temp
							};

							let mut start = 23;
							while line[start..start + 4] != PDB::EMPTYID {
								let residue = &data.residues[Pos(atov(&line[start..start + 4]), line[start + 4] as char).index()];
								Rc::get_mut(site).unwrap().residues.push(
									Rc::clone(residue)
								);
								residue.borrow_mut().sites.push(Rc::clone(&site));
								start += 6;
								if start + 4 >= line_size { break; }
							}
						},

						// Not Yet Implemented (Crystalline Structure) TODO
						"CRYST1" | "ORIGX1" | "ORIGX2" | "ORIGX3" | 
						"SCALE1" | "SCALE2" | "SCALE3" |
						"MTRIX1" | "MTRIX2" | "MTRIX3" => {},

						"MODEL " => {
                            cur_model = atov(&line[10..14]);
                            if data.models.len() < cur_model {
                                data.models.push(Model::bare());
                            }
                        },

                        "ATOM  " => {
                            let atom = Atom {
                                serial: atov(&line[6..11]),
                                name: *array_ref!(line, 12, 4),
                                altloc: line[16] as char,
                                coord: Coordinates { 
                                    x: atov(&line[30..38]), 
                                    y: atov(&line[38..46]), 
                                    z: atov(&line[46..54]),                                
                                },
                                occupancy: atov(&line[54..60]),
                                temp: atov(&line[60..66]),
                                element: *array_ref!(line, 76, 2),
                                charge: if line_size >= 80 {
									*array_ref!(line, 78, 2)
								} else { [b' ', b' '] },
                                connections: Vec::new(),
								links: Vec::new(),
                                anisotropic: [0; 6]
                            };
                            let residue: usize = Pos(atov(&line[22..26]), line[26] as char).index();
                            let atom_ref = Rc::new(RefCell::new(atom));
                            // Atom is added to base object for ease of reference.
                            data.atoms.push(Rc::clone(&atom_ref));
                            // Atom is also added to the specific residue.
							// Don't need a second push since the underlying (residue) data is already updated.
                            data.residues[residue].borrow_mut().atoms.push(atom_ref);
							//data.models[cur_model].residues[residue].borrow_mut().atoms.push(Rc::clone(&atom_ref));
                        },

                        "ANISOU" => {
							let atom: usize =  atov(&line[6..11]);
                            data.atoms[atom].borrow_mut().anisotropic = [
                                atov(&line[28..35]),
                                atov(&line[35..42]),
                                atov(&line[42..49]),
                                atov(&line[49..56]),
                                atov(&line[56..63]),
                                atov(&line[63..70]),
                            ];
                            
                        },

                        "CONECT" => {
                            let mut start: usize = 11;
							let atom: usize = atov(&line[6..11]);
							while line[start..start + 4] != PDB::EMPTYID {
                                let conn: usize = atov(&line[start..start + 4]);
                                // Increment reference count 
                                let connected = Rc::clone(&data.atoms[conn]);
                                data.atoms[atom].borrow_mut().connections.push(connected); // Damn Borrow Checker.  Fun with RefCell Next.
								start += 5;
								if start + 4 >= line_size { break; }
							}
                        },

                        "HETATM" => {}, // TODO

                        "TER   " | "ENDMDL" => {}, // Do Nothing - We Don't Need Them the Way We're Parsing

						"MASTER" => {} // TODO - Use for Verification
						_ => return Err(io::Error::new(io::ErrorKind::InvalidData, format!("Invalid Section {} Found", str::from_utf8(&line[0..6]).unwrap())))
					}
					line_size = buf_reader.read_until(b'\n', &mut file_data)?;
				}
				if str::from_utf8(&file_data[file_data.len() - 4..file_data.len() - 1]).unwrap() == "END" {
					// Get through Sheet Buffer
					for line in sheet_buffer {
						print!("{}", str::from_utf8(&line).unwrap());
						let sheet_id = str::from_utf8(&line[11..14]).unwrap().trim_start().to_string();
						let strand = (
							Pos(atov(&line[22..26]), line[26] as char),
							Pos(atov(&line[33..37]), line[37] as char),
							match atov::<i32>(&line[38..40]) { 
								0 => None, // First thing in a strand has no residue links.
								x => if line_size > 41 {  // Sometimes the others don't either (??)
										Some(Strand {
											sense: x > 0,
											bonded: data.find_atom_in_residue(
													Pos(atov(&line[50..54]), line[54] as char).index(), 
													*array_ref!(&line, 41, 4)
												).unwrap(),
											remote: data.find_atom_in_residue(
													Pos(atov(&line[50..54]), line[54] as char).index(), 
													*array_ref!(&line, 41, 4)
												).unwrap(),
										})
									} else { None }
							}
						);
						if let Some(mut maybe_ss) = data.structures.last_mut() {
							if let SecondaryStructure::Sheet(old_name, strands) = Rc::get_mut(&mut maybe_ss).unwrap() {
								if *old_name == sheet_id {
									strands.push(strand);
								} else {
									data.structures.push(Rc::new(SecondaryStructure::Sheet(sheet_id, vec![strand])));
								}
							} else {
								data.structures.push(Rc::new(SecondaryStructure::Sheet(sheet_id, vec![strand])));
							}
						} else { data.structures.push(Rc::new(SecondaryStructure::Sheet(sheet_id, vec![strand]))); }
					}

					// Get Through Link Buffer
					for line in link_buffer {
						print!("{}", str::from_utf8(&line).unwrap());
						let first = data.find_atom_in_residue(Pos(atov(&line[22..26]), line[26] as char).index(), *array_ref!(&line, 12, 4)).unwrap();
						let second = data.find_atom_in_residue(Pos(atov(&line[52..56]), line[56] as char).index(), *array_ref!(&line, 42, 4)).unwrap();
						let connection = Rc::new(Connectivity::Link(
							(Rc::clone(&first), atov(&line[59..65])),
							(Rc::clone(&second), atov(&line[66..72])),
							atov(&line[73..78])
						));
						first.borrow_mut().links.push(Rc::clone(&connection));
						second.borrow_mut().links.push(Rc::clone(&connection));
						data.connectivity.push(connection);
					}

					// Assign Secondary Structures to Atoms
					// Update Secondary Structure Atom References
					println!("Found End.");
				}

				Ok(data)
			}, Err(e) => Err(e)
		}
	}
}

// Body Mass Table N. Lartillot