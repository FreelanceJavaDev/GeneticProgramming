use std::collections::HashMap;
use std::ffi::{OsString, OsStr};
pub(crate) use std::fs::{File, OpenOptions};
use std::io::{BufReader, BufWriter, Error, ErrorKind};
use std::path::{Path, PathBuf};

const REQ_ARG_NUM: usize = 3;
#[derive(Eq, PartialEq)]
pub enum InputFileFormat {
	Fna,
	Fasta,
	Unsupported,
}

#[derive(Eq, PartialEq, Copy, Clone)]
pub enum EncodingMethod {
	INVALID = -2,
	M1CS = -1, //out = (a_bp & 0x06) >> 1
	M2A0C1,    //A=00, C=01, G=10, T=11
	M2A1C0,    //A=01, C=00, G=00, T=10
	M2A2C3,    //A=10, C=11, G=00, T=01
	M2A3C2,    //A=11, C=10, G=01, T=00
	M3A0C1,    //A=0, C=1, G=1, T=0
	M3A1C0,    //A=1, C=0, G=0, T=1
	M3A0C1XOR, //A=0, C=1, G=0, T=1
	M3A1C0XOR, //A=1, C=0, G=1, T=0

}


#[derive(Eq, PartialEq, Copy, Clone)]
#[repr(u8)]
pub enum BitFlags {
	NoFlags,
	CodeSectionPadding = 0b00_001,
	FixedLeftShift = 0b00_010,
	UseCodons = 0b00_100,
	CodonData = 0b01_000,
	InPlace = 0b10_000
}


/***
 * A main argument parser and validator.
 * This should never be instantiated directly.
 */
pub struct ParamsCK {
    pub reader: Option<BufReader<File>>,
    pub writer: Option<BufWriter<File>>,
	pub code_fp: Option<OsString>,
	pub data_fp: Option<OsString>,
	pub compile_encode: EncodingMethod,
    pub ifile_ext: InputFileFormat,
	pub flags: u8, // b0: 0= lsb offset shifting, 1= simple fixed bit shift, b1=codons?, b2=in-place(1) or contiguous(0) data/code for codons. b3=codon inside = code(0) or data(1)
}
impl ParamsCK {
	fn with_flags_code_only(i_file: File, o_file: File, encode: EncodingMethod, if_ext: InputFileFormat, cond: u8, code_bin: &OsStr) -> ParamsCK {
		ParamsCK {
			reader: Option::Some(BufReader::new(i_file)),
			writer: Option::Some(BufWriter::new(o_file)),
			code_fp: Option::Some(code_bin.to_os_string()),
			data_fp: None,
			compile_encode: encode,
			ifile_ext: if_ext,
			flags: cond
		}
	}

	fn with_flags_code_and_data(i_file: File, o_file: File, encode: EncodingMethod, if_ext: InputFileFormat, cond: u8, code_bin: &OsStr, data_bin: &OsStr) -> ParamsCK {
		ParamsCK {
			reader: Option::Some(BufReader::new(i_file)),
			writer: Option::Some(BufWriter::new(o_file)),
			code_fp: Option::Some(code_bin.to_os_string()),
			data_fp: Option::Some(data_bin.to_os_string()),
			compile_encode: encode,
			ifile_ext: if_ext,
			flags: cond
		}
	}
	/***
	 * Checks if the file reader and file writer actually exist.
	*/
	pub fn check_files(&self) -> bool { self.reader.is_some() && self.writer.is_some() }

	/***
	 * Checks all critical parameters are valid.
	 */
	pub fn is_req_valid(&self) -> bool {
		self.reader.is_some() && self.writer.is_some() && self.compile_encode != EncodingMethod::INVALID && self.ifile_ext != InputFileFormat::Unsupported
	}
}
impl Default for ParamsCK {
	/***
	 * Constructs a default initialized parameter checker.
	 * This is only used when the `--help` or `-h` argument is used with no other arguments.
	 * It is designed to print help information and then terminate.
	 */
    fn default() -> ParamsCK {
        ParamsCK {
            reader: None,
            writer: None,
			code_fp: None,
			data_fp: None,
            compile_encode: EncodingMethod::INVALID,
            ifile_ext: InputFileFormat::Unsupported,
			flags: 0
        }
    }
}

fn check_input_file(in_path: &str) -> std::io::Result<File> {
    if in_path.is_empty() {
        return Err(Error::from(ErrorKind::InvalidFilename));
    }
    let path = Path::new(in_path);
    if !path.exists() || !path.is_file() || path.file_stem().is_none() {
        return Err(Error::from(ErrorKind::InvalidFilename));
    }
    match path.extension().expect("No file extension detected.").to_str().unwrap_or_default(){
        "fna" | "fasta" => return File::open(path),
        _ => { return Err(Error::new(ErrorKind::Unsupported,"Unsupported File format.")); },
    }
}

fn check_output_file(dir_path: &str, out_file: &str, encode: &EncodingMethod) -> std::io::Result<File> {
    if dir_path.is_empty() || out_file.is_empty() {
        return Err(Error::from(ErrorKind::InvalidFilename));
    }
    if out_file.rfind('.').is_none() {
        return Err(Error::new(ErrorKind::InvalidFilename, "No file extension found!"));
    }
    let mut path: PathBuf = PathBuf::from(dir_path);
    if !path.exists() || !path.is_dir() {
        return Err(Error::from(ErrorKind::InvalidFilename));
    }
    let mut t = match encode {
        EncodingMethod::M2A0C1 => "M2A0C1_".to_string(),
        EncodingMethod::M2A1C0 => "M2A1C0_".to_string(),
        EncodingMethod::M2A2C3 => "M2A2C3_".to_string(),
        EncodingMethod::M2A3C2 => "M2A3C2_".to_string(),
		EncodingMethod::M3A0C1 => "M3A0C1_".to_string(),
		EncodingMethod::M3A1C0 => "M3A1C0_".to_string(),
		EncodingMethod::M3A0C1XOR => "M3A0C1XOR_".to_string(),
		EncodingMethod::M3A1C0XOR => "M3A1C0XOR_".to_string(),
        EncodingMethod::M1CS => "M1CS_".to_string(),
        EncodingMethod::INVALID => { return Err(Error::new(ErrorKind::Other, "Invalid encoding method")); }
    };
    t.push_str(out_file);
    path.push(t);

    let file_out: File = (match path.exists() {
        true => OpenOptions::new().write(true).truncate(true).read(true).open(path.as_path()),
        false => OpenOptions::new().write(true).create(true).read(true).open(path.as_path())
	})?;
    Ok(file_out)
}

fn check_flags(flag_list: Vec<&str>, dir_path: &Path) -> std::io::Result<u8> {
	if !dir_path.exists() || !dir_path.is_dir() {
        return Err(Error::from(ErrorKind::NotFound));
    }
	if flag_list.is_empty() { return Ok(0); }
	let mut ret: u8 = BitFlags::NoFlags as u8;
	for flag_str in flag_list {
		ret |= match flag_str {
			"-p"|"--padding" => BitFlags::CodeSectionPadding as u8,
			"-fls"|"--fixed-lshift" => BitFlags::FixedLeftShift as u8,
			"--use-codon"|"-uc3" => BitFlags::UseCodons as u8,
			"--in-place"|"-inp" => BitFlags::InPlace as u8,
			"--codon-data"|"-ctd" => BitFlags::CodonData as u8,
			_ => { return Err(Error::new(ErrorKind::InvalidInput,"Unknown flag provided.")); },
		};
	}

	match ret {
		0b01000|0b10000|0b11000|0b01010|0b10010|0b11010 => { return Err(Error::new(ErrorKind::InvalidData, "In-place and codon data options are only valid with use-codons")); },
		0b01001|0b10001|0b11001|0b01011|0b10011|0b11011 => { return Err(Error::new(ErrorKind::InvalidData, "In-place and codon data options are only valid with use-codons")); },
		0b01100|0b10100|0b11100|0b01110|0b10110|0b11110|0b00110=> { return Err(Error::new(ErrorKind::Unsupported, "Alternate Codon paths are currently unimplemented.")); },
		0b01101|0b10101|0b11101|0b01111|0b10111|0b11111|0b00111=> { return Err(Error::new(ErrorKind::Unsupported, "Alternate Codon paths are currently unimplemented.")); },
		_ => return Ok(ret),
	};
}

pub fn parse_main_args(arg_list: &Vec<String>) -> std::io::Result<(ParamsCK, String)> {
	let mut ret_str: String = String::new();
    if arg_list.len() >= 4 {
        let mut a_flags: Vec<&str> = Vec::new();
        let mut a_params: Vec<&str> = Vec::new();
        for i in 1..arg_list.len() {
            if arg_list[i].starts_with("--") {
                if arg_list[i].contains('=') {
                    a_params.push(arg_list[i].as_str());
                } else {
                    a_flags.push(arg_list[i].as_str());
                }
            } else if arg_list[i].starts_with('-') {
                if arg_list[i].contains('=') {
                    return Err(Error::new(ErrorKind::InvalidInput,"Invalid Input. Only flags start with -"));
                } else {
                    a_flags.push(arg_list[i].as_str());
                }
            } else {
                return Err(Error::new(ErrorKind::InvalidInput,"Invalid Input. All user arguments start with -- or -"));
            }
        }
        let mut param_map: HashMap<usize, &str> = HashMap::new();
        let mut st: &str = "M2A0C1";
        for id_val in 0..a_params.len() {
            let val_parse: Vec<&str> = a_params[id_val].split('=').collect();
            match val_parse[0] {
                "--input" => {
                    param_map.insert(0, val_parse[1]);
                }
                "--outputDir" => {
					ret_str = String::from(val_parse[1]);
                    param_map.insert(1, val_parse[1]);
                }
                "--outputFile" => {
                    param_map.insert(2, val_parse[1]);
                }
                "--encoding" => {
                    st = val_parse[1];
                    st.to_string().make_ascii_uppercase();
                    param_map.insert(3, st);
                }
                _ => {
                    return Err(Error::new(ErrorKind::InvalidInput,"Unknown parameter tag."));
                }
            }
        }
        if param_map.len() < REQ_ARG_NUM {
            return Err(Error::new(ErrorKind::InvalidInput, "Missing required arguments."));
        }
        if !param_map.contains_key(&3) {
            param_map.insert(3, st);
        }
        let encoding: EncodingMethod = match param_map.get(&3).unwrap() {
            &"M2A0C1" => EncodingMethod::M2A0C1,
            &"M2A1C0" => EncodingMethod::M2A1C0,
            &"M2A2C3" => EncodingMethod::M2A2C3,
            &"M2A3C2" => EncodingMethod::M2A3C2,
			&"M3A0C1" => EncodingMethod::M3A0C1,
			&"M3A1C0" => EncodingMethod::M3A1C0,
			&"M3A0C1XOR" => EncodingMethod::M3A0C1XOR,
            &"M3A1C0XOR" => EncodingMethod::M3A1C0XOR,
			&"M1CS" => EncodingMethod::M1CS,
            _ => EncodingMethod::INVALID,
        };
        for i in 0..REQ_ARG_NUM {
            if !param_map.contains_key(&i) {
                return Err(Error::new(ErrorKind::InvalidInput,"Missing required arguments."));
            }
        }
        let in_file: File = check_input_file(param_map[&0])?;
        let out_file: File = check_output_file(param_map[&1], param_map[&2], &encoding)?;
        let path_in = Path::new(param_map[&0]);
        let in_format: InputFileFormat = match path_in.extension().expect("No file extension detected.").to_str().unwrap_or_default()
        {
            "fna" => InputFileFormat::Fna,
            "fasta" => InputFileFormat::Fasta,
            _ => InputFileFormat::Unsupported,
        };

		let flags: u8 = check_flags(a_flags,Path::new(&param_map[&(1 as usize)]))?;
		let mut c_path: PathBuf = PathBuf::from(param_map[&1]);
		c_path.push("code.bin");
		let code_fp =  c_path.as_os_str();
		if (flags & BitFlags::UseCodons as u8) > 0 {
			let mut d_path: PathBuf = PathBuf::from(param_map[&1]);
			d_path.push("data.bin");
			let data_fp =  d_path.as_os_str();
			Ok((ParamsCK::with_flags_code_and_data(in_file, out_file, encoding, in_format, flags, code_fp, data_fp), ret_str))
		}
		else {
        	Ok((ParamsCK::with_flags_code_only(in_file, out_file, encoding, in_format, flags, code_fp), ret_str))
		}
    }
	else if arg_list.len() == 2 {
        match arg_list[1].as_str() {
            "-h" | "--help" => println!("Usage: <prog_path> --input=<path to fastfa file> --outDir=<output file directory path> --outFile=<file name>"),
            _ => {
                return Err(Error::new(ErrorKind::InvalidInput,"Only valid 2 argument rin is for help function."));
            }
        }
        Ok((ParamsCK::default(), ret_str))
    } else {
        return Err(Error::new(ErrorKind::InvalidInput,"Not enough arguments! Use --help for details."));
    }
}