use std::collections::HashMap;
use std::fs::File;
use std::fs::OpenOptions;
// use std::io::BufRead;
use std::io::BufReader;
use std::io::BufWriter;
// use std::io::Write;
use std::io::{Error, ErrorKind};
use std::path::Path;
use std::path::PathBuf;
const REQ_ARG_NUM: usize = 3;
pub enum InputFileFormat {
	Fna,
	Fasta,
	Unsupported,
}

pub enum EncodingMethod {
	INVALID = -2,
	M1CS = -1, //out = (a_bp & 0x06) >> 1
	M2A0C1,    //A=00, C=01, G=10, T=11
	M2A1C0,    //A=01, C=00, G=00, T=10
	M2A2C3,    //A=10, C=11, G=00, T=01
	M2A3C2,    //A=11, C=10, G=01, T=00
}
// #[repr(usize)]
// #[derive(Debug, PartialEq, Eq, Clone, Copy)]
// pub enum BaseValIndex {
// 	A,
// 	C,
// 	G,
// 	T,
// 	NumBases,
// }

pub struct ParamsCK {
    pub reader: Option<BufReader<File>>,
    pub writer: Option<BufWriter<File>>,
	pub compile_encode: EncodingMethod,
    pub ifile_ext: InputFileFormat,
}
impl ParamsCK {
    pub fn new(i_file: File, o_file: File,encode: EncodingMethod,if_ext: InputFileFormat) -> ParamsCK {
        ParamsCK {
            reader: Option::Some(BufReader::new(i_file)),
            writer: Option::Some(BufWriter::new(o_file)),
            compile_encode: encode,
            ifile_ext: if_ext,
        }
    }
}
impl Default for ParamsCK {
    fn default() -> ParamsCK {
        ParamsCK {
            reader: None,
            writer: None,
            compile_encode: EncodingMethod::INVALID,
            ifile_ext: InputFileFormat::Unsupported,
        }
    }
}

pub fn check_input_file(in_path: &str) -> std::io::Result<File> {
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

pub fn check_output_file(dir_path: &str, out_file: &str, encode: &EncodingMethod) -> std::io::Result<File> {
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

pub fn parse_main_args(arg_list: &Vec<String>) -> std::io::Result<ParamsCK> {
    if arg_list.len() >= 4 {
        let mut a_flags: Vec<String> = Vec::new();
        let mut a_params: Vec<String> = Vec::new();
        for i in 1..arg_list.len() {
            if arg_list[i].starts_with("--") {
                if arg_list[i].contains('=') {
                    a_params.push(arg_list[i].clone());
                } else {
                    a_flags.push(arg_list[i].clone());
                }
            } else if arg_list[i].starts_with('-') {
                if arg_list[i].contains('=') {
                    return Err(Error::new(ErrorKind::InvalidInput,"Invalid Input. Only flags start with -"));
                } else {
                    a_flags.push(arg_list[i].clone());
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

        Ok(ParamsCK::new(in_file, out_file, encoding, in_format))
    } else if arg_list.len() == 2 {
        match arg_list[1].as_str() {
            "-h" | "--help" => println!(
                "Usage: <prog_path> --input=<path to fastfa file> --outDir=<output file directory path> --outFile=<file name>"
            ),
            _ => {
                return Err(Error::new(ErrorKind::InvalidInput,"Only valid 2 argument rin is for help function."));
            }
        }
        Ok(ParamsCK::default())
    } else {
        return Err(Error::new(ErrorKind::InvalidInput,"Not enough arguments! Use --help for details."));
    }
}