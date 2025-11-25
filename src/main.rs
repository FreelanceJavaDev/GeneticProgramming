use std::env;
use std::io::Error;
use std::vec;
mod arg_parser;
use arg_parser::*;
mod dna_compiler;
use dna_compiler::*;
mod executable_header;
use executable_header::*;
mod file_handler;
use file_handler::*;



fn main() -> std::io::Result<()> {
	let main_args: Vec<String> = env::args().collect();
	let (check, tmp_path) = parse_main_args(&main_args)?;
	if !check.check_files() { return Ok(()); }
	if !check.is_req_valid() { Error::new(std::io::ErrorKind::InvalidData, "Parameters are invalid!"); }

	let mut params: Params = Params::from(check);
	let table: [u8; ASCII_TABLE_SIZE] = gen_binary_table(params.bp_val_map());
	match params.flags {
		0b00_000|0b00_001 => ls_bp_first(&mut params, &table)?,
		0b00_010|0b00_011 => fixed_left_shift_first(&mut params, &table)?,
		0b00_100|0b00_101 => codon_ls_bp_shift(&mut params, &table)?,
		0b00_110|0b00_111 => codon_bp_fr_shift(&mut params, &table)?,
		_ => unreachable!(),
	}
	cleanup(&mut params, tmp_path.as_str())
}

fn cleanup(params: &mut Params, tmp_path: &str) -> std::io::Result<()> {
	params.sync_files()?;
	generate_executable(params)?;
	return remove_temp_files(tmp_path,vec!["code.bin", "data.bin"]);
}
