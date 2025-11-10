use std::collections::HashMap;
use std::env;
use std::fs::File;
use std::io::BufRead;
use std::fs::OpenOptions;
use std::io::BufReader;
use std::io::BufWriter;
use std::io::Error;
use std::io::Read;
use std::io::Seek;
use std::io::Write;
use std::path::PathBuf;
use std::vec;
mod arg_parser;
mod executable_header;
use executable_header::*;
#[cfg(target_os="windows")]
use executable_header::win::*;
use arg_parser::*;

#[cfg(target_os="linux")]
use crate::executable_header::linux::*;

// #[cfg(target_os ="windows")]
// const EOL: String = "\r\n";
// #[cfg(not(target_os = "windows"))]
// const EOL: String = "\n";

//TODO: bool codon

pub struct CodonOnly {
	pub codon_map: HashMap<u8, u8>,
	pub bp_map: [u8; 256],
	stop_codons: [u8; 3],
	start_codon: u8,
	pub bp_shift: u8,
	pub read_frame: u8,
	pub read_frame_used: u8,
	pub tmp: u8,
	pub tmp_bits_used: u8,
	pub codon_buff: u16,
	pub codon_buff_used: u16,
	pub flags: u8,
	pub inside_codon: bool
}

impl CodonOnly {
	pub fn new(table: &[u8; 256], encoding: EncodingMethod, b_flags: u8) ->Self {
		CodonOnly { codon_map: HashMap::with_capacity(64), bp_map: *table, stop_codons: [0; 3], start_codon: 0,
			bp_shift: (match encoding {
		EncodingMethod::M3A0C1 | EncodingMethod::M3A1C0 | EncodingMethod::M3A0C1XOR | EncodingMethod::M3A1C0XOR => 1,
		_ => 2,
	}), read_frame: 0, read_frame_used: 0, tmp: 0, tmp_bits_used: 0, codon_buff: 0, codon_buff_used: 0, flags: b_flags, inside_codon: false }
	}

	pub fn dyn_init(&mut self) {
		const BP_MAP: [usize; 4] = ['A' as usize, 'C' as usize, 'G' as usize, 'T' as usize];
    	const BP_MAP_STR: [u8; 4] = [b'A', b'C', b'G', b'T'];
		let bp_shifts: [u8; 2] = [1*self.bp_shift, 2*self.bp_shift];
		let mut str_codon: HashMap<String, u8> = HashMap::with_capacity(64);
		match self.flags & 0x1 {
			1 =>  {
				self.start_codon = (self.bp_map[BP_MAP[0]] << bp_shifts[1]) + (self.bp_map[BP_MAP[3]] << bp_shifts[0]) + self.bp_map[BP_MAP[2]]; // ATG
				self.stop_codons = [
					(self.bp_map[BP_MAP[3]] << bp_shifts[1]) + (self.bp_map[BP_MAP[0]] << bp_shifts[0]) + self.bp_map[BP_MAP[0]], //TAA
					(self.bp_map[BP_MAP[3]] << bp_shifts[1]) + (self.bp_map[BP_MAP[2]] << bp_shifts[0]) + self.bp_map[BP_MAP[0]], //TGA
					(self.bp_map[BP_MAP[3]] << bp_shifts[1]) + (self.bp_map[BP_MAP[0]] << bp_shifts[0]) + self.bp_map[BP_MAP[2]], //TAG
				];
				for first in 0..BP_MAP.len() {
					for second in 0..BP_MAP.len() {
						for third in 0..BP_MAP.len() {
							let key_next: u8 = (self.bp_map[BP_MAP[first]] << bp_shifts[1]) | (self.bp_map[BP_MAP[second]] << bp_shifts[0]) | self.bp_map[BP_MAP[third]];
							str_codon.insert(String::from_utf8(vec![BP_MAP_STR[first],BP_MAP_STR[second],BP_MAP_STR[third],]).unwrap(),key_next);
						}
					}
				}
			},
			0 =>  {
				self.start_codon = (self.bp_map[BP_MAP[2]] << bp_shifts[1]) + (self.bp_map[BP_MAP[3]] << bp_shifts[0]) + self.bp_map[BP_MAP[0]]; // ATG
				self.stop_codons = [
					(self.bp_map[BP_MAP[0]] << bp_shifts[1]) + (self.bp_map[BP_MAP[0]] << bp_shifts[0]) + self.bp_map[BP_MAP[3]], //TAA
					(self.bp_map[BP_MAP[0]] << bp_shifts[1]) + (self.bp_map[BP_MAP[2]] << bp_shifts[0]) + self.bp_map[BP_MAP[3]], //TGA
					(self.bp_map[BP_MAP[2]] << bp_shifts[1]) + (self.bp_map[BP_MAP[0]] << bp_shifts[0]) + self.bp_map[BP_MAP[3]], //TAG
				];
				for first in 0..BP_MAP.len() {
					for second in 0..BP_MAP.len() {
						for third in 0..BP_MAP.len() {
							let key_next: u8 = (self.bp_map[BP_MAP[third]] << bp_shifts[1]) | (self.bp_map[BP_MAP[second]] << bp_shifts[0]) | self.bp_map[BP_MAP[first]];
							str_codon.insert(String::from_utf8(vec![BP_MAP_STR[first],BP_MAP_STR[second],BP_MAP_STR[third],]).unwrap(),key_next);
						}
					}
				}
			},
			_ => unreachable!(),
		}
		self.gen_codon_map_raw(str_codon);

	}

	fn gen_codon_map_raw(&mut self, str_map: HashMap<String, u8>) {
		for item in str_map {
			if self.codon_map.contains_key(&item.1) { continue; }
			match item.0.as_str() {
				_ => { self.codon_map.insert(item.1, item.1); },
			}
		}
	}

	pub fn is_stop_codon(&self, codon_val: u8) -> bool {
		let ck: u8 = self.codon_map[&codon_val];
    	ck == self.stop_codons[0] || ck == self.stop_codons[1] || ck == self.stop_codons[2]
	}

	pub fn is_start_codon(&self, codon_val: u8) -> bool {
		let ck: u8 = self.codon_map[&codon_val];
    	ck == self.start_codon
	}
	pub fn reset_tmp_buf(&mut self) {
		self.tmp = 0;
		self.tmp_bits_used = 0;
	}

	pub fn reset_read_frame(&mut self) {
		self.read_frame = 0;
		self.read_frame_used = 0;
	}
	pub fn update_codon_buf(&mut self) {
		let tmp = self.read_frame as u16;
		self.codon_buff += tmp << self.codon_buff_used;
		self.codon_buff_used += self.read_frame_used as u16;
	}

}

pub struct Params {
    reader: BufReader<File>,
    writer: BufWriter<File>,
	code_fd: File,
	data_fd: Option<File>,
	ln_read_buf: String,
    compile_encode: EncodingMethod,
    ifile_ext: InputFileFormat,
	flags: u8, // b0: 0= lsb offset shifting, 1= simple fixed bit shift, b1=codons?, b2=in-place(1) or contiguous(0) data/code for codons. b3=codon inside = code(0) or data(1)
}
impl Params {
    pub fn from(ck: ParamsCK) -> Params {
		if ck.data_fp.is_none() {
        	Params {reader: ck.reader.unwrap(), writer: ck.writer.unwrap(),
				code_fd: OpenOptions::new().write(true).create(true).read(true).open(ck.code_fp.unwrap().as_os_str()).unwrap(),
				data_fd: None, ln_read_buf: String::new(), compile_encode: ck.compile_encode, ifile_ext: ck.ifile_ext, flags: ck.flags }
		}
		else {
			Params {reader: ck.reader.unwrap(), writer: ck.writer.unwrap(),
				code_fd: OpenOptions::new().write(true).create(true).read(true).open(ck.code_fp.unwrap().as_os_str()).unwrap(),
				data_fd: Some(OpenOptions::new().write(true).create(true).read(true).open(ck.data_fp.unwrap().as_os_str()).unwrap()),
				ln_read_buf: String::new(), compile_encode: ck.compile_encode, ifile_ext: ck.ifile_ext, flags: ck.flags }
		}
    }
    /**
     * Generates Base Pair value map based on encoding scheme.
	 *
     * Returns a 4-element byte array that stores the BP values in the following base pair order: [A, C, G, T]
     */
    pub fn bp_val_map(&self) -> [u8; 4] {
        match self.compile_encode {
            EncodingMethod::M1CS => [ method1_ascii_to_2bit_convert(b'A'), method1_ascii_to_2bit_convert(b'C'), method1_ascii_to_2bit_convert(b'G'), method1_ascii_to_2bit_convert(b'T')],
            EncodingMethod::M2A0C1 => [0, 1, 2, 3],
            EncodingMethod::M2A1C0 => [1, 0, 3, 2],
            EncodingMethod::M2A2C3 => [2, 3, 0, 1],
            EncodingMethod::M2A3C2 => [3, 2, 1, 0],
			EncodingMethod::M3A0C1 => [0, 1, 1, 0],
			EncodingMethod::M3A1C0 => [1, 0, 0, 1],
			EncodingMethod::M3A0C1XOR => [0, 1, 0, 1],
			EncodingMethod::M3A1C0XOR => [1, 0, 1, 0],
			EncodingMethod::INVALID => [255, 255, 255, 255],
        }
    }
	pub fn is_read_buf_empty(&self) -> bool { self.ln_read_buf.is_empty() }

}

const BUFFER_SIZE: usize = 1024;
const CODON_BIT_MASK: u8 = 0b00111111;
const CODON_BITS_LEN: u8 = 6;
const TEMP_BUFF_BITS: u8 = u8::BITS as u8;

fn is_stop_codon(codon_val: u8, stop_codons: [u8; 3]) -> bool {
    codon_val == stop_codons[0] || codon_val == stop_codons[1] || codon_val == stop_codons[2]
}

fn gen_codon_map(bp_tri: &str, bp_val: u8) -> (u8, u8) {
    let codon: u8 = match bp_tri {
        "AAA" | "AAG" => todo!(),                                 //Lys AAR
        "AAT" | "AAC" => todo!(),                                 //Asn AAY
        "GCT" | "GCC" | "GCA" | "GCG" => todo!(),                 //Ala GCN
        "CGT" | "CGC" | "CGA" | "CGG" | "AGA" | "AGG" => todo!(), //Arg CGN, AGR; or CGY, MGR
        "ATT" | "ATC" | "ATA" => todo!(),                         //Ile ATH
        "CTT" | "CTC" | "CTA" | "CTG" | "TTA" | "TTG" => todo!(), //Leu CTN, TTR; or CTY, YTR
        "GAT" | "GAC" => todo!(),                                 //Asp GAY
        "ATG" => todo!(), //Met (start codon) Alternates (if needed) TTG, GTG, CTG => NTG
        "TTT" | "TTC" => todo!(), // Phe TTY
        "TGT" | "TGC" => todo!(), //Cys TGY
        "CCT" | "CCC" | "CCA" | "CCG" => todo!(), //Pro CCN
        "CAA" | "CAG" => todo!(), //Gln CAR
        "TCT" | "TCC" | "TCA" | "TCG" | "AGT" | "AGC" => todo!(), //Ser TCN, AGY
        "GAA" | "GAG" => todo!(), //Glu GAR
        "ACT" | "ACC" | "ACA" | "ACG" => todo!(), //Thr ACN
        "TGG" => todo!(), //Trp
        "GGT" | "GGC" | "GGA" | "GGG" => todo!(), //Gly GCN
        "TAT" | "TAC" => todo!(), //Tyr TAY
        "CAT" | "CAC" => todo!(), //His CAY
        "GTT" | "GTC" | "GTA" | "GTG" => todo!(), //Val GTN
        "TAA" | "TGA" | "TAG" => todo!(), //Stop codon TRA TAR
        _ => 0b01111111,
    };
    (bp_val, codon)
}

fn gen_codon_map_raw(bp_tri: &str, bp_val: u8) -> (u8, u8) {
    let codon: u8 = match bp_tri {
        // "ATG" => 0xFF, //Met (start codon), 0xFF is a call opcode Alternates (if needed) TTG, GTG, CTG => NTG
        // "TAA" | "TGA" | "TAG" => 0xCB, //Stop codon, 0xCB = far return to calling procedure, 0xC3 is near return TRA TAR
        _ => bp_val,
    };
    (bp_val, codon)
}

fn generate_codon_map(bp_table: &mut [u8; 256]) -> HashMap<u8, u8> {
    const BP_MAP: [usize; 4] = ['A' as usize, 'C' as usize, 'G' as usize, 'T' as usize];
    const BP_MAP_STR: [u8; 4] = [b'A', b'C', b'G', b'T'];
    let mut codon_map: HashMap<String, u8> = HashMap::with_capacity(64);
    let mut h_map: HashMap<u8, u8> = HashMap::with_capacity(64);
    for first in 0..BP_MAP.len() {
        for second in 0..BP_MAP.len() {
            for third in 0..BP_MAP.len() {
                let key_next: u8 = (bp_table[BP_MAP[first]] << 4) | (bp_table[BP_MAP[second]] << 2) | bp_table[BP_MAP[third]];
                codon_map.insert(String::from_utf8(vec![BP_MAP_STR[first],BP_MAP_STR[second],BP_MAP_STR[third],]).unwrap(),key_next);
            }
        }
    }
    for item in codon_map {
        let k_v = gen_codon_map_raw(item.0.as_str(), item.1);
        h_map.insert(k_v.0, k_v.1);
    }
    return h_map;
}

fn codon_lsb_shift(mut params: &mut Params, mut table: &mut [u8; 256]) -> std::io::Result<()> {
	let mut codon_dt = CodonOnly::new(table, params.compile_encode, params.flags);
	codon_dt.dyn_init();
	let bit_mask: u8 = match codon_dt.bp_shift {
		1 => 0b01,
		2 => 0b11,
		_ => unreachable!(),
	};
	let mut code_bw: BufWriter<&File> = BufWriter::new(&params.code_fd);
	let codon_bits: u8 = codon_dt.bp_shift*3;
	let mut data_bw :BufWriter<&File> = if params.data_fd.is_some() {
		BufWriter::new(params.data_fd.as_ref().unwrap())
	}
	else { BufWriter::new(&params.code_fd) };
	buf_line_reader(&mut params.reader, &mut params.ln_read_buf)?;
	let mut other_buf: [u8; BUFFER_SIZE] = [0; BUFFER_SIZE];
	let mut codon_write_buf: [u8; BUFFER_SIZE] = [0; BUFFER_SIZE];
    let mut codon_write_buf_size: usize = 0;
	let mut other_buf_size: usize = 0;
	let mut total_bytes_written: usize = 0;

	while !params.ln_read_buf.is_empty() {
		for c in params.ln_read_buf.bytes() {
			codon_dt.read_frame += (codon_dt.bp_map[c as usize] << codon_dt.read_frame_used);
			codon_dt.read_frame_used += codon_dt.bp_shift;
			if codon_dt.read_frame_used < codon_bits { continue; }
			else if !codon_dt.inside_codon {
				codon_dt.inside_codon = codon_dt.is_start_codon(codon_dt.read_frame);
				if !codon_dt.inside_codon {
					codon_dt.tmp += (codon_dt.read_frame & bit_mask) << codon_dt.tmp_bits_used;
					codon_dt.tmp_bits_used += codon_dt.bp_shift;
					codon_dt.read_frame >>= codon_dt.bp_shift;
					codon_dt.read_frame_used -= codon_dt.bp_shift;
					if codon_dt.tmp_bits_used == TEMP_BUFF_BITS {
						other_buf[other_buf_size] = codon_dt.tmp;
						codon_dt.reset_tmp_buf();
						other_buf_size += 1;
						if other_buf_size == BUFFER_SIZE {
							total_bytes_written += buf_writer(&mut data_bw, &other_buf, other_buf_size)?;
							other_buf_size = 0;
						}
					}
				}
				else {
					codon_dt.codon_buff += (codon_dt.read_frame as u16) << codon_dt.codon_buff_used;
					codon_dt.codon_buff_used += (codon_dt.read_frame_used as u16);
					codon_dt.reset_read_frame();
					if codon_dt.codon_buff_used >= 8 {
						codon_write_buf[codon_write_buf_size] = codon_dt.codon_buff as u8;
						codon_dt.codon_buff >>= 8;
						codon_dt.codon_buff_used -= 8;
						codon_write_buf_size += 1;
						if codon_write_buf_size == BUFFER_SIZE {
							total_bytes_written += buf_writer(&mut code_bw, &codon_write_buf, codon_write_buf_size)?;
							codon_write_buf_size = 0;
						}
					}
					//codon_dt.update_codon_buf();
				}
			}
			else {
				codon_dt.inside_codon = !codon_dt.is_stop_codon(codon_dt.read_frame);
				codon_dt.codon_buff += (codon_dt.read_frame as u16) << codon_dt.codon_buff_used;
				codon_dt.codon_buff_used += (codon_dt.read_frame_used as u16);
				codon_dt.reset_read_frame();
				//codon_dt.update_codon_buf();

				if codon_dt.codon_buff_used >= 8 {
					codon_write_buf[codon_write_buf_size] = codon_dt.codon_buff as u8;
					codon_dt.codon_buff >>= 8;
					codon_dt.codon_buff_used -= 8;
					codon_write_buf_size += 1;
					if codon_write_buf_size == BUFFER_SIZE {
						total_bytes_written += buf_writer(&mut code_bw, &codon_write_buf, codon_write_buf_size)?;
						codon_write_buf_size = 0;
					}
				}
			}
		}
		params.ln_read_buf.clear();
		buf_line_reader(&mut params.reader, &mut params.ln_read_buf)?;
	}
	//TODO: write buffers and temp holding variables to appropriate files. Write only full codons to code file.
	if codon_dt.read_frame_used < codon_bits || !codon_dt.inside_codon {
		let tmp_remain = TEMP_BUFF_BITS - codon_dt.read_frame;
		let overflow = tmp_remain < codon_dt.read_frame_used;

		codon_dt.tmp += codon_dt.read_frame << codon_dt.tmp_bits_used;
		codon_dt.tmp_bits_used += tmp_remain;
		other_buf[other_buf_size] = codon_dt.tmp;
		codon_dt.reset_tmp_buf();
		other_buf_size += 1;
		if other_buf_size == BUFFER_SIZE {
			total_bytes_written += buf_writer(&mut data_bw, &other_buf, other_buf_size)?;
			other_buf_size = 0;
		}
		if overflow {
			codon_dt.read_frame >>= tmp_remain;
			codon_dt.read_frame_used -= tmp_remain;
			codon_dt.tmp += codon_dt.read_frame;
			codon_dt.tmp_bits_used += codon_dt.read_frame_used;
		}
		codon_dt.reset_read_frame();
	}
	else {
		codon_dt.codon_buff += (codon_dt.read_frame as u16) << codon_dt.codon_buff_used;
		codon_dt.codon_buff_used += (codon_dt.read_frame_used as u16);
		codon_dt.reset_read_frame();
	}

	//TODO: check if inside codon or not, handle appropriately.
	//TODO: check other buffer and tmp var, write to file


	if codon_dt.tmp_bits_used > 0 {
		other_buf[other_buf_size] = codon_dt.tmp;
		codon_dt.reset_tmp_buf();
		other_buf_size += 1;
	}
	if other_buf_size > 0 {
		total_bytes_written += buf_writer(&mut data_bw, &other_buf, other_buf_size)?;
		other_buf_size = 0;
	}
	while codon_dt.codon_buff_used >= 8 {
		codon_write_buf[codon_write_buf_size] = codon_dt.codon_buff as u8;
		codon_dt.codon_buff >>= 8;
		codon_dt.codon_buff_used -= 8;
		codon_write_buf_size += 1;
		if codon_write_buf_size == BUFFER_SIZE {
			total_bytes_written += buf_writer(&mut code_bw, &codon_write_buf, codon_write_buf_size)?;
			codon_write_buf_size = 0;
		}
	}
	if codon_dt.codon_buff_used > 0 {
		codon_write_buf[codon_write_buf_size] = codon_dt.codon_buff as u8;
		codon_dt.codon_buff >>= codon_dt.codon_buff_used;
		codon_dt.codon_buff_used = 0;
		codon_write_buf_size += 1;
	}
	if codon_write_buf_size > 0 {
		total_bytes_written += buf_writer(&mut code_bw, &codon_write_buf, codon_write_buf_size)?;
		codon_write_buf_size = 0;
	}
	Ok(())
}


fn codon_msb_simple_left_shift(mut params: &mut Params, mut table: &mut [u8; 256]) -> std::io::Result<()> {
    const BP_MAP: [usize; 8] = ['A' as usize, 'a' as usize, 'C' as usize, 'c' as usize, 'G' as usize, 'g' as usize, 'T' as usize, 't' as usize];
	let mut codon_dt = CodonOnly::new(table, params.compile_encode, params.flags);
	codon_dt.dyn_init();
	let mut code_bw: BufWriter<&File> = BufWriter::new(&params.code_fd);
	let mut data_bw = if params.data_fd.is_some() {
		BufWriter::new((params.data_fd.as_ref().unwrap()))
	}
	else { BufWriter::new(&params.code_fd) };
    buf_line_reader(&mut params.reader, &mut params.ln_read_buf)?;
    let mut buff_2bit: [u8; BUFFER_SIZE] = [0; BUFFER_SIZE];
	let mut codon_write_buf: [u8; BUFFER_SIZE] = [0; BUFFER_SIZE];
    let mut codon_write_buf_size: usize = 0;
	let mut buff_2bit_temp: u8 = 0;
    let mut buff_2bit_size: usize = 0;
    let mut read_frame : u8 = 0; // 6bits
	let mut read_frame_used : u8 = 0;
	let mut codon_buff: u64 = 0;
	let mut codon_buff_remain: u8 = 64;
	let mut temp_bits_used: u8 = 0;
    let mut total_bytes_written: u64 = 0;
    while !params.ln_read_buf.is_empty() {
        for c in params.ln_read_buf.bytes() {
			read_frame <<= 2;
			read_frame += codon_dt.bp_map[c as usize];
			read_frame_used += 2;
			if read_frame_used < CODON_BITS_LEN { continue; }
            if read_frame_used == TEMP_BUFF_BITS {
                read_frame &= CODON_BIT_MASK;
                read_frame_used -= 2;
            }
			else if !codon_dt.inside_codon {
				codon_dt.inside_codon = codon_dt.is_start_codon(read_frame);
				if !codon_dt.inside_codon {
					buff_2bit_temp <<= 2;
					buff_2bit_temp += read_frame >> 4;
					temp_bits_used += 2;
					 if temp_bits_used == TEMP_BUFF_BITS {
						buff_2bit[buff_2bit_size] = buff_2bit_temp;
						temp_bits_used = 0;
						buff_2bit_temp = 0;
						buff_2bit_size += 1;
						if buff_2bit_size == BUFFER_SIZE {
							total_bytes_written += buf_writer(&mut params.writer, &buff_2bit, buff_2bit_size)? as u64;
							buff_2bit_size = 0;
						}
					}
				}
				else {
					codon_buff = read_frame as u64;
					read_frame = 0;
					codon_buff_remain -= 8;
					read_frame_used = 0;
					//total_bytes_written += buf_writer_align(&mut code_bw, &mut buff_2bit, buff_2bit_size, buff_2bit_temp, temp_bits_used, total_bytes_written, b'\0')?;
				}
			}
			else {
                // if read_frame_used == TEMP_BUFF_MAX_BITS {
				// 	read_frame &= CODON_BIT_MASK;
				// 	read_frame_used -= 2;
				// }
                codon_dt.inside_codon = codon_dt.is_stop_codon(read_frame);
                if codon_dt.inside_codon {
                    if codon_buff_remain >= read_frame_used {
                        codon_buff <<= read_frame_used;
                        codon_buff += read_frame as u64;
                        codon_buff_remain -= read_frame_used;
                        read_frame_used = 0;
                        read_frame = 0;
                	}
                    else {
                        //total_bytes_written += buf_64_mod_6_align(&mut code_bw, &mut codon_write_buf, &mut codon_write_buf_size, read_frame, read_frame_used, &mut codon_buff, &mut codon_buff_remain)?;
                    }
                }
                else {
                    // println!("{:>width$b} ->", codon_buff, width=(64-codon_buff_remain) as usize);
                    // if codon_buff_remain > 0 {
                    //     codon_buff <<= codon_buff_remain;
                    //     // println!("{:b}", codon_buff);

                	// }
                    if codon_write_buf_size+8 > BUFFER_SIZE  {
                        //total_bytes_written += buf_writer(&mut params.writer, &mut codon_write_buf, codon_write_buf_size)? as u64;
                        codon_write_buf_size = 0;
                    }
                    for off in 7..=0 {
                        codon_write_buf[codon_write_buf_size+off] = codon_buff as u8;
                        codon_buff >>= 8;
                    }
                    codon_buff = 0;
                    codon_buff_remain = 64;
                    codon_write_buf_size += 8;
                    let align = if (codon_write_buf_size+1) % 16 != 0 { (codon_write_buf_size+1) % 16 } else { 15 };
                    if codon_write_buf_size+align > BUFFER_SIZE {
                        //total_bytes_written += buf_writer(&mut params.writer, &mut codon_write_buf, codon_write_buf_size)? as u64;
                        codon_write_buf_size = 0;
                    }

                    codon_write_buf[codon_write_buf_size] = read_frame;
                    codon_write_buf_size += 1;
                    read_frame_used = 0;
                    read_frame = 0;
                    // for a in 0..align {
                    //     codon_write_buf[codon_write_buf_size+a] = NOOP_CODE;
                    // }
                    // codon_write_buf_size += align;
                    //total_bytes_written += buf_writer(&mut params.writer, &mut codon_write_buf, codon_write_buf_size)? as u64;
                    codon_write_buf_size = 0;
                }
            }
        }
        params.ln_read_buf.clear();
        buf_line_reader(&mut params.reader, &mut params.ln_read_buf)?;
    }
    // if read_frame_used < CODON_BITS_LEN { todo!(); }
    if read_frame_used == TEMP_BUFF_BITS {
        read_frame &= CODON_BIT_MASK;
        read_frame_used -= 2;
    }
    if !codon_dt.inside_codon {

    //   (buff_2bit_temp, temp_bits_used, read_frame,read_frame_used) = eof_8_mod_6(buff_2bit_temp, temp_bits_used, read_frame, read_frame_used);

        if temp_bits_used == TEMP_BUFF_BITS {
            buff_2bit[buff_2bit_size] = buff_2bit_temp;
            temp_bits_used = 0;
            buff_2bit_temp = 0;
            buff_2bit_size += 1;
            if read_frame_used > 0 {
                buff_2bit_temp = read_frame;
                temp_bits_used = read_frame_used;
                read_frame = 0;
                read_frame_used = 0;
            }
            if buff_2bit_size == BUFFER_SIZE {
                //total_bytes_written += buf_writer(&mut params.writer, &buff_2bit, buff_2bit_size)? as u64;
                buff_2bit_size = 0;
            }
        }
        //total_bytes_written += buf_writer_align(&mut params.writer, &mut buff_2bit, buff_2bit_size, buff_2bit_temp, temp_bits_used, total_bytes_written, b'\0')? as u64;
    }
    else {
        let align_fill: u8 = if  read_frame_used == 6 && codon_dt.is_stop_codon(read_frame) { read_frame } else { CODE_PADDING };

        if align_fill != CODE_PADDING {
            if codon_buff_remain >= read_frame_used {
                codon_buff <<= read_frame_used;
                codon_buff += read_frame as u64;
                codon_buff_remain -= read_frame_used;
                read_frame_used = 0;
                read_frame = 0;
            }
            else {
                //total_bytes_written += buf_64_mod_6_align(&mut params.writer, &mut codon_write_buf, &mut codon_write_buf_size, read_frame, read_frame_used, &mut codon_buff, &mut codon_buff_remain)?;

            }
            // if codon_buff_remain > 0 {
            //     codon_buff <<= codon_buff_remain;
            // }
            if codon_write_buf_size+8 > BUFFER_SIZE  {
                //total_bytes_written += buf_writer(&mut params.writer, &mut codon_write_buf, codon_write_buf_size)? as u64;
                codon_write_buf_size = 0;
            }
            for off in 7..=0 {
                //codon_write_buf[codon_write_buf_size+off] = codon_buff as u8;
                codon_buff >>= 8;
            }
            codon_buff = 0;
            codon_buff_remain = 64;
            codon_write_buf_size += 8;
        }
        else {
            // if codon_buff_remain > 0 {
            //     codon_buff <<= codon_buff_remain;
            // }
            if codon_write_buf_size+8 > BUFFER_SIZE  {
                //total_bytes_written += buf_writer(&mut params.writer, &mut codon_write_buf, codon_write_buf_size)? as u64;
                codon_write_buf_size = 0;
            }
            for off in 7..=0 {
                codon_write_buf[codon_write_buf_size+off] = codon_buff as u8;
                codon_buff >>= 8;
            }
            codon_buff = 0;
            codon_buff_remain = 64;
            codon_write_buf_size += 8;
            // total_bytes_written += buf_64_mod_6_align(&mut params.writer, &mut codon_write_buf, codon_write_buf_size, read_frame, read_frame_used, &mut codon_buff, &mut codon_buff_remain)?;
            // codon_write_buf_size = 0;
        }

        let align = if (codon_write_buf_size+1) % 16 != 0 { (codon_write_buf_size+1) % 16 } else { 15 } ;
        if codon_write_buf_size+align > BUFFER_SIZE {
            //total_bytes_written += buf_writer(&mut params.writer, &mut codon_write_buf, codon_write_buf_size)? as u64;
            codon_write_buf_size = 0;
        }

        codon_write_buf[codon_write_buf_size] = align_fill;
        codon_write_buf_size += 1;
        read_frame_used = 0;
        read_frame = 0;
        for a in 0..align {
            codon_write_buf[codon_write_buf_size+a] = CODE_PADDING;
        }
        codon_write_buf_size += align;
        //total_bytes_written += buf_writer(&mut params.writer, &mut codon_write_buf, codon_write_buf_size)? as u64;
        codon_write_buf_size = 0;
    }

    // if temp_bits_used < 8 {
    //     buff_2bit[buff_2bit_size] = buff_2bit_temp;
    //     temp_bits_used = 0;
    //     buff_2bit_temp = 0;
    //     buff_2bit_size += 1;

    // }
    // if buff_2bit_size > 0 {
    //     total_bytes_written += buf_writer(&mut params.writer, &buff_2bit, buff_2bit_size)? as u64;
    //     buff_2bit_size = 0;
    // }
    // let align = total_bytes_written % (size_of::<usize>() as u64); //align to system default word bit size. Expect 64 most of the time
    // if align != 0 {
    //     println!("padding: {} bytes", align);
    //     for i in 0..align as usize {
    //         buff_2bit[i] = NOOP_CODE;
    //         buff_2bit_size += 1;
    //     }
    //     total_bytes_written += buf_writer(&mut params.writer, &buff_2bit, buff_2bit_size)? as u64;
    //     buff_2bit_size = 0;
    // }
    // println!("");
    code_bw.flush()?;

    Ok(())
}

fn ms_bit_first_main(mut params: &mut Params, &mut table: &mut [u8; 256]) -> std::io::Result<()> {
	let bp_shift: u8 = match params.compile_encode {
		EncodingMethod::M3A0C1 | EncodingMethod::M3A1C0 | EncodingMethod::M3A0C1XOR | EncodingMethod::M3A1C0XOR => 1,
		_ => 2,
	};
	let mut code_bw: BufWriter<&File> = BufWriter::new(&params.code_fd);
    buf_line_reader(&mut params.reader, &mut params.ln_read_buf).unwrap_or_else(|e| { panic!("Failed to read line in open input file: {}", e); });
    let mut buff_2bit: [u8; BUFFER_SIZE] = [0; BUFFER_SIZE];
	let mut buff_2bit_temp: u8 = 0;
    let mut buff_2bit_size: usize = 0;
    let mut temp_bits_used: u8 = 0;
    let mut total_bytes_written: u64 = 0;
    while !params.is_read_buf_empty() {
        for c in params.ln_read_buf.bytes() {
			buff_2bit_temp <<= bp_shift;
			buff_2bit_temp += table[c as usize];
			temp_bits_used += bp_shift;
			if temp_bits_used == TEMP_BUFF_BITS {
				buff_2bit[buff_2bit_size] = buff_2bit_temp;
				temp_bits_used = 0;
				buff_2bit_temp = 0;
				buff_2bit_size += 1;
				if buff_2bit_size == BUFFER_SIZE {
					total_bytes_written += buf_writer(&mut code_bw, &buff_2bit, buff_2bit_size)? as u64;
                    buff_2bit_size = 0;
				}
			}
        }
        params.ln_read_buf.clear();
        buf_line_reader(&mut params.reader, &mut params.ln_read_buf)?;
    }
    if buff_2bit_size == BUFFER_SIZE {
        total_bytes_written += buf_writer(&mut code_bw, &buff_2bit, buff_2bit_size)? as u64;
        buff_2bit_size = 0;
    }
    if temp_bits_used > 0 {
        buff_2bit[buff_2bit_size] = buff_2bit_temp;
        temp_bits_used = 0;
        buff_2bit_temp = 0;
        buff_2bit_size += 1;

    }
    if buff_2bit_size > 0 {
        total_bytes_written += buf_writer(&mut code_bw, &buff_2bit, buff_2bit_size)? as u64;
        buff_2bit_size = 0;
    }
    let align = total_bytes_written % 4; //align to 32-bit words
    if align != 0 {
        for i in 0..align as usize {
            buff_2bit[i] = CODE_PADDING;
            buff_2bit_size += 1;
        }
        total_bytes_written += buf_writer(&mut code_bw, &buff_2bit, buff_2bit_size)? as u64;
        buff_2bit_size = 0;
    }
    code_bw.flush()?;

    Ok(())
}

fn ls_bit_first_main(mut params: &mut Params, &mut table: &mut [u8; 256]) -> std::io::Result<()> {
	let bp_shift: u8 = match params.compile_encode {
		EncodingMethod::M3A0C1 | EncodingMethod::M3A1C0 | EncodingMethod::M3A0C1XOR | EncodingMethod::M3A1C0XOR => 1,
		_ => 2,
	};
	let mut code_bw: BufWriter<&File> = BufWriter::new(&params.code_fd);
    buf_line_reader(&mut params.reader, &mut params.ln_read_buf).unwrap_or_else(|e| { panic!("Failed to read line in open input file: {}", e); });
    let mut buff_2bit: [u8; BUFFER_SIZE] = [0; BUFFER_SIZE];
	let mut buff_2bit_temp: u8 = 0;
    let mut buff_2bit_size: usize = 0;
    let mut temp_bits_used: u8 = 0;
    let mut total_bytes_written: u64 = 0;
    while !params.is_read_buf_empty() {
        for c in params.ln_read_buf.bytes() {
			buff_2bit_temp += (table[c as usize]) << temp_bits_used;
			temp_bits_used += bp_shift;
			if temp_bits_used == TEMP_BUFF_BITS {
				buff_2bit[buff_2bit_size] = buff_2bit_temp;
				temp_bits_used = 0;
				buff_2bit_temp = 0;
				buff_2bit_size += 1;
				if buff_2bit_size == BUFFER_SIZE {
					total_bytes_written += buf_writer(&mut code_bw, &buff_2bit, buff_2bit_size)? as u64;
                    buff_2bit_size = 0;
				}
			}
        }
        params.ln_read_buf.clear();
        buf_line_reader(&mut params.reader, &mut params.ln_read_buf).unwrap_or_else(|e| { panic!("Failed to read line in open input file: {}", e); });
    }
    if buff_2bit_size == BUFFER_SIZE {
        total_bytes_written += buf_writer(&mut code_bw, &buff_2bit, buff_2bit_size)? as u64;
        buff_2bit_size = 0;
    }
    if temp_bits_used > 0 {
        buff_2bit[buff_2bit_size] = buff_2bit_temp;
        temp_bits_used = 0;
        buff_2bit_temp = 0;
        buff_2bit_size += 1;

    }
    if buff_2bit_size > 0 {
        total_bytes_written += buf_writer(&mut code_bw, &buff_2bit, buff_2bit_size)? as u64;
        buff_2bit_size = 0;
    }
    let align = total_bytes_written % 4; //align to 32-bit word size
    if align != 0 {
        for i in 0..align as usize {
            buff_2bit[i] = CODE_PADDING;
            buff_2bit_size += 1;
        }
        total_bytes_written += buf_writer(&mut code_bw, &buff_2bit, buff_2bit_size)? as u64;
        buff_2bit_size = 0;
    }
    code_bw.flush()?;
    Ok(())
}



#[cfg(target_os="windows")]
fn calc_alignment(len: u32, align_val: u32) -> u32 {
	(match len % align_val {
		0 => len/align_val,
		_ => (len/align_val) + 1,
	})*align_val
}

#[cfg(target_os="windows")]
fn gen_sect_ent_data(sec_sizes: Vec<u32>, num_sections: u16) {
	let mut sec_data: Vec<u32> = Vec::with_capacity((num_sections as usize+1)*3);
	sec_data.push(win::get_header_size(num_sections as usize) as u32); //header data size in bytes
	sec_data.push(calc_alignment(sec_data[0], FILE_ALIGNMENT)); //Header size aligned to file alignment
	sec_data.push(0); //va address section offset

	for i in 1..(sec_sizes.len()+1) {
		sec_data.push(sec_sizes[i-1]);
		sec_data.push(calc_alignment(sec_data[i*3], FILE_ALIGNMENT));
		sec_data.push(sec_data[(i-1)+2] + calc_alignment(sec_data[i*3], SECTION_ALIGNMENT));
	}
}

#[cfg(target_os="windows")]
fn generate_win_executable_header(sec_sizes: Vec<u32>, num_sections: u16) -> HeaderPE {
	//let code_len = sec_sizes[0];
	if num_sections == 1 {
		return generate_win_executable_code_only(sec_sizes[0]);
	}
	let data_dir_list: [ImageDataDirectory; 16] = [ImageDataDirectory::default(); 16];

	let header_data_len: u32 = win::get_header_size(num_sections as usize) as u32;

	let header_file_align: u32 = calc_alignment(header_data_len, FILE_ALIGNMENT);
	let va_code_offset: u32 = calc_alignment(header_data_len,SECTION_ALIGNMENT);
	let code_size_align: u32 = calc_alignment(sec_sizes[0], FILE_ALIGNMENT);
	let data_size_align: u32 = calc_alignment(sec_sizes[1], FILE_ALIGNMENT);
	let va_data_offset: u32 = va_code_offset + calc_alignment(sec_sizes[0], SECTION_ALIGNMENT);
	let code_file_align: u32 = header_file_align + calc_alignment(sec_sizes[0], FILE_ALIGNMENT);
	let img_sz: u32 = va_data_offset + calc_alignment(sec_sizes[1],SECTION_ALIGNMENT);

	let opt_file_header: OptionalImageFileHeader32 = OptionalImageFileHeader32::with_code_and_data(code_size_align, data_size_align, va_code_offset, va_code_offset, va_data_offset, img_sz, data_dir_list);
	let sect_header_list: Vec<ImageSectionHeaderEnt> = vec![
		ImageSectionHeaderEnt::with_perms(*b".text\0\0\0", sec_sizes[0], va_code_offset, code_size_align, header_file_align, 0x60000020),
		ImageSectionHeaderEnt::with_perms(*b".data\0\0\0", sec_sizes[1], va_data_offset, data_size_align, code_file_align, 0xC0000040)
	];
	HeaderPE::with_defaults(num_sections, opt_file_header, sect_header_list)


}

#[cfg(target_os="windows")]
fn generate_win_executable_code_only(code_len: u32) -> HeaderPE {
	let data_dir_list: [ImageDataDirectory; 16] = [ImageDataDirectory::default(); 16];
	let header_data_len: u32 = win::get_header_size(1) as u32;
	let header_file_align: u32 = calc_alignment(header_data_len, FILE_ALIGNMENT);
	let sec_file_align: u32 = calc_alignment(code_len, FILE_ALIGNMENT);
	let header_va_offset: u32 = calc_alignment(header_data_len,SECTION_ALIGNMENT);
	let va_data_addr_offset: u32 = header_va_offset + calc_alignment(code_len, SECTION_ALIGNMENT);
	let opt_file_header: OptionalImageFileHeader32 = OptionalImageFileHeader32::with_code_only(sec_file_align, header_va_offset, header_va_offset, va_data_addr_offset, data_dir_list);
	let nt_header_32: ImageNTHeaders32 = ImageNTHeaders32::with_defaults(1, opt_file_header);
	HeaderPE {
		init_header: InitFileHeader::default(), nt_header: nt_header_32,
		section_header: vec![
			ImageSectionHeaderEnt::with_perms(*b".text\0\0\0", code_len, header_va_offset, sec_file_align, header_file_align, 0x60000020)
			]
	}

}

#[cfg(target_os="linux")]
fn generate_elf_executable_code_only(code_len: u32) -> FileHeader {
	let phead_size: usize = linux::get_header_size(2, 0);
	let p_head_list: Vec<ELF32ProgHeaderEnt> = vec![ELF32ProgHeaderEnt::new(phead_size as u32), ELF32ProgHeaderEnt::with_params_flags(1, 0x1000, code_len, 0x05, 0x1000)];
	//let sec_head_list: Vec<ELF32SectHeaderEnt> = vec![ELF32SectHeaderEnt::default(), ELF32SectHeaderEnt::with_params_flags(1, 1, 0x06, 0x1000, code_len, 0, 0, 16, 0)];
	FileHeader::with_defaults(ELF32Header::new(0x1000,2), p_head_list)


}

fn generate_executable(mut params: &mut Params) -> std::io::Result<()> {
	const ALIGN_FILE: usize = FILE_ALIGNMENT as usize;
	let mut f_sizes: Vec<u32> = Vec::new();
	let code_size: u64 = params.code_fd.metadata()?.len();
	let mut temp_reader: BufReader<&File> = BufReader::new(&params.code_fd);
	temp_reader.rewind()?;
	f_sizes.push(code_size as u32);
	let f_data_valid = params.flags >= 2 && params.data_fd.is_some();
	if f_data_valid {
		f_sizes.push(params.data_fd.as_ref().unwrap().metadata()?.len() as u32);
	}

	#[cfg(target_os="windows")]
	let header = generate_win_executable_header(f_sizes.clone(), f_sizes.len() as u16);
	#[cfg(target_os="linux")]
	let header = generate_elf_executable_code_only(f_sizes[0]);
	let header_buff = header.serialize();
	let mut bytes_written: usize = buf_writer(&mut params.writer, &header_buff.as_slice(), header_buff.len())?;
	let buff_align: usize = ALIGN_FILE - (bytes_written % ALIGN_FILE);
	let mut write_buff: [u8; ALIGN_SECTION] = [HEADER_PADDING_BYTE; ALIGN_SECTION];
	if buff_align < ALIGN_FILE {
		bytes_written += buf_writer(&mut params.writer, &write_buff, buff_align)?;
	}
	params.writer.flush()?;
	for _full_reads in 0..f_sizes[0]/SECTION_ALIGNMENT {
		temp_reader.read_exact(&mut write_buff).unwrap_or_else(|e| { panic!("code.bin file reader error: {}", e); });
		bytes_written += buf_writer(&mut params.writer, &write_buff,ALIGN_SECTION)?;
	}
	let mut end_vec: Vec<u8> = Vec::with_capacity(ALIGN_SECTION);
	let code_data_end = temp_reader.read_to_end(&mut end_vec)?;
	bytes_written += buf_writer(&mut params.writer, &end_vec.as_slice(), end_vec.len())?;
	if code_data_end % ALIGN_SECTION > 0 {
		write_buff = [CODE_PADDING; ALIGN_SECTION];
		bytes_written += buf_writer(&mut params.writer, &write_buff, ALIGN_SECTION -(code_data_end % ALIGN_SECTION))?;
	}
	params.writer.flush()?;
	if f_data_valid {
		let mut data_reader: BufReader<&File> = BufReader::new(params.data_fd.as_mut().unwrap());
		data_reader.rewind()?;
		write_buff = [CODE_PADDING; ALIGN_SECTION];
		for _full_buf in 0..f_sizes[1]/SECTION_ALIGNMENT {
			data_reader.read_exact(&mut write_buff).unwrap_or_else(|e| { panic!("data.bin file reader error: {}", e); });
			bytes_written += buf_writer(&mut params.writer, &write_buff,ALIGN_SECTION)?;
		}
		let mut end_vec: Vec<u8> = Vec::with_capacity(ALIGN_SECTION);
		let data_end = data_reader.read_to_end(&mut end_vec)?;
		bytes_written += buf_writer(&mut params.writer, &end_vec.as_slice(), end_vec.len())?;
		if data_end % ALIGN_SECTION > 0 {
			write_buff = [CODE_PADDING; ALIGN_SECTION];
			bytes_written += buf_writer(&mut params.writer, &write_buff, ALIGN_SECTION -(data_end % ALIGN_SECTION))?;
		}
		params.writer.flush()?;
	}


	Ok(())
}
fn remove_temp_files(temp_dir_p: &str, temp_files: Vec<&str>) -> std::io::Result<()> {
	for item in temp_files {
		let mut fp = PathBuf::from(temp_dir_p);
		fp.push(item);
		if fp.exists() && fp.is_file() { std::fs::remove_file(fp.as_path())?; }
	}
	Ok(())
}
fn main() -> std::io::Result<()> {
    const BP_MAP: [usize; 8] = ['A' as usize, 'a' as usize, 'C' as usize, 'c' as usize, 'G' as usize, 'g' as usize, 'T' as usize, 't' as usize,];
    let mut table: [u8; 256] = [255; 256];
    let main_args: Vec<String> = env::args().collect();
    let (check, temp_path) = parse_main_args(&main_args)?;
    if !check.check_files() { return Ok(()); }
	if !check.is_req_valid() { Error::new(std::io::ErrorKind::InvalidData, "Parameters are invalid!"); }

	let mut params: Params = Params::from(check);
    let bp_map_val: [u8; 4] = params.bp_val_map();
    table[BP_MAP[0]] = bp_map_val[0];
    table[BP_MAP[1]] = bp_map_val[0];
    table[BP_MAP[2]] = bp_map_val[1];
    table[BP_MAP[3]] = bp_map_val[1];
    table[BP_MAP[4]] = bp_map_val[2];
    table[BP_MAP[5]] = bp_map_val[2];
    table[BP_MAP[6]] = bp_map_val[3];
    table[BP_MAP[7]] = bp_map_val[3];
	match params.flags {
		0b0 => ls_bit_first_main(&mut params, &mut table)?,
		0b1 => ms_bit_first_main(&mut params, &mut table)?,
		0b10 => codon_lsb_shift(&mut params, &mut table)?,
		_ => unreachable!(),
	}
	params.code_fd.sync_all()?;
	if params.flags >= 2 && params.data_fd.is_some() { //codons
		params.data_fd.as_mut().unwrap().sync_all()?;
	}
	generate_executable(&mut params)?;
	return remove_temp_files(temp_path.as_str(),vec!["code.bin", "data.bin"]);


}

#[inline(always)]
fn method1_ascii_to_2bit_convert(c: u8) -> u8 { (c & 0x06) >> 1 }

fn buf_writer<W: Write>(writer: &mut W, buf: &[u8], buff_len: usize) -> std::io::Result<usize> {
	let mut written_bytes: usize = 0;
	while written_bytes < buff_len {
		let n: usize = writer.write(&buf[written_bytes..buff_len])?;
		written_bytes += n;
	}
	Ok(written_bytes)
}

fn buf_line_reader<R: BufRead>(reader: &mut R, buf: &mut String) -> std::io::Result<usize> {
	let mut len: usize = reader.read_line(buf)?;
	if buf.chars().next().unwrap_or_default() == '>' {
		buf.clear();
		len = reader.read_line(buf)?;
	}

	if !buf.is_empty() && buf.chars().last().unwrap() == '\n' {
		buf.pop(); // remove \n character
		len -= 1;
		if len > 0 && buf.chars().last().unwrap() == '\r' {
			buf.pop();
			len -= 1;
		}
	}
	Ok(len)
}
