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
#[repr(u8)]
pub enum CodonState {
	OutRead = 0b00,
	ToInRd = 0b01,
	ToOutRd = 0b10,
	InRead = 0b11
}
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
		self.codon_buff += (self.read_frame as u16) << self.codon_buff_used;
		self.codon_buff_used += self.read_frame_used as u16;
		self.reset_read_frame();
	}

	//pub fn flush_tmp(&mut self, )


}

pub struct Params {
    reader: BufReader<File>,
    writer: BufWriter<File>,
	code_fd: File,
	data_fd: Option<File>,
	ln_read_buf: String,
    compile_encode: EncodingMethod,
    #[allow(dead_code)]
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
            EncodingMethod::M1CS => [ self.method1_ascii_to_2bit_convert(b'A'), self.method1_ascii_to_2bit_convert(b'C'), self.method1_ascii_to_2bit_convert(b'G'), self.method1_ascii_to_2bit_convert(b'T')],
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
	pub fn read_buf_not_empty(&self) -> bool { !self.ln_read_buf.is_empty() }
	//pub fn clear_read_buf(&mut self) { self.ln_read_buf.clear(); }
	pub fn flag_check(&self, mask : BitFlags) -> bool { self.flags & mask as u8 > 0 }

	#[inline(always)]
	fn method1_ascii_to_2bit_convert(&self, c: u8) -> u8 { (c & 0x06) >> 1 }


}

const BUFFER_SIZE: usize = 1024;
const TEMP_BUFF_BITS: u8 = u8::BITS as u8;
const PADDING_BUFF: [u8; 4] = [CODE_PADDING; 4];

#[allow(dead_code)]
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

#[allow(dead_code)]
fn gen_codon_map_raw(bp_tri: &str, bp_val: u8) -> (u8, u8) {
    let codon: u8 = match bp_tri {
        // "ATG" => 0xFF, //Met (start codon), 0xFF is a call opcode Alternates (if needed) TTG, GTG, CTG => NTG
        // "TAA" | "TGA" | "TAG" => 0xCB, //Stop codon, 0xCB = far return to calling procedure, 0xC3 is near return TRA TAR
        _ => bp_val,
    };
    (bp_val, codon)
}

#[allow(dead_code)]
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



fn codon_lsb_shift(params: &mut Params, table: &[u8; 256]) -> std::io::Result<()> {
	let mut codon_dt = CodonOnly::new(table, params.compile_encode, params.flags);
	codon_dt.dyn_init();
	let bit_mask: u8 = match codon_dt.bp_shift {
		1 => 0b01,
		2 => 0b11,
		_ => unreachable!(),
	};
	let mut code_bw: BufWriter<&File> = BufWriter::new(&params.code_fd);
	let add_padding: bool = params.flag_check(BitFlags::CodeSectionPadding);
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
	let mut _total_bytes_written: usize = 0;
	if add_padding && params.read_buf_not_empty() {
		code_bw.write_all(&PADDING_BUFF)?;
		_total_bytes_written = 4;
	}

	while params.read_buf_not_empty() {
		for c in params.ln_read_buf.bytes() {
			codon_dt.read_frame += codon_dt.bp_map[c as usize] << codon_dt.read_frame_used;
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
							_total_bytes_written += buf_writer(&mut data_bw, &other_buf, other_buf_size)?;
							other_buf_size = 0;
						}
					}
				}
				else {
					codon_dt.update_codon_buf();
					if codon_dt.codon_buff_used >= 8 {
						codon_write_buf[codon_write_buf_size] = codon_dt.codon_buff as u8;
						codon_dt.codon_buff >>= 8;
						codon_dt.codon_buff_used -= 8;
						codon_write_buf_size += 1;
						if codon_write_buf_size == BUFFER_SIZE {
							_total_bytes_written += buf_writer(&mut code_bw, &codon_write_buf, codon_write_buf_size)?;
							codon_write_buf_size = 0;
						}
					}
				}
			}
			else {
				codon_dt.inside_codon = !codon_dt.is_stop_codon(codon_dt.read_frame);
				codon_dt.update_codon_buf();

				if codon_dt.codon_buff_used >= 8 {
					codon_write_buf[codon_write_buf_size] = codon_dt.codon_buff as u8;
					codon_dt.codon_buff >>= 8;
					codon_dt.codon_buff_used -= 8;
					codon_write_buf_size += 1;
					if codon_write_buf_size == BUFFER_SIZE {
						_total_bytes_written += buf_writer(&mut code_bw, &codon_write_buf, codon_write_buf_size)?;
						codon_write_buf_size = 0;
					}
				}
			}
		}
		params.ln_read_buf.clear();
		buf_line_reader(&mut params.reader, &mut params.ln_read_buf)?;
	}
	if codon_dt.read_frame_used > 0 {
		if codon_dt.read_frame_used < codon_bits || !codon_dt.inside_codon {
			let tmp_remain = TEMP_BUFF_BITS - codon_dt.tmp_bits_used;
			let overflow = tmp_remain < codon_dt.read_frame_used;

			codon_dt.tmp += codon_dt.read_frame << codon_dt.tmp_bits_used;
			codon_dt.tmp_bits_used += tmp_remain;
			other_buf[other_buf_size] = codon_dt.tmp;
			codon_dt.reset_tmp_buf();
			other_buf_size += 1;
			if other_buf_size == BUFFER_SIZE {
				_total_bytes_written += buf_writer(&mut data_bw, &other_buf, other_buf_size)?;
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
			codon_dt.update_codon_buf();
		}
	}

	if codon_dt.tmp_bits_used > 0 {
		other_buf[other_buf_size] = codon_dt.tmp;
		codon_dt.reset_tmp_buf();
		other_buf_size += 1;
	}

	if other_buf_size > 0 {
		_total_bytes_written += buf_writer(&mut data_bw, &other_buf, other_buf_size)?;
		other_buf_size = 0;
	}
	while codon_dt.codon_buff_used >= 8 {
		codon_write_buf[codon_write_buf_size] = codon_dt.codon_buff as u8;
		codon_dt.codon_buff >>= 8;
		codon_dt.codon_buff_used -= 8;
		codon_write_buf_size += 1;
		if codon_write_buf_size == BUFFER_SIZE {
			_total_bytes_written += buf_writer(&mut code_bw, &codon_write_buf, codon_write_buf_size)?;
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
		_total_bytes_written += buf_writer(&mut code_bw, &codon_write_buf, codon_write_buf_size)?;
		codon_write_buf_size = 0;
	}
	if add_padding {
		code_bw.write_all(&PADDING_BUFF)?;
		_total_bytes_written += 4;
	}
	if params.data_fd.is_some() { data_bw.flush()?; }

	code_bw.flush()?;
	Ok(())
}

#[allow(dead_code)]
fn codon_simple_left_shift(params: &mut Params, table: &[u8; 256]) -> std::io::Result<()> {
    let mut codon_dt = CodonOnly::new(table, params.compile_encode, params.flags);
	codon_dt.dyn_init();
	let bit_mask: u8 = match codon_dt.bp_shift {
		1 => 0b00000_011,
		2 => 0b00_001111,
		_ => unreachable!(),
	};
	let add_padding: bool = params.flag_check(BitFlags::CodeSectionPadding);
	let codon_bits: u8 = codon_dt.bp_shift*3;
	let bp_lsb_offset: u8 = codon_dt.bp_shift*2;
	let mut code_bw: BufWriter<&File> = BufWriter::new(&params.code_fd);
	let mut data_bw = if params.data_fd.is_some() {
		BufWriter::new(params.data_fd.as_ref().unwrap())
	}
	else { BufWriter::new(&params.code_fd) };
    buf_line_reader(&mut params.reader, &mut params.ln_read_buf)?;
    let mut other_buf: [u8; BUFFER_SIZE] = [0; BUFFER_SIZE];
	let mut codon_write_buf: [u8; BUFFER_SIZE] = [0; BUFFER_SIZE];
    let mut codon_write_buf_size: usize = 0;
	let mut other_buf_size: usize = 0;
	let mut _total_bytes_written: usize = 0;

	if add_padding && params.read_buf_not_empty() {
		code_bw.write_all(&PADDING_BUFF)?;
		_total_bytes_written = 4;
	}
    while params.read_buf_not_empty() {
        for c in params.ln_read_buf.bytes() {
			codon_dt.read_frame <<= codon_dt.bp_shift;
			codon_dt.read_frame += codon_dt.bp_map[c as usize];
			codon_dt.read_frame_used += codon_dt.bp_shift;
			if codon_dt.read_frame_used < codon_bits { continue; }
			if !codon_dt.inside_codon {
				codon_dt.inside_codon = codon_dt.is_start_codon(codon_dt.read_frame);
				if !codon_dt.inside_codon {
					codon_dt.tmp <<= codon_dt.bp_shift;
					codon_dt.tmp += codon_dt.read_frame >> bp_lsb_offset;
					codon_dt.tmp_bits_used += codon_dt.bp_shift;
					codon_dt.read_frame &= bit_mask;
					codon_dt.read_frame_used -= codon_dt.bp_shift;
					 if codon_dt.tmp_bits_used == TEMP_BUFF_BITS {
						other_buf[other_buf_size] = codon_dt.tmp;
						codon_dt.reset_tmp_buf();
						other_buf_size += 1;
						if other_buf_size == BUFFER_SIZE {
							_total_bytes_written += buf_writer(&mut data_bw, &other_buf, other_buf_size)?;
							other_buf_size = 0;
						}
					}
				}
				else {
					codon_dt.update_codon_buf();
					if codon_dt.codon_buff_used >= 8 {
						codon_write_buf[codon_write_buf_size] = codon_dt.codon_buff as u8;
						codon_dt.codon_buff >>= 8;
						codon_dt.codon_buff_used -= 8;
						codon_write_buf_size += 1;
						if codon_write_buf_size == BUFFER_SIZE {
							_total_bytes_written += buf_writer(&mut code_bw, &codon_write_buf, codon_write_buf_size)?;
							codon_write_buf_size = 0;
						}
					}
				}
			}
			else {
                codon_dt.inside_codon = !codon_dt.is_stop_codon(codon_dt.read_frame);
				codon_dt.update_codon_buf();
				if codon_dt.codon_buff_used >= 8 {
					codon_write_buf[codon_write_buf_size] = codon_dt.codon_buff as u8;
					codon_dt.codon_buff >>= 8;
					codon_dt.codon_buff_used -= 8;
					codon_write_buf_size += 1;
					if codon_write_buf_size == BUFFER_SIZE {
						_total_bytes_written += buf_writer(&mut code_bw, &codon_write_buf, codon_write_buf_size)?;
						codon_write_buf_size = 0;
					}
				}
            }
        }
        params.ln_read_buf.clear();
        buf_line_reader(&mut params.reader, &mut params.ln_read_buf)?;
    }
	if codon_dt.read_frame_used > 0 {
		if codon_dt.read_frame_used < codon_bits || !codon_dt.inside_codon {
			let tmp_remain = TEMP_BUFF_BITS - codon_dt.tmp_bits_used;
			let overflow = tmp_remain < codon_dt.read_frame_used;
			let mut bp_mask: u8 = match codon_dt.bp_shift {
				1 => match codon_dt.read_frame_used {
						3 => 0b00000_111,
						2 => 0b00000_011,
						1 => 0b00000_001,
						_ => unreachable!(),
					},
				2 => match codon_dt.read_frame_used {
					6 => 0b00_111111,
					4 => 0b00_001111,
					2 => 0b00_000011,
					_ => unreachable!(),
				},
				_ => unreachable!(),
			};
			let mut rd_frame: u8 = codon_dt.read_frame;
			let mut offset = if overflow { tmp_remain } else { codon_dt.read_frame_used};
			while offset > 0 {
				rd_frame &= bp_mask;
				offset -= codon_dt.bp_shift;
				codon_dt.tmp += (rd_frame >> offset) << codon_dt.tmp_bits_used;
				codon_dt.tmp_bits_used += codon_dt.bp_shift;
				bp_mask >>= codon_dt.bp_shift;

			}

			other_buf[other_buf_size] = codon_dt.tmp;
			codon_dt.reset_tmp_buf();
			other_buf_size += 1;
			if other_buf_size == BUFFER_SIZE {
				_total_bytes_written += buf_writer(&mut data_bw, &other_buf, other_buf_size)?;
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
			codon_dt.update_codon_buf();
			//codon_dt.reset_read_frame();
		}
	}

	if codon_dt.tmp_bits_used > 0 {
		other_buf[other_buf_size] = codon_dt.tmp;
		codon_dt.reset_tmp_buf();
		other_buf_size += 1;
	}

	if other_buf_size > 0 {
		_total_bytes_written += buf_writer(&mut data_bw, &other_buf, other_buf_size)?;
		other_buf_size = 0;
	}

	while codon_dt.codon_buff_used >= 8 {
		codon_write_buf[codon_write_buf_size] = codon_dt.codon_buff as u8;
		codon_dt.codon_buff >>= 8;
		codon_dt.codon_buff_used -= 8;
		codon_write_buf_size += 1;
		if codon_write_buf_size == BUFFER_SIZE {
			_total_bytes_written += buf_writer(&mut code_bw, &codon_write_buf, codon_write_buf_size)?;
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
		_total_bytes_written += buf_writer(&mut code_bw, &codon_write_buf, codon_write_buf_size)?;
		codon_write_buf_size = 0;
	}
	if add_padding {
		code_bw.write_all(&PADDING_BUFF)?;
		_total_bytes_written += 4;
	}
	if params.data_fd.is_some() { data_bw.flush()?; }

	code_bw.flush()?;

    Ok(())
}

fn ms_bit_first_main(params: &mut Params, table:  &[u8; 256]) -> std::io::Result<()> {
	let bp_shift: u8 = match params.compile_encode {
		EncodingMethod::M3A0C1 | EncodingMethod::M3A1C0 | EncodingMethod::M3A0C1XOR | EncodingMethod::M3A1C0XOR => 1,
		_ => 2,
	};
	let add_padding:bool= params.flag_check(BitFlags::CodeSectionPadding);
	let mut code_bw: BufWriter<&File> = BufWriter::new(&params.code_fd);
    buf_line_reader(&mut params.reader, &mut params.ln_read_buf).unwrap_or_else(|e| { panic!("Failed to read line in open input file: {}", e); });

	let mut buff_2bit: [u8; BUFFER_SIZE] = [0; BUFFER_SIZE];
	let mut _buff_2bit_temp: u8 = 0;
    let mut _buff_2bit_size: usize = 0;
    let mut _temp_bits_used: u8 = 0;
    let mut _total_bytes_written: usize = 0;
	if add_padding && params.read_buf_not_empty() {
		code_bw.write_all(&PADDING_BUFF)?;
		_total_bytes_written = 4;
	}
    while params.read_buf_not_empty() {
        for c in params.ln_read_buf.bytes() {
			_buff_2bit_temp <<= bp_shift;
			_buff_2bit_temp += table[c as usize];
			_temp_bits_used += bp_shift;
			if _temp_bits_used == TEMP_BUFF_BITS {
				buff_2bit[_buff_2bit_size] = _buff_2bit_temp;
				_temp_bits_used = 0;
				_buff_2bit_temp = 0;
				_buff_2bit_size += 1;
				if _buff_2bit_size == BUFFER_SIZE {
					_total_bytes_written += buf_writer(&mut code_bw, &buff_2bit, _buff_2bit_size)?;
                    _buff_2bit_size = 0;
				}
			}
        }
        params.ln_read_buf.clear();
        buf_line_reader(&mut params.reader, &mut params.ln_read_buf)?;
    }
    if _buff_2bit_size == BUFFER_SIZE {
        _total_bytes_written += buf_writer(&mut code_bw, &buff_2bit, _buff_2bit_size)?;
        _buff_2bit_size = 0;
    }
    if _temp_bits_used > 0 {
        buff_2bit[_buff_2bit_size] = _buff_2bit_temp;
        _temp_bits_used = 0;
        _buff_2bit_temp = 0;
        _buff_2bit_size += 1;

    }
    if _buff_2bit_size > 0 {
        _total_bytes_written += buf_writer(&mut code_bw, &buff_2bit, _buff_2bit_size)?;
        _buff_2bit_size = 0;
    }
    let align = _total_bytes_written % 4; //align to 32-bit words
    if align != 0 {
        for i in 0..align as usize {
            buff_2bit[i] = CODE_PADDING;
            _buff_2bit_size += 1;
        }
        _total_bytes_written += buf_writer(&mut code_bw, &buff_2bit, _buff_2bit_size)?;
        _buff_2bit_size = 0;
    }
	if add_padding {
		code_bw.write_all(&PADDING_BUFF)?;
		_total_bytes_written += 4;

	}
    code_bw.flush()?;

    Ok(())
}

fn ls_bit_first_main(params: &mut Params, table: &[u8; 256]) -> std::io::Result<()> {
	let bp_shift: u8 = match params.compile_encode {
		EncodingMethod::M3A0C1 | EncodingMethod::M3A1C0 | EncodingMethod::M3A0C1XOR | EncodingMethod::M3A1C0XOR => 1,
		_ => 2,
	};
	let add_padding: bool = params.flag_check(BitFlags::CodeSectionPadding);

	let mut code_bw: BufWriter<&File> = BufWriter::new(&params.code_fd);
    buf_line_reader(&mut params.reader, &mut params.ln_read_buf).unwrap_or_else(|e| { panic!("Failed to read line in open input file: {}", e); });
    let mut _buff_2bit: [u8; BUFFER_SIZE] = [0; BUFFER_SIZE];
	let mut _buff_2bit_temp: u8 = 0;
    let mut _buff_2bit_size: usize = 0;
    let mut _temp_bits_used: u8 = 0;
    let mut _total_bytes_written: u64 = 0;
	if add_padding && params.read_buf_not_empty() {
		code_bw.write_all(&PADDING_BUFF)?;
		_total_bytes_written = 4;
	}

    while params.read_buf_not_empty() {
        for c in params.ln_read_buf.bytes() {
			_buff_2bit_temp += (table[c as usize]) << _temp_bits_used;
			_temp_bits_used += bp_shift;
			if _temp_bits_used == TEMP_BUFF_BITS {
				_buff_2bit[_buff_2bit_size] = _buff_2bit_temp;
				_temp_bits_used = 0;
				_buff_2bit_temp = 0;
				_buff_2bit_size += 1;
				if _buff_2bit_size == BUFFER_SIZE {
					_total_bytes_written += buf_writer(&mut code_bw, &_buff_2bit, _buff_2bit_size)? as u64;
                    _buff_2bit_size = 0;
				}
			}
        }
        params.ln_read_buf.clear();
        buf_line_reader(&mut params.reader, &mut params.ln_read_buf).unwrap_or_else(|e| { panic!("Failed to read line in open input file: {}", e); });
    }
    if _buff_2bit_size == BUFFER_SIZE {
        _total_bytes_written += buf_writer(&mut code_bw, &_buff_2bit, _buff_2bit_size)? as u64;
        _buff_2bit_size = 0;
    }
    if _temp_bits_used > 0 {
        _buff_2bit[_buff_2bit_size] = _buff_2bit_temp;
        _temp_bits_used = 0;
        _buff_2bit_temp = 0;
        _buff_2bit_size += 1;

    }
    if _buff_2bit_size > 0 {
        _total_bytes_written += buf_writer(&mut code_bw, &_buff_2bit, _buff_2bit_size)? as u64;
        _buff_2bit_size = 0;
    }
    let align = _total_bytes_written % 4; //align to 32-bit word size
    if align != 0 {
        for i in 0..align as usize {
            _buff_2bit[i] = CODE_PADDING;
            _buff_2bit_size += 1;
        }
        _total_bytes_written += buf_writer(&mut code_bw, &_buff_2bit, _buff_2bit_size)? as u64;
        _buff_2bit_size = 0;
    }

	if add_padding {
		code_bw.write_all(&PADDING_BUFF)?;
		_total_bytes_written += 4;
	}
    code_bw.flush()?;
    Ok(())
}




#[cfg(target_os="windows")]
#[allow(dead_code)]
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

#[cfg(target_os="linux")]
fn generate_elf_executable_code_only(code_len: u32, entry_byte_offset: u32) -> FileHeader {
	let phead_size: usize = linux::get_header_size(3, 0);
	let code_align = calc_alignment(phead_size as u32, SECTION_ALIGNMENT);

	let p_head_list: Vec<ELF32ProgHeaderEnt> = vec![ELF32ProgHeaderEnt::new(phead_size as u32), ELF32ProgHeaderEnt::with_params(1, 0, code_align, SECTION_ALIGNMENT), ELF32ProgHeaderEnt::with_params_flags(1, code_align, code_len+(2*entry_byte_offset), 0x05, SECTION_ALIGNMENT)];
	//let sec_head_list: Vec<ELF32SectHeaderEnt> = vec![ELF32SectHeaderEnt::default(), ELF32SectHeaderEnt::as_code(1, 0x1000, 0, code_len), ELF32SectHeaderEnt::as_shrtrtab(2, 0, )];
	FileHeader::with_defaults(ELF32Header::new(code_align+entry_byte_offset,3), p_head_list)
	//FileHeader::with_target_arch_abi(ELF32Header::new(code_align+entry_byte_offset, 3), p_head_list, 97, 40) //ARM RISC isa


}

#[cfg(target_os="linux")]
fn generate_elf_executable_header(sec_sizes: Vec<u32>, num_sections: u16, entry_byte_offset: u32) -> FileHeader {
	if num_sections == 1 {
		generate_elf_executable_code_only(sec_sizes[0], entry_byte_offset)
	}
	else {
		let phead_size: usize = linux::get_header_size(4, 0);
		let code_len = sec_sizes[0]+ 2*entry_byte_offset;
		let code_sec_offset = calc_alignment(phead_size as u32, SECTION_ALIGNMENT);
		let data_sec_offset = code_sec_offset + calc_alignment(code_len, SECTION_ALIGNMENT);
		let p_head_list: Vec<ELF32ProgHeaderEnt> = vec![ELF32ProgHeaderEnt::new(phead_size as u32), ELF32ProgHeaderEnt::with_params(1, 0, code_sec_offset, SECTION_ALIGNMENT), ELF32ProgHeaderEnt::with_params_flags(1, code_sec_offset, code_len, 0x05, SECTION_ALIGNMENT), ELF32ProgHeaderEnt::with_params_flags(1, data_sec_offset, sec_sizes[1], 0x06, SECTION_ALIGNMENT)];
		//let sec_head_list: Vec<ELF32SectHeaderEnt> = vec![ELF32SectHeaderEnt::default(), ELF32SectHeaderEnt::with_params_flags(1, 1, 0x06, 0x1000, code_len, 0, 0, 16, 0)];
		FileHeader::with_defaults(ELF32Header::new(code_sec_offset+entry_byte_offset,4), p_head_list)
		//FileHeader::with_target_arch_abi(ELF32Header::new(code_sec_offset+entry_byte_offset, 4), p_head_list, 97, 40) //ARM RISC isa
	}



}

#[cfg(target_os="linux")]
fn linux_elf_writer(params: &mut Params)  -> std::io::Result<()> {
	const ALIGN_SECTION: usize = SECTION_ALIGNMENT as usize;
	let mut f_sizes: Vec<u32> = Vec::new();
	let code_size: u64 = params.code_fd.metadata()?.len();
	let mut temp_reader: BufReader<&File> = BufReader::new(&params.code_fd);
	let entry_offset: u32 = if params.flag_check(BitFlags::CodeSectionPadding) { WORD_ALIGN_BYTES } else { 0 };
	temp_reader.rewind()?;
	f_sizes.push(code_size as u32);
	let f_data_valid = params.flag_check(BitFlags::UseCodons) && params.data_fd.is_some();
	if f_data_valid {
		f_sizes.push(params.data_fd.as_ref().unwrap().metadata()?.len() as u32);
	}

	let header = generate_elf_executable_header(f_sizes.clone(), f_sizes.len() as u16, entry_offset);
	let header_buff = header.serialize();
	let mut _bytes_written: usize = buf_writer(&mut params.writer, &header_buff.as_slice(), header_buff.len())?;
	let buff_align: usize = ALIGN_SECTION - (_bytes_written % ALIGN_SECTION);
	let mut write_buff: [u8; ALIGN_SECTION] = [HEADER_PADDING_BYTE; ALIGN_SECTION];
	if buff_align < ALIGN_SECTION {
		_bytes_written += buf_writer(&mut params.writer, &write_buff, buff_align)?;
	}
	params.writer.flush()?;
	for _full_reads in 0..f_sizes[0]/SECTION_ALIGNMENT {
		temp_reader.read_exact(&mut write_buff).unwrap_or_else(|e| { panic!("code.bin file reader error: {}", e); });
		_bytes_written += buf_writer(&mut params.writer, &write_buff,ALIGN_SECTION)?;
	}
	let mut end_vec: Vec<u8> = Vec::new();
	let code_data_end = temp_reader.read_to_end(&mut end_vec)?;
	_bytes_written += buf_writer(&mut params.writer, end_vec.as_slice(), end_vec.len())?;
	let align = ALIGN_SECTION - (code_data_end % ALIGN_SECTION);

	if align > 0 {
		write_buff = [CODE_PADDING; ALIGN_SECTION];
		_bytes_written += buf_writer(&mut params.writer, &write_buff, align)?;
	}
	params.writer.flush()?;
	if f_data_valid {
		let mut data_reader: BufReader<&File> = BufReader::new(params.data_fd.as_mut().unwrap());
		data_reader.rewind()?;
		write_buff = [CODE_PADDING; ALIGN_SECTION];
		for _full_buf in 0..f_sizes[1]/SECTION_ALIGNMENT {
			data_reader.read_exact(&mut write_buff).unwrap_or_else(|e| { panic!("data.bin file reader error: {}", e); });
			_bytes_written += buf_writer(&mut params.writer, &write_buff,ALIGN_SECTION)?;
		}
		params.writer.flush()?;
		let mut end_vec: Vec<u8> = Vec::new();
		let data_end = data_reader.read_to_end(&mut end_vec)?;
		_bytes_written += buf_writer(&mut params.writer, end_vec.as_slice(), end_vec.len())?;
		let align_data = ALIGN_SECTION - (data_end % ALIGN_SECTION);
		if align_data > 0 {
			write_buff = [CODE_PADDING; ALIGN_SECTION];
			_bytes_written += buf_writer(&mut params.writer, &write_buff, align_data)?;
		}
		params.writer.flush()?;
	}
	Ok(())
}

#[cfg(target_os="windows")]
fn win_exe_writer(params: &mut Params)  -> std::io::Result<()> {
	const ALIGN_FILE: usize = FILE_ALIGNMENT as usize;
	let mut f_sizes: Vec<u32> = Vec::new();
	let code_size: u64 = params.code_fd.metadata()?.len();
	let mut temp_reader: BufReader<&File> = BufReader::new(&params.code_fd);
	let entry_offset: u32 = if params.flag_check(BitFlags::CodeSectionPadding) { WORD_ALIGN_BYTES } else { 0 };
	temp_reader.rewind()?;
	f_sizes.push(code_size as u32);
	let f_data_valid = params.flag_check(BitFlags::UseCodons) && params.data_fd.is_some();
	if f_data_valid {
		f_sizes.push(params.data_fd.as_ref().unwrap().metadata()?.len() as u32);
	}

	let header = win::generate_win_executable_header(f_sizes.clone(), f_sizes.len() as u16, entry_offset);
	let header_buff = header.serialize();
	let mut _bytes_written: usize = buf_writer(&mut params.writer, &header_buff.as_slice(), header_buff.len())?;
	let buff_align: usize = ALIGN_FILE - (_bytes_written % ALIGN_FILE);
	let mut write_buff: [u8; ALIGN_FILE] = [HEADER_PADDING_BYTE; ALIGN_FILE];
	if buff_align < ALIGN_FILE {
		_bytes_written += buf_writer(&mut params.writer, &write_buff, buff_align)?;
	}
	params.writer.flush()?;
	for _full_reads in 0..f_sizes[0]/FILE_ALIGNMENT {
		temp_reader.read_exact(&mut write_buff).unwrap_or_else(|e| { panic!("code.bin file reader error: {}", e); });
		_bytes_written += buf_writer(&mut params.writer, &write_buff,ALIGN_FILE)?;
	}
	let mut end_vec: Vec<u8> = Vec::new();
	let code_data_end = temp_reader.read_to_end(&mut end_vec)?;
	_bytes_written += buf_writer(&mut params.writer, end_vec.as_slice(), end_vec.len())?;
	let align = ALIGN_FILE - (code_data_end % ALIGN_FILE);

	if align > 0 {
		write_buff = [CODE_PADDING; ALIGN_FILE];
		_bytes_written += buf_writer(&mut params.writer, &write_buff, align)?;
	}
	params.writer.flush()?;
	if f_data_valid {
		let mut data_reader: BufReader<&File> = BufReader::new(params.data_fd.as_mut().unwrap());
		data_reader.rewind()?;
		write_buff = [CODE_PADDING; ALIGN_FILE];
		for _full_buf in 0..f_sizes[1]/FILE_ALIGNMENT {
			data_reader.read_exact(&mut write_buff).unwrap_or_else(|e| { panic!("data.bin file reader error: {}", e); });
			_bytes_written += buf_writer(&mut params.writer, &write_buff,ALIGN_FILE)?;
		}
		params.writer.flush()?;
		let mut end_vec: Vec<u8> = Vec::new();
		let data_end = data_reader.read_to_end(&mut end_vec)?;
		_bytes_written += buf_writer(&mut params.writer, end_vec.as_slice(), end_vec.len())?;
		let align_data = ALIGN_FILE - (data_end % ALIGN_FILE);
		if align_data > 0 {
			write_buff = [CODE_PADDING; ALIGN_FILE];
			_bytes_written += buf_writer(&mut params.writer, &write_buff, align_data)?;
		}
		params.writer.flush()?;
	}

	Ok(())
}
fn generate_executable(params: &mut Params) -> std::io::Result<()> {
	#[cfg(target_os="windows")]
	return win_exe_writer(params);
	#[cfg(target_os="linux")]
	return linux_elf_writer(params);
}
fn remove_temp_files(temp_dir_p: &str, temp_files: Vec<&str>) -> std::io::Result<()> {
	for item in temp_files {
		let mut fp = PathBuf::from(temp_dir_p);
		fp.push(item);
		if fp.exists() && fp.is_file() { std::fs::remove_file(fp.as_path())?; }
	}
	Ok(())
}

fn gen_binary_table(bp_map_val: [u8; 4]) -> [u8; 256] {
	 const BP_MAP: [usize; 8] = ['A' as usize, 'a' as usize, 'C' as usize, 'c' as usize, 'G' as usize, 'g' as usize, 'T' as usize, 't' as usize,];
	let mut ret: [u8; 256] = [255; 256];
	ret[BP_MAP[0]] = bp_map_val[0];
    ret[BP_MAP[1]] = bp_map_val[0];
    ret[BP_MAP[2]] = bp_map_val[1];
    ret[BP_MAP[3]] = bp_map_val[1];
    ret[BP_MAP[4]] = bp_map_val[2];
    ret[BP_MAP[5]] = bp_map_val[2];
    ret[BP_MAP[6]] = bp_map_val[3];
    ret[BP_MAP[7]] = bp_map_val[3];
	ret

}
fn main() -> std::io::Result<()> {
    let main_args: Vec<String> = env::args().collect();
    let (check, temp_path) = parse_main_args(&main_args)?;
    if !check.check_files() { return Ok(()); }
	if !check.is_req_valid() { Error::new(std::io::ErrorKind::InvalidData, "Parameters are invalid!"); }

	let mut params: Params = Params::from(check);
    let table: [u8; 256] = gen_binary_table(params.bp_val_map());
    match params.flags {
		0b00_000|0b00_001 => ls_bit_first_main(&mut params, &table)?,
		0b00_010|0b00_011 => ms_bit_first_main(&mut params, &table)?,
		0b00_100|0b00_101 => codon_lsb_shift(&mut params, &table)?,
		0b00_110|0b00_111 => codon_simple_left_shift(&mut params, &table)?,
		_ => unreachable!(),
	}
	params.code_fd.sync_all()?;
	if params.flag_check(BitFlags::UseCodons) && params.data_fd.is_some() { //codons
		params.data_fd.as_mut().unwrap().sync_all()?;
	}
	generate_executable(&mut params)?;
	return remove_temp_files(temp_path.as_str(),vec!["code.bin", "data.bin"]);


}

fn buf_writer<W: Write>(writer: &mut W, buf: &[u8], buff_len: usize) -> std::io::Result<usize> {
	let mut _written_bytes: usize = 0;
	while _written_bytes < buff_len {
		let n: usize = writer.write(&buf[_written_bytes..buff_len])?;
		_written_bytes += n;
	}
	Ok(_written_bytes)
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
