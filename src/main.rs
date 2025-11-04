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
use std::path::{Path, PathBuf};
mod arg_parser;
pub mod executable_header;
use arg_parser::*;

use crate::executable_header::FILE_ALIGNMENT;
use crate::executable_header::HeaderPE;
use crate::executable_header::ImageDataDirectory;
use crate::executable_header::ImageNTHeaders32;
use crate::executable_header::ImageSectionHeaderEnt;
use crate::executable_header::InitFileHeader;
use crate::executable_header::OptionalImageFileHeader32;
use crate::executable_header::SECTION_ALIGNMENT;
use crate::executable_header::ToBytes;
/*
enum bp_map_default{ // Y = (X & 0x06) >> 1
    A=0b00, // 'A'=01000|00|1
    C=0b01, // 'C'=01000|01|1
    T=0b10, // 'T'=01010|10|0
    G=0b11, // 'G'=01000|11|1
} */
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
	pub codon_buff: u32,
	pub codon_buff_remain: u8,
	pub flags: u8,
	pub inside_codon: bool
}

impl CodonOnly {
	pub fn new(table: &[u8; 256], encoding: EncodingMethod, b_flags: u8) ->Self {
		CodonOnly { codon_map: HashMap::with_capacity(64), bp_map: *table, stop_codons: [0; 3], start_codon: 0,
			bp_shift: (match encoding {
		EncodingMethod::M3A0C1 | EncodingMethod::M3A1C0 | EncodingMethod::M3A0C1XOR | EncodingMethod::M3A1C0XOR => 1,
		_ => 2,
	}), read_frame: 0, read_frame_used: 0, tmp: 0, tmp_bits_used: 0, codon_buff: 0, codon_buff_remain: 0, flags: b_flags, inside_codon: false }
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
}


/***
 *
 */
pub struct Params {
    reader: BufReader<File>,
    writer: BufWriter<File>,
	code_fd: File,
	ln_read_buf: String,
    compile_encode: EncodingMethod,
    ifile_ext: InputFileFormat,
	flags: u8, // b0: 0= lsb offset shifting, 1= simple fixed bit shift, b1=codons?, b2=in-place(1) or contiguous(0) data/code for codons. b3=codon inside = code(0) or data(1)
}
impl Params {
    pub fn from(ck: ParamsCK) -> Params {
        Params {reader: ck.reader.unwrap(), writer: ck.writer.unwrap(), code_fd: OpenOptions::new().write(true).create(true).read(true).open(ck.code_fp.unwrap().as_os_str()).unwrap(), ln_read_buf: String::new(), compile_encode: ck.compile_encode, ifile_ext: ck.ifile_ext, flags: ck.flags }
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
const NOOP_CODE: u8 = 0x90;
const INT3_CODE: u8 = 0xCC;
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

	Ok(())
}


fn codon_msb_simple_left_shift(mut params: &mut Params, mut table: &mut [u8; 256]) -> std::io::Result<()> {
    const BP_MAP: [usize; 8] = ['A' as usize, 'a' as usize, 'C' as usize, 'c' as usize, 'G' as usize, 'g' as usize, 'T' as usize, 't' as usize];
	let mut codon_dt = CodonOnly::new(table, params.compile_encode, params.flags);
	codon_dt.dyn_init();
	let mut code_bw: BufWriter<&File> = BufWriter::new(&params.code_fd);
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

       (buff_2bit_temp, temp_bits_used, read_frame,read_frame_used) = eof_8_mod_6(buff_2bit_temp, temp_bits_used, read_frame, read_frame_used);

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
        let align_fill: u8 = if  read_frame_used == 6 && codon_dt.is_stop_codon(read_frame) { read_frame } else { INT3_CODE };

        if align_fill != INT3_CODE {
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
            codon_write_buf[codon_write_buf_size+a] = INT3_CODE;
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
    let align = total_bytes_written % (size_of::<usize>() as u64); //align to system default word bit size. Expect 64 most of the time
    if align != 0 {
        for i in 0..align as usize {
            buff_2bit[i] = NOOP_CODE;
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
    let align = total_bytes_written % (size_of::<usize>() as u64); //align to system default word bit size. Expect 64 most of the time
    if align != 0 {
        for i in 0..align as usize {
            buff_2bit[i] = NOOP_CODE;
            buff_2bit_size += 1;
        }
        total_bytes_written += buf_writer(&mut code_bw, &buff_2bit, buff_2bit_size)? as u64;
        buff_2bit_size = 0;
    }
    code_bw.flush()?;
    Ok(())
}

const fn get_header_size(num_sections: usize) -> usize {
	use std::mem::size_of;
	size_of::<InitFileHeader>() + size_of::<ImageNTHeaders32>() + (size_of::<ImageSectionHeaderEnt>() * num_sections)
}

fn generate_executable(mut params: &mut Params) -> std::io::Result<()> {
	let code_size: u64 = params.code_fd.metadata()?.len();
	let mut temp_reader: BufReader<&File> = BufReader::new(&params.code_fd);
	temp_reader.rewind()?;
	let m32_code_size = code_size as u32;
	//println!("code size: {code_size} == {m32_code_size}");
	let header_data_len: u32 = get_header_size(1) as u32;
	let header_align: u32 = (match header_data_len % FILE_ALIGNMENT {
		0 => header_data_len/FILE_ALIGNMENT,
		_ => (header_data_len/FILE_ALIGNMENT) + 1,
	}) * FILE_ALIGNMENT;
	let va_header_offset: u32 = (match header_data_len % SECTION_ALIGNMENT {
		0 => header_data_len/SECTION_ALIGNMENT,
		_ => (header_data_len/SECTION_ALIGNMENT) + 1,
	}) * SECTION_ALIGNMENT;
	let va_data_addr_offset: u32 = va_header_offset + (match m32_code_size % SECTION_ALIGNMENT {
		0 => m32_code_size/SECTION_ALIGNMENT,
		_ => (m32_code_size/SECTION_ALIGNMENT) + 1,
	}) * SECTION_ALIGNMENT;

	let data_dir_list: [ImageDataDirectory; 16] = [
		ImageDataDirectory::default(), ImageDataDirectory::default(), ImageDataDirectory::default(),
		ImageDataDirectory::default(), ImageDataDirectory::default(),
		/*ImageDataDirectory::with_parameters(va_reloc_addr_offset, SECTION_ALIGNMENT),*/ ImageDataDirectory::default(), //base addr relocation table
		ImageDataDirectory::default(),
		ImageDataDirectory::default(), //Architecture Specific Directory
		ImageDataDirectory::default(), //X86 specific field
		ImageDataDirectory::default(), ImageDataDirectory::default(),
		ImageDataDirectory::default(), //load configuration directory
		ImageDataDirectory::default(), ImageDataDirectory::default(), ImageDataDirectory::default(), ImageDataDirectory::default(),
		];
	let opt_file_header: OptionalImageFileHeader32 = OptionalImageFileHeader32::with_code_only(m32_code_size, va_header_offset, va_header_offset, va_data_addr_offset, data_dir_list);
	let nt_header_32: ImageNTHeaders32 = ImageNTHeaders32::with_defaults(1, opt_file_header);
	let header = HeaderPE { init_header: InitFileHeader::default(), nt_header: nt_header_32, section_header: vec![ImageSectionHeaderEnt::with_perms(*b".text\0\0\0",va_data_addr_offset-va_header_offset, va_header_offset, m32_code_size, va_header_offset, 0x60000020)]};
	let header_buff = header.bufferize();
	let mut bytes_written: usize = buf_writer(&mut params.writer, &header_buff.as_slice(), header_buff.len())?;
	let buff_align: usize = FILE_ALIGNMENT as usize - (bytes_written % FILE_ALIGNMENT as usize);
	let mut write_buff: [u8; SECTION_ALIGNMENT as usize] = [0x00; SECTION_ALIGNMENT as usize];
	if buff_align < FILE_ALIGNMENT as usize {
		bytes_written += buf_writer(&mut params.writer, &write_buff, buff_align)?;
	}
	params.writer.flush()?;
	for full_reads in 0..m32_code_size/SECTION_ALIGNMENT {
		temp_reader.read_exact(&mut write_buff)?;
		bytes_written += buf_writer(&mut params.writer, &write_buff, SECTION_ALIGNMENT as usize)?
	}
	let mut end_vec: Vec<u8> = Vec::with_capacity(SECTION_ALIGNMENT as usize);
	let code_data_end = temp_reader.read_to_end(&mut end_vec)?;
	bytes_written += buf_writer(&mut params.writer, &end_vec.as_slice(), end_vec.len())?;
	if code_data_end % SECTION_ALIGNMENT as usize > 0 {
		write_buff = [0x00; SECTION_ALIGNMENT as usize];
		bytes_written += buf_writer(&mut params.writer, &write_buff, SECTION_ALIGNMENT as usize -(code_data_end % SECTION_ALIGNMENT as usize))?;
	}
	params.writer.flush()?;

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
    if !check.check_files() {
        return Ok(());
    }
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
		0 => ls_bit_first_main(&mut params, &mut table)?,
		1 => ms_bit_first_main(&mut params, &mut table)?,
		_ => unreachable!(),
	}
	params.code_fd.sync_all()?;
	generate_executable(&mut params)?;
	return remove_temp_files(temp_path.as_str(),vec!["code.bin"]);


}

#[inline(always)]
fn method1_ascii_to_2bit_convert(c: u8) -> u8 { (c & 0x06) >> 1 }

fn eof_8_mod_6(mut tmp: u8, mut tmp_bits_used: u8, mut frame: u8, mut frame_bits: u8) -> (u8, u8, u8, u8) {
    match tmp_bits_used {
        6 => { tmp = (tmp << 2) + (frame >> 4);
            tmp_bits_used += 2;
            frame &= CODON_BIT_MASK >> 2;
            frame_bits -= 2;
        },
        4 => { tmp = (tmp << 4) + (frame >> 2);
            tmp_bits_used += 4;
            frame &= CODON_BIT_MASK >> 4;
            frame_bits -= 4;
        },
        2 => { tmp = (tmp << 6) + frame;
            tmp_bits_used += 6;
            frame = 0;
            frame_bits = 0;
        },
        0 => { tmp = frame;
            tmp_bits_used = 6;
            frame = 0;
            frame_bits = 0;
        },
        _ => unreachable!(), //vals >=8, 7, 5, 3 and 1 are unreachable
    }
    (tmp, tmp_bits_used, frame, frame_bits)
}

// fn buf_writer_eof_out_codon<W: Write>(writer: &mut W, buf: &mut [u8], mut buf_used: usize, mut temp: u8, mut temp_bits: u8, total_bytes: u64) -> std::io::Result<u64>{
//     let mut last_buff_write: u64 = 0;

//     Ok(last_buff_write)
// }

// fn buf_writer_eof_in_codon<W: Write>(writer: &mut W, buf: &mut [u8], mut buf_used: usize, mut frame: u8, mut frame_bits: u8, mut codon_word: u64, mut codon_bits_remain: u8) -> std::io::Result<u64>{
//     let mut last_buff_write: u64 = 0;

//     Ok(last_buff_write)
// }

fn buf_64_mod_6_align<W: Write>(writer: &mut W, buf: &mut [u8], buff_len: &mut usize, mut frame: u8, mut frame_bits: u8, codon_word: &mut u64, codon_bits_remain: &mut u8) -> std::io::Result<u64>{
	let mut ret: u64 = 0;
	let mut len: usize = *buff_len;
	if *codon_bits_remain >= frame_bits {
		*codon_word <<= frame_bits;
		*codon_word += frame as u64;
		*codon_bits_remain -= frame_bits;
		frame_bits = 0;
		frame = 0;
	}
	else {
		*codon_word <<= *codon_bits_remain;
		*codon_word += (frame >> (frame_bits-*codon_bits_remain)) as u64;
		let mask: u8 = (!((frame >> (frame_bits-*codon_bits_remain)) << (frame_bits-*codon_bits_remain))) & CODON_BIT_MASK;
		frame &= mask;
		frame_bits -= *codon_bits_remain;
		*codon_bits_remain = 0;
	}
	if len+8 > BUFFER_SIZE  {
		ret = buf_writer(writer, buf, len)? as u64;
		len = 0;
	}
	if *codon_bits_remain == 0 {
		for off in 7..=0 {
			buf[len+off] = *codon_word as u8;
			*codon_word >>= 8;
		}
		*codon_word = 0;
		len += 8;
		*codon_bits_remain = 64;
	}
	if frame_bits > 0 {
		*codon_word = frame as u64;
		frame = 0;
		*codon_bits_remain -= frame_bits;
		frame_bits = 0;
	}
	*buff_len = len;
	Ok(ret)
}

fn buf_writer_align<W: Write>(writer: &mut W, buf: &mut [u8], mut buff_len: usize, temp: u8, tmp_len: u8, total_bytes: u64, align_val: u8) -> std::io::Result<u64> {
	let mut write_align_cnt: u64 = 0;
	if tmp_len > 0 {
		if buff_len >= BUFFER_SIZE {
			write_align_cnt += buf_writer(writer, buf, buff_len)? as u64;
			buff_len = 0;
		}
		buf[buff_len] = temp;
		buff_len += 1;
	}
	// let align = (total_bytes + write_align_cnt + buff_len as u64) % (size_of::<usize>() as u64);
	// if (align + buff_len as u64) > (BUFFER_SIZE as u64) {
	// 	write_align_cnt += buf_writer(writer, buf, buff_len)? as u64;
	// 	buff_len = 0;
	// }
	// for _ in 0..align {
	// 	buf[buff_len] = align_val;//data gets \0 char, code gets NO OP
	// 	buff_len += 1;
	// }
	write_align_cnt += buf_writer(writer, buf, buff_len)? as u64;
	Ok(write_align_cnt)
}

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
