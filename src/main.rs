use std::collections::HashMap;
use std::env;
use std::fs::File;
use std::io::BufWriter;
use std::io::Error;
use std::io::Write;
use std::vec;
mod arg_parser;
use arg_parser::*;
mod dna_compiler;
use dna_compiler::*;
mod executable_header;
use executable_header::*;
mod file_handler;
use file_handler::*;
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

	pub fn is_stop_codon(&self) -> bool {
		let ck: u8 = self.codon_map[&self.read_frame];
    	ck == self.stop_codons[0] || ck == self.stop_codons[1] || ck == self.stop_codons[2]
	}

	pub fn is_start_codon(&self) -> bool {
		let ck: u8 = self.codon_map[&self.read_frame];
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

}


const BUFFER_SIZE: usize = 1024;
const TEMP_BUFF_BITS: u8 = u8::BITS as u8;
const PADDING_BUFF: [u8; 4] = [CODE_PADDING; 4];

#[repr(C)]
pub struct StoreBuf {
	pub buf: [u8; BUFFER_SIZE],
	pub buf_used: usize
}

impl StoreBuf {
	pub fn push(&mut self, val: u8) {
		self.buf[self.buf_used] = val;
		self.buf_used += 1;
	}
	pub fn is_full(&self) -> bool { self.buf_used == BUFFER_SIZE }
	pub fn not_empty(&self) -> bool { self.buf_used > 0 }
	pub fn write_buf<W: Write>(&mut self, writer: &mut W) -> std::io::Result<usize> {
		let ret = buf_writer(writer, &self.buf, self.buf_used)?;
		self.buf_used = 0;
		Ok(ret)
	}
}

impl Default for StoreBuf {
	fn default() -> Self {
		StoreBuf { buf: [0; BUFFER_SIZE], buf_used: 0 }
	}
}


fn codon_ls_bp_shift(params: &mut Params, table: &[u8; 256]) -> std::io::Result<()> {
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
	let mut _other_buf: [u8; BUFFER_SIZE] = [0; BUFFER_SIZE];
	let mut _codon_write_buf: [u8; BUFFER_SIZE] = [0; BUFFER_SIZE];
    let mut _codon_write_buf_size: usize = 0;
	let mut _other_buf_size: usize = 0;
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
				codon_dt.inside_codon = codon_dt.is_start_codon();
				if !codon_dt.inside_codon {
					codon_dt.tmp += (codon_dt.read_frame & bit_mask) << codon_dt.tmp_bits_used;
					codon_dt.tmp_bits_used += codon_dt.bp_shift;
					codon_dt.read_frame >>= codon_dt.bp_shift;
					codon_dt.read_frame_used -= codon_dt.bp_shift;
					if codon_dt.tmp_bits_used == TEMP_BUFF_BITS {
						_other_buf[_other_buf_size] = codon_dt.tmp;
						codon_dt.reset_tmp_buf();
						_other_buf_size += 1;
						if _other_buf_size == BUFFER_SIZE {
							_total_bytes_written += buf_writer(&mut data_bw, &_other_buf, _other_buf_size)?;
							_other_buf_size = 0;
						}
					}
				}
				else {
					codon_dt.update_codon_buf();
					if codon_dt.codon_buff_used >= 8 {
						_codon_write_buf[_codon_write_buf_size] = codon_dt.codon_buff as u8;
						codon_dt.codon_buff >>= 8;
						codon_dt.codon_buff_used -= 8;
						_codon_write_buf_size += 1;
						if _codon_write_buf_size == BUFFER_SIZE {
							_total_bytes_written += buf_writer(&mut code_bw, &_codon_write_buf, _codon_write_buf_size)?;
							_codon_write_buf_size = 0;
						}
					}
				}
			}
			else {
				codon_dt.inside_codon = !codon_dt.is_stop_codon();
				codon_dt.update_codon_buf();

				if codon_dt.codon_buff_used >= 8 {
					_codon_write_buf[_codon_write_buf_size] = codon_dt.codon_buff as u8;
					codon_dt.codon_buff >>= 8;
					codon_dt.codon_buff_used -= 8;
					_codon_write_buf_size += 1;
					if _codon_write_buf_size == BUFFER_SIZE {
						_total_bytes_written += buf_writer(&mut code_bw, &_codon_write_buf, _codon_write_buf_size)?;
						_codon_write_buf_size = 0;
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
			_other_buf[_other_buf_size] = codon_dt.tmp;
			codon_dt.reset_tmp_buf();
			_other_buf_size += 1;
			if _other_buf_size == BUFFER_SIZE {
				_total_bytes_written += buf_writer(&mut data_bw, &_other_buf, _other_buf_size)?;
				_other_buf_size = 0;
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
		_other_buf[_other_buf_size] = codon_dt.tmp;
		codon_dt.reset_tmp_buf();
		_other_buf_size += 1;
	}

	if _other_buf_size > 0 {
		_total_bytes_written += buf_writer(&mut data_bw, &_other_buf, _other_buf_size)?;
		_other_buf_size = 0;
	}
	while codon_dt.codon_buff_used >= 8 {
		_codon_write_buf[_codon_write_buf_size] = codon_dt.codon_buff as u8;
		codon_dt.codon_buff >>= 8;
		codon_dt.codon_buff_used -= 8;
		_codon_write_buf_size += 1;
		if _codon_write_buf_size == BUFFER_SIZE {
			_total_bytes_written += buf_writer(&mut code_bw, &_codon_write_buf, _codon_write_buf_size)?;
			_codon_write_buf_size = 0;
		}
	}
	if codon_dt.codon_buff_used > 0 {
		_codon_write_buf[_codon_write_buf_size] = codon_dt.codon_buff as u8;
		codon_dt.codon_buff >>= codon_dt.codon_buff_used;
		codon_dt.codon_buff_used = 0;
		_codon_write_buf_size += 1;
	}
	if _codon_write_buf_size > 0 {
		_total_bytes_written += buf_writer(&mut code_bw, &_codon_write_buf, _codon_write_buf_size)?;
		_codon_write_buf_size = 0;
	}
	if add_padding {
		code_bw.write_all(&PADDING_BUFF)?;
		_total_bytes_written += 4;
	}
	if params.data_fd.is_some() { data_bw.flush()?; }

	code_bw.flush()?;
	Ok(())
}

fn codon_bp_fr_shift(params: &mut Params, table: &[u8; 256]) -> std::io::Result<()> {
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
    let mut _other_buf: [u8; BUFFER_SIZE] = [0; BUFFER_SIZE];
	let mut _codon_write_buf: [u8; BUFFER_SIZE] = [0; BUFFER_SIZE];
    let mut _codon_write_buf_size: usize = 0;
	let mut _other_buf_size: usize = 0;
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
				codon_dt.inside_codon = codon_dt.is_start_codon();
				if !codon_dt.inside_codon {
					codon_dt.tmp <<= codon_dt.bp_shift;
					codon_dt.tmp += codon_dt.read_frame >> bp_lsb_offset;
					codon_dt.tmp_bits_used += codon_dt.bp_shift;
					codon_dt.read_frame &= bit_mask;
					codon_dt.read_frame_used -= codon_dt.bp_shift;
					 if codon_dt.tmp_bits_used == TEMP_BUFF_BITS {
						_other_buf[_other_buf_size] = codon_dt.tmp;
						codon_dt.reset_tmp_buf();
						_other_buf_size += 1;
						if _other_buf_size == BUFFER_SIZE {
							_total_bytes_written += buf_writer(&mut data_bw, &_other_buf, _other_buf_size)?;
							_other_buf_size = 0;
						}
					}
				}
				else {
					codon_dt.update_codon_buf();
					if codon_dt.codon_buff_used >= 8 {
						_codon_write_buf[_codon_write_buf_size] = codon_dt.codon_buff as u8;
						codon_dt.codon_buff >>= 8;
						codon_dt.codon_buff_used -= 8;
						_codon_write_buf_size += 1;
						if _codon_write_buf_size == BUFFER_SIZE {
							_total_bytes_written += buf_writer(&mut code_bw, &_codon_write_buf, _codon_write_buf_size)?;
							_codon_write_buf_size = 0;
						}
					}
				}
			}
			else {
                codon_dt.inside_codon = !codon_dt.is_stop_codon();
				codon_dt.update_codon_buf();
				if codon_dt.codon_buff_used >= 8 {
					_codon_write_buf[_codon_write_buf_size] = codon_dt.codon_buff as u8;
					codon_dt.codon_buff >>= 8;
					codon_dt.codon_buff_used -= 8;
					_codon_write_buf_size += 1;
					if _codon_write_buf_size == BUFFER_SIZE {
						_total_bytes_written += buf_writer(&mut code_bw, &_codon_write_buf, _codon_write_buf_size)?;
						_codon_write_buf_size = 0;
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

			_other_buf[_other_buf_size] = codon_dt.tmp;
			codon_dt.reset_tmp_buf();
			_other_buf_size += 1;
			if _other_buf_size == BUFFER_SIZE {
				_total_bytes_written += buf_writer(&mut data_bw, &_other_buf, _other_buf_size)?;
				_other_buf_size = 0;
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
		_other_buf[_other_buf_size] = codon_dt.tmp;
		codon_dt.reset_tmp_buf();
		_other_buf_size += 1;
	}

	if _other_buf_size > 0 {
		_total_bytes_written += buf_writer(&mut data_bw, &_other_buf, _other_buf_size)?;
		_other_buf_size = 0;
	}

	while codon_dt.codon_buff_used >= 8 {
		_codon_write_buf[_codon_write_buf_size] = codon_dt.codon_buff as u8;
		codon_dt.codon_buff >>= 8;
		codon_dt.codon_buff_used -= 8;
		_codon_write_buf_size += 1;
		if _codon_write_buf_size == BUFFER_SIZE {
			_total_bytes_written += buf_writer(&mut code_bw, &_codon_write_buf, _codon_write_buf_size)?;
			_codon_write_buf_size = 0;
		}
	}

	if codon_dt.codon_buff_used > 0 {
		_codon_write_buf[_codon_write_buf_size] = codon_dt.codon_buff as u8;
		codon_dt.codon_buff >>= codon_dt.codon_buff_used;
		codon_dt.codon_buff_used = 0;
		_codon_write_buf_size += 1;
	}
	if _codon_write_buf_size > 0 {
		_total_bytes_written += buf_writer(&mut code_bw, &_codon_write_buf, _codon_write_buf_size)?;
		_codon_write_buf_size = 0;
	}
	if add_padding {
		code_bw.write_all(&PADDING_BUFF)?;
		_total_bytes_written += 4;
	}
	if params.data_fd.is_some() { data_bw.flush()?; }

	code_bw.flush()?;

    Ok(())
}

fn ms_bp_first(params: &mut Params, table:  &[u8; 256]) -> std::io::Result<()> {
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

fn ls_bp_first(params: &mut Params, table: &[u8; 256]) -> std::io::Result<()> {
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
    let mut _total_bytes_written: usize = 0;
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
					_total_bytes_written += buf_writer(&mut code_bw, &_buff_2bit, _buff_2bit_size)?;
                    _buff_2bit_size = 0;
				}
			}
        }
        params.ln_read_buf.clear();
        buf_line_reader(&mut params.reader, &mut params.ln_read_buf).unwrap_or_else(|e| { panic!("Failed to read line in open input file: {}", e); });
    }
    if _buff_2bit_size == BUFFER_SIZE {
        _total_bytes_written += buf_writer(&mut code_bw, &_buff_2bit, _buff_2bit_size)?;
        _buff_2bit_size = 0;
    }
    if _temp_bits_used > 0 {
        _buff_2bit[_buff_2bit_size] = _buff_2bit_temp;
        _temp_bits_used = 0;
        _buff_2bit_temp = 0;
        _buff_2bit_size += 1;

    }
    if _buff_2bit_size > 0 {
        _total_bytes_written += buf_writer(&mut code_bw, &_buff_2bit, _buff_2bit_size)?;
        _buff_2bit_size = 0;
    }
    let align = _total_bytes_written % 4; //align to 32-bit word size
    if align != 0 {
        for i in 0..align as usize {
            _buff_2bit[i] = CODE_PADDING;
            _buff_2bit_size += 1;
        }
        _total_bytes_written += buf_writer(&mut code_bw, &_buff_2bit, _buff_2bit_size)?;
        _buff_2bit_size = 0;
    }

	if add_padding {
		code_bw.write_all(&PADDING_BUFF)?;
		_total_bytes_written += 4;
	}
    code_bw.flush()?;
    Ok(())
}



fn generate_executable(params: &mut Params) -> std::io::Result<()> {
	#[cfg(target_os="windows")]
	return win::exe_writer(params);
	#[cfg(target_os="linux")]
	return linux::elf_writer(params);
}

fn cleanup(params: &mut Params, tmp_path: &str) -> std::io::Result<()> {
	params.sync_files()?;
	generate_executable(params)?;
	return remove_temp_files(tmp_path,vec!["code.bin", "data.bin"]);
}

fn main() -> std::io::Result<()> {
    let main_args: Vec<String> = env::args().collect();
    let (check, tmp_path) = parse_main_args(&main_args)?;
    if !check.check_files() { return Ok(()); }
	if !check.is_req_valid() { Error::new(std::io::ErrorKind::InvalidData, "Parameters are invalid!"); }

	let mut params: Params = Params::from(check);
    let table: [u8; 256] = gen_binary_table(params.bp_val_map());
    match params.flags {
		0b00_000|0b00_001 => ls_bp_first(&mut params, &table)?,
		0b00_010|0b00_011 => ms_bp_first(&mut params, &table)?,
		0b00_100|0b00_101 => codon_ls_bp_shift(&mut params, &table)?,
		0b00_110|0b00_111 => codon_bp_fr_shift(&mut params, &table)?,
		_ => unreachable!(),
	}
	cleanup(&mut params, tmp_path.as_str())
}
