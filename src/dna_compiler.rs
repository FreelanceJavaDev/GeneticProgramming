use std::collections::HashMap;
use crate::executable_header::CODE_PADDING;
use crate::arg_parser::{File, BitFlags};
use crate::file_handler::{BufWriter, Write, Params, buf_line_reader};

pub const BUFFER_SIZE: usize = 1024;
pub const ASCII_TABLE_SIZE: usize = 256;

const CS_FF: u8 = 0b000;
const CS_FT: u8 = 0b001;
const CS_TF: u8 = 0b010;
const CS_TT: u8 = 0b011;
//const CS_OK: u8 = 0b100;

pub trait BufferOperations {
	fn is_full(&self) -> bool;
	fn not_empty(&self) -> bool;
	fn push(&mut self, val: u8);
	fn write_buf<W: std::io::Write>(&mut self, writer: &mut W) -> std::io::Result<usize>;
}

pub trait TranslateSharedOperation {
	const TEMP_BUFF_BITS: u8 = u8::BITS as u8;
	const PADDING_BUFF: [u8; 4] = [CODE_PADDING; 4];
	const WORD_ALIGN_N_BYTES: usize = 4;
	fn ls_bp_add(&mut self, bp: u8);
	fn fixed_ls_add(&mut self, bp: u8);
	fn word_align(&mut self);
	fn flush_writer(&mut self) -> std::io::Result<()>;
	fn has_pad(&mut self, pred: bool);
	fn has_pad_end(&mut self);
}

pub mod ls {
	pub trait BP {
		fn state_fn_table(&mut self);
		fn rf_to_inside(&mut self);
		fn rf_to_outside(&mut self);
		fn flush_rf(&mut self);
		fn update_codon_buf(&mut self);

	}

	pub trait Fixed {
		fn state_fn_table(&mut self);
		fn rf_to_inside(&mut self);
		fn rf_to_outside(&mut self);
		fn flush_rf(&mut self);
		fn update_codon_buf(&mut self);
	}
}

#[repr(C)]
pub struct DNAtoBinary {
	pub buf: StoreBuf,
	pub f_out: BufWriter<File>,
	bp_map: [u8; ASCII_TABLE_SIZE],
	_total_bytes_written: usize,
	bp_bits: u8,
	bp_bits_used: u8,
	bp_shift: u8,
	add_padding: bool
}


pub struct CodonOnly {
	pub buffers: [StoreBuf; 2], // [0] is inside codon, [1] is outside codon
	pub f_writers: Vec<BufWriter<File>>, // idx 0 is inside codon, idx 1 is outside codon
	codon_map: HashMap<u8, u8>,
	bp_map: [u8; ASCII_TABLE_SIZE],
	stop_codons: [u8; 3],
	start_codon: u8,
	bytes_written: [usize; 3],
	bp_shift: u8,
	n_bits: u8,
	bp_oldest_offset: u8,
	bit_mask: u8,
	read_frame: u8,
	read_frame_used: u8,
	tmp: u8,
	tmp_bits_used: u8,
	codon_buff: u16,
	codon_buff_used: u16,
	flags: u8,
	pub codon_state: u8,
	add_padding: bool
}


impl CodonOnly {
	const INSIDE_CODON_IDX: usize = 0;
	const OUTSIDE_CODON_IDX: usize = 1;
	pub fn from(params: &Params, table: &[u8; ASCII_TABLE_SIZE]) -> Self {
		CodonOnly { codon_map: HashMap::with_capacity(64), buffers: [StoreBuf::default(), StoreBuf::default()], f_writers: Vec::with_capacity(2), bp_map: *table, stop_codons: [0; 3],
			start_codon: 0, bytes_written: [0; 3], bp_shift: params.get_bp_bit_len(), n_bits: params.get_bp_bit_len()*3, bp_oldest_offset: 0,
			bit_mask: 0, read_frame: 0, read_frame_used: 0, tmp: 0, tmp_bits_used: 0, codon_buff: 0, codon_buff_used: 0,
			flags: params.flags, codon_state: CS_FF, add_padding: params.flag_check(BitFlags::CodeSectionPadding)
		}
	}
	pub fn is_too_short(&self) -> bool { self.read_frame_used < self.n_bits }
	pub fn flag_check(&self, mask: BitFlags) -> bool { self.flags & (mask as u8) > 0 }
	pub fn dyn_init(&mut self) {
		const BP_MAP: [usize; 4] = ['A' as usize, 'C' as usize, 'G' as usize, 'T' as usize];
    	const BP_MAP_STR: [u8; 4] = [b'A', b'C', b'G', b'T'];
		let bp_shifts: [u8; 2] = [1*self.bp_shift, 2*self.bp_shift];
		let mut str_codon: HashMap<String, u8> = HashMap::with_capacity(64);
		match self.flag_check(BitFlags::FixedLeftShift) {
			true =>  {
				self.start_codon = (self.bp_map[BP_MAP[0]] << bp_shifts[1]) + (self.bp_map[BP_MAP[3]] << bp_shifts[0]) + self.bp_map[BP_MAP[2]]; // ATG
				self.stop_codons = [
					(self.bp_map[BP_MAP[3]] << bp_shifts[1]) + (self.bp_map[BP_MAP[0]] << bp_shifts[0]) + self.bp_map[BP_MAP[0]], //TAA
					(self.bp_map[BP_MAP[3]] << bp_shifts[1]) + (self.bp_map[BP_MAP[2]] << bp_shifts[0]) + self.bp_map[BP_MAP[0]], //TGA
					(self.bp_map[BP_MAP[3]] << bp_shifts[1]) + (self.bp_map[BP_MAP[0]] << bp_shifts[0]) + self.bp_map[BP_MAP[2]], //TAG
				];
				self.bp_oldest_offset = 2*self.bp_shift;
				for first in 0..BP_MAP.len() {
					for second in 0..BP_MAP.len() {
						for third in 0..BP_MAP.len() {
							let key_next: u8 = (self.bp_map[BP_MAP[first]] << bp_shifts[1]) | (self.bp_map[BP_MAP[second]] << bp_shifts[0]) | self.bp_map[BP_MAP[third]];
							str_codon.insert(String::from_utf8(vec![BP_MAP_STR[first],BP_MAP_STR[second],BP_MAP_STR[third]]).unwrap(),key_next);
						}
					}
				}
				self.bit_mask = match self.bp_shift {
					1 => 0b00000_011,
					2 => 0b00_001111,
					_ => unreachable!(),
				};
			},
			false =>  {
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
							str_codon.insert(String::from_utf8(vec![BP_MAP_STR[first],BP_MAP_STR[second],BP_MAP_STR[third]]).unwrap(),key_next);
						}
					}
				}
				self.bit_mask = match self.bp_shift {
					1 => 0b01,
					2 => 0b11,
					_ => unreachable!(),
				};
			},
		}
		self.gen_codon_map_raw(str_codon);

	}

	pub fn update_state(&mut self) {
		self.codon_state = (self.codon_state << 1) & CS_TF;
		self.codon_state |= match self.codon_state {
			CS_TF => self.state_at_stop_codon(),
			CS_FF => self.state_in_start_codon(),
			_ => unreachable!(),
		};
	}
	pub fn is_inside_codon(&self) -> bool { self.codon_state == CS_FT || self.codon_state == CS_TT }

	pub fn reset_tmp_buf(&mut self) {
		self.tmp = 0;
		self.tmp_bits_used = 0;
	}

	pub fn flush_buff(&mut self, id: usize) {
		let n_bytes = self.buffers[id].write_buf(&mut self.f_writers[id]).expect("Failed to write buffered data to file");
		self.bytes_written[id] += n_bytes;
		self.bytes_written[2] += n_bytes;

	}

	pub fn flush_codon_buf_not_empty(&mut self) {
		while self.codon_buff_used >= 8 {
			self.buffers[Self::INSIDE_CODON_IDX].push( self.codon_buff as u8);
			self.codon_buff >>= 8;
			self.codon_buff_used -= 8;
			if self.buffers[Self::INSIDE_CODON_IDX].is_full() { self.flush_buff(Self::INSIDE_CODON_IDX); }
		}
		if self.codon_buff_used > 0 {
			self.buffers[Self::INSIDE_CODON_IDX].push( self.codon_buff as u8);
			self.codon_buff >>= self.codon_buff_used;
			self.codon_buff_used = 0;
		}
	}

	pub fn flush_tmp_not_empty(&mut self) {
		if self.tmp_bits_used > 0 {
			self.buffers[Self::OUTSIDE_CODON_IDX].push(self.tmp);
			self.reset_tmp_buf();
		}
	}

	pub fn read_frame_not_empty(&self) -> bool { self.read_frame_used > 0 }

	pub fn flush_buffers_not_empty(&mut self) {
		if self.buffers[Self::INSIDE_CODON_IDX].not_empty() { self.flush_buff(Self::INSIDE_CODON_IDX); }
		if self.buffers[Self::OUTSIDE_CODON_IDX].not_empty() { self.flush_buff(Self::OUTSIDE_CODON_IDX);}
	}
	pub fn reset_read_frame(&mut self) {
		self.read_frame = 0;
		self.read_frame_used = 0;
	}

	pub fn init_file_writers(&mut self, params: &Params) {
		let f_code: File = params.code_fd.try_clone().expect("Failed clone code temp file");
		let f_data: File = if params.data_fd.is_some() {
			params.data_fd.as_ref().unwrap().try_clone().expect("Failed to clone data temp file")
		}
		else { params.code_fd.try_clone().expect("Failed copy code temp file as data temp file") };
		match self.is_codon_data() {
			true => { self.f_writers.push(BufWriter::new(f_data)); self.f_writers.push(BufWriter::new(f_code)); },
			false => { self.f_writers.push(BufWriter::new(f_code)); self.f_writers.push(BufWriter::new(f_data)); },
		}
	}

}


#[repr(C)]
pub struct StoreBuf {
	pub buf: [u8; BUFFER_SIZE],
	pub buf_used: usize
}

impl DNAtoBinary {
	pub fn from(params: &Params, table: &[u8; ASCII_TABLE_SIZE]) -> Self {
		let cpy: File = params.code_fd.try_clone().expect("Failed clone output file");
		DNAtoBinary { buf: StoreBuf::default(), f_out: BufWriter::new(cpy), bp_map: *table,
			bp_bits: 0, bp_bits_used: 0, _total_bytes_written: 0,
			bp_shift: params.get_bp_bit_len(),
			add_padding: params.flag_check(BitFlags::CodeSectionPadding)
		}
	}

	pub fn reset_bp_bits(&mut self) {
		self.bp_bits = 0;
		self.bp_bits_used = 0;
	}

	pub fn flush_if_full_buf(&mut self) {
		if self.buf.is_full() {
			self._total_bytes_written += self.buf.write_buf(&mut self.f_out).expect("Failed to write full buffer to file");
		}
	}

	pub fn flush_partial(&mut self) {
		if self.buf.not_empty() {
			self._total_bytes_written += self.buf.write_buf(&mut self.f_out).expect("Failed to write buffered data to file");
		}
	}

	pub fn flush_bp_bits(&mut self) {
		if self.bp_bits_used > 0 { self.push_bp_bits(); }
	}

	pub fn push_bp_bits(&mut self) {
		self.buf.push(self.bp_bits);
		self.reset_bp_bits();

	}
}

impl TranslateSharedOperation for DNAtoBinary {
	fn ls_bp_add(&mut self, bp: u8) {
		self.bp_bits += (self.bp_map[bp as usize]) << self.bp_bits_used;
		self.bp_bits_used += self.bp_shift;
		if self.bp_bits_used == Self::TEMP_BUFF_BITS {
			self.push_bp_bits();
			self.flush_if_full_buf();
		}
	}

	fn fixed_ls_add(&mut self, bp: u8) {
		self.bp_bits <<= self.bp_shift;
		self.bp_bits += self.bp_map[bp as usize];
		self.bp_bits_used += self.bp_shift;
		if self.bp_bits_used == Self::TEMP_BUFF_BITS {
			self.push_bp_bits();
			self.flush_if_full_buf();
		}
	}

	fn word_align(&mut self)  {
		let align: usize = self._total_bytes_written % Self::WORD_ALIGN_N_BYTES; //align to 32-bit word size
		if align != 0 {
			for _ in 0..align as usize {
				self.buf.push(CODE_PADDING);
			}
			self._total_bytes_written += self.buf.write_buf(&mut self.f_out).expect("BufWriter failed to write to file");
		}
	}
	fn flush_writer(&mut self) -> std::io::Result<()> { self.f_out.flush() }
	fn has_pad(&mut self, pred: bool) {
		if self.add_padding && pred {
			self.f_out.write_all(&Self::PADDING_BUFF).expect("Failed write padding bytes output file");
			self._total_bytes_written += Self::WORD_ALIGN_N_BYTES;
		}
	}

	fn has_pad_end(&mut self) {
		if self.add_padding && self._total_bytes_written > 0 {
			self.f_out.write_all(&Self::PADDING_BUFF).expect("Failed write padding bytes output file");
			self._total_bytes_written += Self::WORD_ALIGN_N_BYTES;
		}
	}

}

impl TranslateSharedOperation for CodonOnly {
	fn ls_bp_add(&mut self, bp: u8) {
		self.read_frame += self.bp_map[bp as usize] << self.read_frame_used;
		self.read_frame_used += self.bp_shift;
	}

	fn fixed_ls_add(&mut self, bp: u8) {
		self.read_frame <<= self.bp_shift;
		self.read_frame += self.bp_map[bp as usize];
		self.read_frame_used += self.bp_shift;
	}

	fn word_align(&mut self) {
		let idx: usize = self.get_code_idx();
		let align: usize = self.bytes_written[idx] % Self::WORD_ALIGN_N_BYTES; //align to 32-bit word size
		if align != 0 {
			for _ in 0..align as usize {
				self.buffers[idx].push(CODE_PADDING);
			}
			self.flush_buff(idx);
		}
	}


	fn has_pad(&mut self, pred: bool) {
		if self.add_padding && pred {
			let idx: usize = self.get_code_idx();
			self.f_writers[idx].write_all(&Self::PADDING_BUFF).expect("Failed write padding bytes output to file");
			self.bytes_written[idx] += Self::WORD_ALIGN_N_BYTES;

			self.bytes_written[2] += Self::WORD_ALIGN_N_BYTES;
		}
	}

	fn has_pad_end(&mut self) {
		self.has_pad(self.bytes_written[self.get_code_idx()] > 0);
	}

	fn flush_writer(&mut self) -> std::io::Result<()> {
		self.f_writers[Self::OUTSIDE_CODON_IDX].flush()?;
		self.f_writers[Self::INSIDE_CODON_IDX].flush()
	}
}

impl BufferOperations for StoreBuf {
	fn push(&mut self, val: u8) {
		self.buf[self.buf_used] = val;
		self.buf_used += 1;
	}

	fn is_full(&self) -> bool { self.buf_used == BUFFER_SIZE }

	fn not_empty(&self) -> bool { self.buf_used > 0 }

	fn write_buf<W: std::io::Write>(&mut self, writer: &mut W) -> std::io::Result<usize> {
		let mut _written_bytes: usize = 0;
		while _written_bytes < self.buf_used {
			_written_bytes += writer.write(&self.buf[_written_bytes..self.buf_used])?;
		}
		self.buf_used = 0;
		Ok(_written_bytes)
	}
}

//Left shift Base Pair before adding to read frame ls_bp
impl ls::BP for CodonOnly {
	fn state_fn_table(&mut self) {
		self.update_state();
		match self.codon_state {
			CS_FF => <CodonOnly as ls::BP>::rf_to_outside(self), //outside codon
			CS_FT => <CodonOnly as ls::BP>::rf_to_inside(self), //Transition from outside to inside (Start codon)
			CS_TF => <CodonOnly as ls::BP>::rf_to_inside(self), //Transition from inside to outside (Stop codon reached)
			CS_TT => <CodonOnly as ls::BP>::rf_to_inside(self), //Inside codon
			_ => unreachable!()
		}
	}

	fn rf_to_inside(&mut self) {
		<CodonOnly as ls::BP>::update_codon_buf(self);
		if self.codon_buff_used >= 8 {
			self.buffers[Self::INSIDE_CODON_IDX].push( self.codon_buff as u8);
			self.codon_buff >>= 8;
			self.codon_buff_used -= 8;
			if self.buffers[Self::INSIDE_CODON_IDX].is_full() { self.flush_buff(Self::INSIDE_CODON_IDX); }
		}
	}
	fn flush_rf(&mut self) {
		if self.is_too_short() || !self.is_inside_codon() {
			let tmp_remain: u8 = Self::TEMP_BUFF_BITS - self.tmp_bits_used;
			let overflow: bool = tmp_remain < self.read_frame_used;

			self.tmp += self.read_frame << self.tmp_bits_used;
			self.tmp_bits_used += tmp_remain;
			self.buffers[Self::OUTSIDE_CODON_IDX].push(self.tmp);
			self.reset_tmp_buf();
			if self.buffers[Self::OUTSIDE_CODON_IDX].is_full() { self.flush_buff(Self::OUTSIDE_CODON_IDX); }

			if overflow {
				self.read_frame >>= tmp_remain;
				self.read_frame_used -= tmp_remain;
				self.tmp += self.read_frame;
				self.tmp_bits_used += self.read_frame_used;
			}
			self.reset_read_frame();
		}
		else { <CodonOnly as ls::BP>::update_codon_buf(self); }
	}

	fn rf_to_outside(&mut self) {
		self.tmp += (self.read_frame & self.bit_mask) << self.tmp_bits_used;
		self.tmp_bits_used += self.bp_shift;
		self.read_frame >>= self.bp_shift;
		self.read_frame_used -= self.bp_shift;
		if self.tmp_bits_used == Self::TEMP_BUFF_BITS {
			self.buffers[Self::OUTSIDE_CODON_IDX].push(self.tmp);
			self.reset_tmp_buf();
			if self.buffers[Self::OUTSIDE_CODON_IDX].is_full() { self.flush_buff(Self::OUTSIDE_CODON_IDX); }
		}
	}

	fn update_codon_buf(&mut self) {
		self.codon_buff += (self.read_frame as u16) << self.codon_buff_used;
		self.codon_buff_used += self.read_frame_used as u16;
		self.reset_read_frame();
	}

}

//Fixed left shift then add
impl ls::Fixed for CodonOnly {
	fn state_fn_table(&mut self) {
		self.update_state();
		match self.codon_state {
			CS_FF => <CodonOnly as ls::Fixed>::rf_to_outside(self), //outside codon
			CS_FT => <CodonOnly as ls::Fixed>::rf_to_inside(self), //Transition from outside to inside (Start codon)
			CS_TF => <CodonOnly as ls::Fixed>::rf_to_inside(self), //Transition from inside to outside (Stop codon reached)
			CS_TT => <CodonOnly as ls::Fixed>::rf_to_inside(self), //Inside codon
			_ => unreachable!()
		}
	}

	fn rf_to_inside(&mut self) {
		<CodonOnly as ls::Fixed>::update_codon_buf(self);
		if self.codon_buff_used >= 8 {
			self.buffers[Self::INSIDE_CODON_IDX].push( self.codon_buff as u8);
			self.codon_buff >>= 8;
			self.codon_buff_used -= 8;
			if self.buffers[Self::INSIDE_CODON_IDX].is_full() { self.flush_buff(Self::INSIDE_CODON_IDX); }
		}
	}

	fn rf_to_outside(&mut self) {
		self.tmp <<= self.bp_shift;
		self.tmp += self.read_frame >> self.bp_oldest_offset;
		self.tmp_bits_used += self.bp_shift;
		self.read_frame &= self.bit_mask;
		self.read_frame_used -= self.bp_shift;
		if self.tmp_bits_used == Self::TEMP_BUFF_BITS {
			self.buffers[Self::OUTSIDE_CODON_IDX].push(self.tmp);
			self.reset_tmp_buf();
			if self.buffers[Self::OUTSIDE_CODON_IDX].is_full() { self.flush_buff(Self::OUTSIDE_CODON_IDX); }
		}
	}

	fn flush_rf(&mut self) {
		if self.is_too_short() || !self.is_inside_codon() {
			let tmp_remain = Self::TEMP_BUFF_BITS - self.tmp_bits_used;
			let overflow = tmp_remain < self.read_frame_used;
			let mut bp_mask: u8 = match self.bp_shift {
				1 => match self.read_frame_used {
						3 => 0b00000_111,
						2 => 0b00000_011,
						1 => 0b00000_001,
						_ => unreachable!(),
					},
				2 => match self.read_frame_used {
					6 => 0b00_111111,
					4 => 0b00_001111,
					2 => 0b00_000011,
					_ => unreachable!(),
				},
				_ => unreachable!(),
			};
			let mut rd_frame: u8 = self.read_frame;
			let mut offset = if overflow { tmp_remain } else { self.read_frame_used};
			while offset > 0 {
				rd_frame &= bp_mask;
				offset -= self.bp_shift;
				self.tmp += (rd_frame >> offset) << self.tmp_bits_used;
				self.tmp_bits_used += self.bp_shift;
				bp_mask >>= self.bp_shift;
			}
			self.buffers[Self::OUTSIDE_CODON_IDX].push(self.tmp);
			self.reset_tmp_buf();
			if self.buffers[Self::OUTSIDE_CODON_IDX].is_full() { self.flush_buff(Self::OUTSIDE_CODON_IDX); }
			if overflow {
				self.read_frame &= self.bit_mask >> tmp_remain;
				self.read_frame_used -= tmp_remain;
				self.tmp += self.read_frame;
				self.tmp_bits_used += self.read_frame_used;
			}
			self.reset_read_frame();

		}
		else { <CodonOnly as ls::Fixed>::update_codon_buf(self); }
	}

	fn update_codon_buf(&mut self) {
		self.codon_buff <<= self.read_frame_used;
		self.codon_buff += self.read_frame as u16;
		self.codon_buff_used += self.read_frame_used as u16;
		self.reset_read_frame();
	}

}

pub fn codon_ls_bp_shift(params: &mut Params, table: &[u8; ASCII_TABLE_SIZE]) -> std::io::Result<()> {
	let mut codon: CodonOnly = CodonOnly::from(params, table);
	codon.dyn_init();
	codon.init_file_writers(params);
	buf_line_reader(&mut params.reader, &mut params.ln_read_buf)?;
	codon.has_pad(params.read_buf_not_empty());
	while params.read_buf_not_empty() {
		for bp in params.ln_read_buf.bytes() {
			codon.ls_bp_add(bp);
			if codon.is_too_short() { continue; }
			else { ls::BP::state_fn_table(&mut codon); }
		}
		params.ln_read_buf.clear();
		buf_line_reader(&mut params.reader, &mut params.ln_read_buf)?;
	}

	if codon.read_frame_not_empty() { ls::BP::flush_rf(&mut codon); }

	codon.flush_tmp_not_empty();
	codon.flush_codon_buf_not_empty();
	codon.flush_buffers_not_empty();
	codon.word_align();
	codon.has_pad_end();
	codon.flush_writer()
}

pub fn codon_bp_fr_shift(params: &mut Params, table: &[u8; ASCII_TABLE_SIZE]) -> std::io::Result<()> {
    let mut codon: CodonOnly = CodonOnly::from(params, table);
	codon.dyn_init();
	codon.init_file_writers(params);
	buf_line_reader(&mut params.reader, &mut params.ln_read_buf)?;
    codon.has_pad(params.read_buf_not_empty());

	while params.read_buf_not_empty() {
		for bp in params.ln_read_buf.bytes() {
			codon.fixed_ls_add(bp);
			if codon.is_too_short() { continue; }
			else { ls::Fixed::state_fn_table(&mut codon); }
		}
		params.ln_read_buf.clear();
		buf_line_reader(&mut params.reader, &mut params.ln_read_buf)?;
	}
	if codon.read_frame_not_empty() { ls::Fixed::flush_rf(&mut codon); }

	codon.flush_tmp_not_empty();
	codon.flush_codon_buf_not_empty();
	codon.flush_buffers_not_empty();
	codon.word_align();
	codon.has_pad_end();
	codon.flush_writer()
}

pub fn fixed_left_shift_first(params: &mut Params, table: &[u8; ASCII_TABLE_SIZE]) -> std::io::Result<()> {
	let mut translator: DNAtoBinary = DNAtoBinary::from(params, table);
	buf_line_reader(&mut params.reader, &mut params.ln_read_buf).expect("Failed to read line in open input file");
	translator.has_pad(params.read_buf_not_empty());
	while params.read_buf_not_empty() {
		for bp in params.ln_read_buf.bytes() {
			translator.fixed_ls_add(bp);
		}
		params.ln_read_buf.clear();
		buf_line_reader(&mut params.reader, &mut params.ln_read_buf).expect("Failed to read line in open input file");
	}

	translator.flush_if_full_buf();
	translator.flush_bp_bits();
	translator.flush_partial();
	translator.word_align();
	translator.has_pad_end();
	translator.flush_writer()
}

pub fn ls_bp_first(params: &mut Params, table: &[u8; ASCII_TABLE_SIZE]) -> std::io::Result<()> {
	let mut translator: DNAtoBinary = DNAtoBinary::from(params, table);
	buf_line_reader(&mut params.reader, &mut params.ln_read_buf).expect("Failed to read line in open input file");
	translator.has_pad(params.read_buf_not_empty());

	while params.read_buf_not_empty() {
		for bp in params.ln_read_buf.bytes() {
			translator.ls_bp_add(bp);
		}
		params.ln_read_buf.clear();
		buf_line_reader(&mut params.reader, &mut params.ln_read_buf).expect("Failed to read line in open input file");
	}
	translator.flush_if_full_buf();
	translator.flush_bp_bits();
	translator.flush_partial();
	translator.word_align();
	translator.has_pad_end();
	translator.flush_writer()
}

pub fn gen_binary_table(bp_map_val: [u8; 4]) -> [u8; ASCII_TABLE_SIZE] {
	const BP_MAP: [usize; 8] = ['A' as usize, 'a' as usize, 'C' as usize, 'c' as usize, 'G' as usize, 'g' as usize, 'T' as usize, 't' as usize,];
	let mut ret: [u8; ASCII_TABLE_SIZE] = [u8::MAX; ASCII_TABLE_SIZE];
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

impl Default for StoreBuf {
	fn default() -> Self { StoreBuf { buf: [0; BUFFER_SIZE], buf_used: 0 } }
}

//Private method implementation
impl CodonOnly {
	fn is_codon_data(&self) -> bool { self.flags & (BitFlags::CodonData as u8) > 0 }

	fn get_code_idx(&self) -> usize {
		match self.is_codon_data() {
			true => Self::OUTSIDE_CODON_IDX,
			false => Self::INSIDE_CODON_IDX,
		}
	}

	fn state_at_stop_codon(&self) -> u8 {
		let ck: u8 = self.codon_map[&self.read_frame];
		!(ck == self.stop_codons[0] || ck == self.stop_codons[1] || ck == self.stop_codons[2]) as u8
	}

	fn state_in_start_codon(&self) -> u8 {
		let ck: u8 = self.codon_map[&self.read_frame];
		(ck == self.start_codon) as u8
	}

	fn gen_codon_map_raw(&mut self, str_map: HashMap<String, u8>) {
		for item in str_map {
			if self.codon_map.contains_key(&item.1) { continue; }
			match item.0.as_str() {
				_ => { self.codon_map.insert(item.1, item.1); },
			}
		}
	}
}


