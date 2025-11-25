pub(crate) use std::io::{BufReader, BufWriter, BufRead, Write};
use std::path::PathBuf;
pub(crate) use crate::arg_parser::*;
pub struct Params {
    pub reader: BufReader<File>,
    pub writer: BufWriter<File>,
	pub code_fd: File,
	pub data_fd: Option<File>,
	pub ln_read_buf: String,
    pub compile_encode: EncodingMethod,
    #[allow(dead_code)]
	pub ifile_ext: InputFileFormat,
	pub flags: u8, // b0: 0= lsb offset shifting, 1= simple fixed bit shift, b1=codons?, b2=in-place(1) or contiguous(0) data/code for codons. b3=codon inside = code(0) or data(1)
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

	pub fn get_bp_bit_len(&self) -> u8 {
		match self.compile_encode {
			EncodingMethod::M3A0C1 | EncodingMethod::M3A1C0 | EncodingMethod::M3A0C1XOR | EncodingMethod::M3A1C0XOR => 1,
			_ => 2,
		}
	}

	pub fn flag_check(&self, mask : BitFlags) -> bool { self.flags & mask as u8 > 0 }

	pub fn sync_files(&mut self) -> std::io::Result<()> {
		self.code_fd.sync_all()?;
		if self.flag_check(BitFlags::UseCodons) && self.data_fd.is_some() {
			self.data_fd.as_mut().unwrap().sync_all()?;
		}
		Ok(())
	}

	#[inline(always)]
	fn method1_ascii_to_2bit_convert(&self, c: u8) -> u8 { (c & 0x06) >> 1 }
}

pub fn remove_temp_files(temp_dir_p: &str, temp_files: Vec<&str>) -> std::io::Result<()> {
	for item in temp_files {
		let mut fp = PathBuf::from(temp_dir_p);
		fp.push(item);
		if fp.exists() && fp.is_file() { std::fs::remove_file(fp.as_path())?; }
	}
	Ok(())
}


pub fn buf_writer<W: Write>(writer: &mut W, buf: &[u8], buff_len: usize) -> std::io::Result<usize> {
	let mut _written_bytes: usize = 0;
	while _written_bytes < buff_len {
		let n: usize = writer.write(&buf[_written_bytes..buff_len])?;
		_written_bytes += n;
	}
	Ok(_written_bytes)
}

pub fn buf_line_reader<R: BufRead>(reader: &mut R, buf: &mut String) -> std::io::Result<usize> {
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