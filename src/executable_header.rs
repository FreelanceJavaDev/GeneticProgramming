
pub trait Serialize {
	fn serialize(&self) -> Vec<u8>;
}

pub const SECTION_ALIGNMENT: u32 = 0x1000;
#[allow(dead_code)]
pub const FILE_ALIGNMENT: u32 = 0x200;
pub const HEADER_PADDING_BYTE: u8 = 0x00;
pub const WORD_ALIGN_BYTES: u32 = 4; ///32-bit alignment = 4 bytes

#[cfg(target_os="windows")]
pub const CODE_PADDING: u8 = 0xCC; // INT3 opcode
#[cfg(not(target_os="windows"))]
pub const CODE_PADDING: u8 = 0x90; // NO OPERATION opcode
#[cfg(target_os="windows")]
pub mod win  {
	use crate::executable_header::Serialize;
	use crate::file_handler::*;
	use std::io::{Read, Seek};

	pub fn exe_writer(params: &mut Params)  -> std::io::Result<()> {
		const ALIGN_FILE: usize = super::FILE_ALIGNMENT as usize;
		let mut f_sizes: Vec<u32> = Vec::new();
		let code_size: u64 = params.code_fd.metadata()?.len();
		let mut temp_reader: BufReader<&File> = BufReader::new(&params.code_fd);
		let entry_offset: u32 = if params.flag_check(BitFlags::CodeSectionPadding) { super::WORD_ALIGN_BYTES } else { 0 };
		temp_reader.rewind()?;
		f_sizes.push(code_size as u32);
		let f_data_valid = params.flag_check(BitFlags::UseCodons) && params.data_fd.is_some();
		if f_data_valid {
			f_sizes.push(params.data_fd.as_ref().unwrap().metadata()?.len() as u32);
		}

		let header = generate_executable_header(f_sizes.clone(), f_sizes.len() as u16, entry_offset);
		let header_buff = header.serialize();
		let mut _bytes_written: usize = buf_writer(&mut params.writer, &header_buff.as_slice(), header_buff.len())?;
		let buff_align: usize = ALIGN_FILE - (_bytes_written % ALIGN_FILE);
		let mut write_buff: [u8; ALIGN_FILE] = [super::HEADER_PADDING_BYTE; ALIGN_FILE];
		if buff_align < ALIGN_FILE {
			_bytes_written += buf_writer(&mut params.writer, &write_buff, buff_align)?;
		}
		params.writer.flush()?;
		for _full_reads in 0..f_sizes[0]/super::FILE_ALIGNMENT {
			temp_reader.read_exact(&mut write_buff).unwrap_or_else(|e| { panic!("code.bin file reader error: {}", e); });
			_bytes_written += buf_writer(&mut params.writer, &write_buff,ALIGN_FILE)?;
		}
		let mut end_vec: Vec<u8> = Vec::new();
		let code_data_end = temp_reader.read_to_end(&mut end_vec)?;
		_bytes_written += buf_writer(&mut params.writer, end_vec.as_slice(), end_vec.len())?;
		let align = ALIGN_FILE - (code_data_end % ALIGN_FILE);

		if align > 0 {
			write_buff = [super::CODE_PADDING; ALIGN_FILE];
			_bytes_written += buf_writer(&mut params.writer, &write_buff, align)?;
		}
		params.writer.flush()?;
		if f_data_valid {
			let mut data_reader: BufReader<&File> = BufReader::new(params.data_fd.as_mut().unwrap());
			data_reader.rewind()?;
			write_buff = [super::CODE_PADDING; ALIGN_FILE];
			for _full_buf in 0..f_sizes[1]/super::FILE_ALIGNMENT {
				data_reader.read_exact(&mut write_buff).unwrap_or_else(|e| { panic!("data.bin file reader error: {}", e); });
				_bytes_written += buf_writer(&mut params.writer, &write_buff,ALIGN_FILE)?;
			}
			params.writer.flush()?;
			let mut end_vec: Vec<u8> = Vec::new();
			let data_end = data_reader.read_to_end(&mut end_vec)?;
			_bytes_written += buf_writer(&mut params.writer, end_vec.as_slice(), end_vec.len())?;
			let align_data = ALIGN_FILE - (data_end % ALIGN_FILE);
			if align_data > 0 {
				write_buff = [super::CODE_PADDING; ALIGN_FILE];
				_bytes_written += buf_writer(&mut params.writer, &write_buff, align_data)?;
			}
			params.writer.flush()?;
		}

		Ok(())
	}

	const fn get_header_size(num_sections: usize) -> usize {
		use std::mem::size_of;
		size_of::<InitFileHeader>() + size_of::<ImageNTHeaders32>() + (size_of::<ImageSectionHeaderEnt>() * num_sections)
	}

	fn calc_alignment(len: u32, align_val: u32) -> u32 {
		(match len % align_val {
			0 => len/align_val,
			_ => (len/align_val) + 1,
		})*align_val
	}

	fn generate_executable_header(sec_sizes: Vec<u32>, num_sections: u16, entry_byte_offset: u32) -> HeaderPE {
		if num_sections == 1 { return generate_executable_code_only(sec_sizes[0], entry_byte_offset); }
		else { return generate_executable_code_data(sec_sizes, entry_byte_offset); }
	}


	fn generate_executable_code_only(code_len: u32, entry_byte_offset: u32) -> HeaderPE {
		let data_dir_list: [ImageDataDirectory; 16] = [ImageDataDirectory::default(); 16];
		let header_data_len: u32 = get_header_size(1) as u32;
		let header_file_align: u32 = calc_alignment(header_data_len, super::FILE_ALIGNMENT);
		let sec_file_align: u32 = calc_alignment(code_len, super::FILE_ALIGNMENT);
		let header_va_offset: u32 = calc_alignment(header_data_len,super::SECTION_ALIGNMENT);
		let va_data_addr_offset: u32 = header_va_offset + calc_alignment(code_len, super::SECTION_ALIGNMENT);
		let opt_file_header: OptionalImageFileHeader32 = OptionalImageFileHeader32::with_code_only(sec_file_align, header_va_offset+entry_byte_offset, header_va_offset, va_data_addr_offset, data_dir_list);
		HeaderPE::with_defaults(1, opt_file_header, vec![ImageSectionHeaderEnt::with_defaults(*b".text\0\0\0", code_len, header_va_offset, sec_file_align, header_file_align)])
	}

	fn generate_executable_code_data(sec_sizes: Vec<u32>, entry_byte_offset: u32) -> HeaderPE {
		let data_dir_list: [ImageDataDirectory; 16] = [ImageDataDirectory::default(); 16];

		let header_data_len: u32 = get_header_size(2) as u32;

		let header_file_align: u32 = calc_alignment(header_data_len, super::FILE_ALIGNMENT);
		let va_code_offset: u32 = calc_alignment(header_data_len,super::SECTION_ALIGNMENT);
		let code_size_align: u32 = calc_alignment(sec_sizes[0], super::FILE_ALIGNMENT);
		let data_size_align: u32 = calc_alignment(sec_sizes[1], super::FILE_ALIGNMENT);
		let va_data_offset: u32 = va_code_offset + calc_alignment(sec_sizes[0], super::SECTION_ALIGNMENT);
		let code_file_align: u32 = header_file_align + calc_alignment(sec_sizes[0], super::FILE_ALIGNMENT);
		let img_sz: u32 = va_data_offset + calc_alignment(sec_sizes[1],super::SECTION_ALIGNMENT);

		let opt_file_header: OptionalImageFileHeader32 = OptionalImageFileHeader32::with_code_and_data(code_size_align, data_size_align, va_code_offset+entry_byte_offset, va_code_offset, va_data_offset, img_sz, data_dir_list);
		let sect_header_list: Vec<ImageSectionHeaderEnt> = vec![
			ImageSectionHeaderEnt::with_defaults(*b".text\0\0\0", sec_sizes[0], va_code_offset, code_size_align, header_file_align),
			ImageSectionHeaderEnt::with_perms(*b".data\0\0\0", sec_sizes[1], va_data_offset, data_size_align, code_file_align, 0xC0000040)
		];
		HeaderPE::with_defaults(2, opt_file_header, sect_header_list)
}

	#[derive(Clone)]
	#[repr(C)]
	pub struct HeaderPE {
		pub init_header: InitFileHeader,
		pub nt_header: ImageNTHeaders32,
		pub section_header: Vec<ImageSectionHeaderEnt> //0x28 bytes per entry
	} //Size is 0x1E8 + 0x28*header count, Min is 0x210, min alignment is 0x400

	#[derive(Copy, Clone)]
	#[repr(C)]
	pub struct InitFileHeader {
		pub h1: [u8; 0x80],
		pub pad: [u8; 0x70 ]

	} //size is 0xF0

	#[derive(Copy, Clone)]
	#[repr(C)]
	pub struct ImageNTHeaders32 {
		pub signature: u32, //Always 0x4550
		pub file_header: ImageFileHeader,
		pub opt_file_header: OptionalImageFileHeader32
	} //Size is 0xF8
	#[derive(Copy, Clone)]
	#[repr(C)]
	pub struct ImageFileHeader {
		pub machine: u16, //0x014c = x86, 0x8664 = amd64/x86_64, 0x0000 = unknown/any
		pub num_sections: u16, //Num entries in section header
		pub td_stamp: u32, // Don't care about time-date stamp
		pub sym_tbl_ptr: u32, /// Must be 0
		pub num_symbols: u32, /// Must be 0
		pub opt_header_size: u16, //Either e0 for 32-bit 0r f0 for 64-bit
		pub characteristics: u16 //0x0102 for 32-bit, 0x22 for 64-bit
	} //size is 0x14

	#[derive(Copy, Clone)]
	#[repr(C)]
	pub struct OptionalImageFileHeader32 {
		pub magic: u16, // 0x10b for 32-bit 0x20b for 64-bit
		pub linker_ver: [u8; 2], // [major,minor]
		pub code_size: u32,
		pub idata_size: u32, //might be 0?
		pub uninitialized_data_size: u32, //0
		pub	entry_point_addr: u32, // start of code
		pub base_of_code: u32, //base address of code section
		pub base_of_data: u32, // Might be 0?
		pub image_base: u32, //0x00400000 = for 32-bit
		pub section_align: u32, //Section alignment 0x1000
		pub file_align: u32, //raw data alignment, default: 0x200
		pub sys_ver: [u8; 12], // Combined os, image and subsystem versions : Set as le bytes: [0x06,0x00, 0x00,0x00, 0x00,0x00, 0x00,0x00, 0x06,0x00, 0x00,0x00]
		pub w32_version_value: u32, // Always 0
		pub image_size: u32, //size of image in bytes rounded up to section alignment
		pub headers_size: u32, //Default 0x400
		pub ck_sum: u32,  //No checksum
		pub subsystem: u16, // 0 = IMAGE_SUBSYSTEM_UNKNOWN, 2= Windows GUI application, 3= Console app,
		pub dll_characteristics: u16, //Default 0x8740; 0x40 = DLL_CAN_MOVE, 0x100 = DLL_NX_COMPATIBLE, 0x200 = DLL_NO_ISOLATION, 0x400 = DLL_NO_SEH, 0x8000 = DLL_TERMINAL_SERVER_AWARE
		pub stack_reserve_size: u32, // Default 0x100000
		pub stack_commit_size: u32, // Default 0x1000
		pub heap_reserve_size: u32, // Default 0x100000
		pub heap_commit_size: u32, // Default 0x1000
		pub loader_flags: u32, //Always 0
		pub rva_num_and_sizes: u32, //always 16
		pub data_directory: [ImageDataDirectory; 16]
	} //size is 0xE0

	#[derive(Copy, Clone)]
	#[repr(C)]
	pub struct ImageSectionHeaderEnt {
		pub name: [u8; 8],
		pub sec_size: u32, //Size of actual data
		pub virtual_address: u32, //virtual address
		pub size_raw_data: u32, // size of raw data aligned to file
		pub ptr_raw_data: u32, // physical address (aligned to file)
		pub legacy: [u8; 12], // ptr_relocs, ptr_line_nums, num_relocs and num_line_nums are always 0
		pub characteristics: u32 //default is: 60000020  (code section: read & execute)
	} //size is 0x28

	#[derive(Copy, Clone)]
	#[repr(C)]
	pub struct ImageDataDirectory {
		pub virtual_address: u32,
		pub size: u32
	}

	impl HeaderPE {
		pub fn with_defaults(sec_num: u16, opt_header: OptionalImageFileHeader32, sect_list: Vec<ImageSectionHeaderEnt>) -> Self {
			HeaderPE { init_header: InitFileHeader::default(),
				nt_header: ImageNTHeaders32::with_defaults(sec_num, opt_header), section_header: sect_list }
		}
	}

	impl ImageNTHeaders32 {
		pub fn with_defaults(sec_num: u16, opt_hdr: OptionalImageFileHeader32) -> ImageNTHeaders32 {
			ImageNTHeaders32 { signature: 0x4550, file_header: ImageFileHeader::with_defaults(sec_num),  opt_file_header: opt_hdr }
		}
	}

	impl ImageFileHeader {
		pub fn with_defaults(sec_num: u16) ->  ImageFileHeader {
			ImageFileHeader { machine: 0x014c, num_sections: sec_num, td_stamp: 0, sym_tbl_ptr: 0, num_symbols: 0, opt_header_size: 0x00e0, characteristics: 0x0102 }
		}
		#[allow(dead_code)]
		pub fn with_params(mach: u16, sec_num: u16, opt_head_size: u16, flags: u16) ->  ImageFileHeader {
			ImageFileHeader { machine: mach, num_sections: sec_num, td_stamp: 0, sym_tbl_ptr: 0, num_symbols: 0, opt_header_size: opt_head_size, characteristics: flags }
		}
	}

	impl OptionalImageFileHeader32 {
		pub fn with_code_only(code_sz: u32, start_addr: u32, code_base: u32, image_sz: u32, data_dir: [ImageDataDirectory; 16]) -> OptionalImageFileHeader32 {
			OptionalImageFileHeader32 { magic: 0x010B, linker_ver: [0x00, 0x00], code_size: code_sz, idata_size: 0, uninitialized_data_size: 0,
				entry_point_addr: start_addr, base_of_code: code_base, base_of_data: ((((code_sz + code_base) / 0x1000) + 1)*0x1000), image_base: 0x00400000,
				section_align: 0x1000, file_align: 0x200, sys_ver: [0x06,0x00, 0x00,0x00, 0x00,0x00 ,0x00,0x00, 0x06,0x00, 0x00,0x00], w32_version_value: 0,
				image_size: image_sz, headers_size: 0x400, ck_sum: 0, subsystem: 3, dll_characteristics: 0x8740, stack_reserve_size: 0x100000, stack_commit_size: 0x1000, heap_reserve_size: 0x100000,
				heap_commit_size: 0x1000, loader_flags: 0, rva_num_and_sizes: 0x10,
				data_directory: data_dir
			}
		}

		pub fn with_code_and_data(code_sz: u32, idata_sz: u32, start_addr: u32, code_base: u32, idata_base: u32, image_size: u32, data_dir: [ImageDataDirectory; 16]) -> OptionalImageFileHeader32 {
			OptionalImageFileHeader32 { magic: 0x010B, linker_ver: [0x0E, 0x2B], code_size: code_sz, idata_size: idata_sz, uninitialized_data_size: 0,
				entry_point_addr: start_addr, base_of_code: code_base, base_of_data: idata_base, image_base: 0x00400000,
				section_align: 0x1000, file_align: 0x200, sys_ver: [0x06,0x00, 0x00,0x00, 0x00,0x00 ,0x00,0x00, 0x06,0x00, 0x00,0x00], w32_version_value: 0,
				image_size, headers_size: 0x400, ck_sum: 0, subsystem: 3, dll_characteristics: 0x8140, stack_reserve_size: 0x100000, stack_commit_size: 0x1000, heap_reserve_size: 0x100000,
				heap_commit_size: 0x1000, loader_flags: 0, rva_num_and_sizes: 0x10,
				data_directory: data_dir
			}
		}
		#[allow(dead_code)]
		pub fn with_params(bit_arch: u16, code_sz: u32, idata_sz: u32, start_addr: u32, code_base: u32, idata_base: u32, image_sz: u32, data_dir: [ImageDataDirectory; 16]) -> OptionalImageFileHeader32 {
			OptionalImageFileHeader32 { magic: bit_arch, linker_ver: [0x0E, 0x2B], code_size: code_sz, idata_size: idata_sz, uninitialized_data_size: 0,
				entry_point_addr: start_addr, base_of_code: code_base, base_of_data: idata_base, image_base: 0x00400000,
				section_align: 0x1000, file_align: 0x200, sys_ver: [0x06,0x00, 0x00,0x00, 0x00,0x00 ,0x00,0x00, 0x06,0x00, 0x00,0x00], w32_version_value: 0,
				image_size: image_sz, headers_size: 0x400, ck_sum: 0, subsystem: 3, dll_characteristics: 0x8140, stack_reserve_size: 0x100000, stack_commit_size: 0x1000, heap_reserve_size: 0x100000,
				heap_commit_size: 0x1000, loader_flags: 0, rva_num_and_sizes: 0x10,
				data_directory: data_dir
			}
		}
	}

	impl ImageSectionHeaderEnt {
		pub fn with_defaults(s : [u8; 8], size: u32, va: u32, raw_size: u32, ptr_raw: u32) -> ImageSectionHeaderEnt {
			ImageSectionHeaderEnt { name: s, sec_size: size, virtual_address: va, size_raw_data: raw_size, ptr_raw_data: ptr_raw,
				legacy: [0x00; 12], characteristics: 0x60000020 }
		}

		pub fn with_perms(s : [u8; 8], size: u32, va: u32, raw_size: u32, ptr_raw: u32, perms: u32) -> ImageSectionHeaderEnt {
			ImageSectionHeaderEnt { name: s, sec_size: size, virtual_address: va, size_raw_data: raw_size, ptr_raw_data: ptr_raw,
				legacy: [0x00; 12], characteristics: perms }
		}
	}

	impl ImageDataDirectory {
		#[allow(dead_code)]
		pub fn with_params(v_addr: u32, sz: u32) -> ImageDataDirectory {
			ImageDataDirectory { virtual_address: v_addr, size: sz}
		}
	}

	impl Default for ImageDataDirectory {
		fn default() -> Self {
			ImageDataDirectory { virtual_address: 0, size: 0 }
		}
	}

	impl Default for InitFileHeader {
		fn default() -> InitFileHeader {
			InitFileHeader {
				h1: [0x4D, 0x5A, 0x90, 0x00, 0x03, 0x00, 0x00, 0x00, 0x04, 0x00, 0x00, 0x00, 0xFF, 0xFF, 0x00, 0x00,
						0xB8, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x40, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
						0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
						0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0xF0, 0x00, 0x00, 0x00,
						0x0E, 0x1F, 0xBA, 0x0E, 0x00, 0xB4, 0x09, 0xCD, 0x21, 0xB8, 0x01, 0x4C, 0xCD, 0x21, 0x54, 0x68,
						0x69, 0x73, 0x20, 0x70, 0x72, 0x6F, 0x67, 0x72, 0x61, 0x6D, 0x20, 0x63, 0x61, 0x6E, 0x6E, 0x6F,
						0x74, 0x20, 0x62, 0x65, 0x20, 0x72, 0x75, 0x6E, 0x20, 0x69, 0x6E, 0x20, 0x44, 0x4F, 0x53, 0x20,
						0x6D, 0x6F, 0x64, 0x65, 0x2E, 0x0D, 0x0D, 0x0A, 0x24, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00],
				pad: [ 0x00; 0x70],
			}
		}
	}


	impl Serialize for HeaderPE {
		fn serialize(&self) -> Vec<u8> {
			let mut ret: Vec<u8> =  self.init_header.serialize();
			ret.append(&mut self.nt_header.serialize());
			for sect in (&self.section_header).iter() {
				ret.append(&mut sect.serialize());
			}
			ret
		}
	}

	impl Serialize for ImageNTHeaders32 {
		fn serialize(&self) -> Vec<u8> {
			let mut buf: Vec<u8> = self.signature.to_le_bytes().to_vec();
			buf.append(&mut self.file_header.serialize());
			buf.append(&mut self.opt_file_header.serialize());
			buf
		}
	}

	impl Serialize for ImageFileHeader {
		fn serialize(&self) -> Vec<u8> {
			let mut buf: Vec<u8> = self.machine.to_le_bytes().to_vec();
			buf.extend_from_slice(&self.num_sections.to_le_bytes());
			buf.extend_from_slice(&self.td_stamp.to_le_bytes());
			buf.extend_from_slice(&self.sym_tbl_ptr.to_le_bytes());
			buf.extend_from_slice(&self.num_symbols.to_le_bytes());
			buf.extend_from_slice(&self.opt_header_size.to_le_bytes());
			buf.extend_from_slice(&self.characteristics.to_le_bytes());
			buf
		}
	}

	impl Serialize for OptionalImageFileHeader32 {
		fn serialize(&self) -> Vec<u8> {
			let mut buf: Vec<u8> = self.magic.to_le_bytes().to_vec();
			buf.extend_from_slice(&self.linker_ver);
			buf.extend_from_slice(&self.code_size.to_le_bytes());
			buf.extend_from_slice(&self.idata_size.to_le_bytes());
			buf.extend_from_slice(&self.uninitialized_data_size.to_le_bytes());
			buf.extend_from_slice(&self.entry_point_addr.to_le_bytes());
			buf.extend_from_slice(&self.base_of_code.to_le_bytes());
			buf.extend_from_slice(&self.base_of_data.to_le_bytes());
			buf.extend_from_slice(&self.image_base.to_le_bytes());
			buf.extend_from_slice(&self.section_align.to_le_bytes());
			buf.extend_from_slice(&self.file_align.to_le_bytes());
			buf.extend_from_slice(&self.sys_ver);
			buf.extend_from_slice(&self.w32_version_value.to_le_bytes());
			buf.extend_from_slice(&self.image_size.to_le_bytes());
			buf.extend_from_slice(&self.headers_size.to_le_bytes());
			buf.extend_from_slice(&self.ck_sum.to_le_bytes());
			buf.extend_from_slice(&self.subsystem.to_le_bytes());
			buf.extend_from_slice(&self.dll_characteristics.to_le_bytes());
			buf.extend_from_slice(&self.stack_reserve_size.to_le_bytes());
			buf.extend_from_slice(&self.stack_commit_size.to_le_bytes());
			buf.extend_from_slice(&self.heap_reserve_size.to_le_bytes());
			buf.extend_from_slice(&self.heap_reserve_size.to_le_bytes());
			buf.extend_from_slice(&self.loader_flags.to_le_bytes());
			buf.extend_from_slice(&self.rva_num_and_sizes.to_le_bytes());
			for sub_struct in &self.data_directory {
				buf.append(&mut sub_struct.serialize());
			}
			buf
		}
	}

	impl Serialize for ImageSectionHeaderEnt {
		fn serialize(&self) -> Vec<u8> {
			let mut buf: Vec<u8> = self.name.to_vec();
			buf.extend_from_slice(&self.sec_size.to_le_bytes());
			buf.extend_from_slice(&self.virtual_address.to_le_bytes());
			buf.extend_from_slice(&self.size_raw_data.to_le_bytes());
			buf.extend_from_slice(&self.ptr_raw_data.to_le_bytes());
			buf.extend_from_slice(&self.legacy);
			buf.extend_from_slice(&self.characteristics.to_le_bytes());
			buf
		}
	}

	impl Serialize for InitFileHeader {
		fn serialize(&self) -> Vec<u8> {
			let mut buf: Vec<u8> = self.h1.to_vec();
			buf.extend_from_slice(&self.pad);
			buf
		}
	}

	impl Serialize for ImageDataDirectory {
		fn serialize(&self) -> Vec<u8> {
			let mut buf: Vec<u8> = self.virtual_address.to_le_bytes().to_vec();
			buf.extend_from_slice(&self.size.to_le_bytes());
			buf
		}
	}

}

#[cfg(target_os="linux")]
pub mod linux {
	use crate::executable_header::Serialize;
	use crate::file_handler::*;
	use std::io::{Read, Seek};

	const ELF32_VADDR_OFFSET: u32 = 0x08048000;
	pub fn elf_writer(params: &mut Params)  -> std::io::Result<()> {
		const ALIGN_SECTION: usize = super::SECTION_ALIGNMENT as usize;
		let mut f_sizes: Vec<u32> = Vec::new();
		let code_size: u64 = params.code_fd.metadata()?.len();
		let mut temp_reader: BufReader<&File> = BufReader::new(&params.code_fd);
		let entry_offset: u32 = if params.flag_check(BitFlags::CodeSectionPadding) { super::WORD_ALIGN_BYTES } else { 0 };
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
		let mut write_buff: [u8; ALIGN_SECTION] = [super::HEADER_PADDING_BYTE; ALIGN_SECTION];
		if buff_align < ALIGN_SECTION {
			_bytes_written += buf_writer(&mut params.writer, &write_buff, buff_align)?;
		}
		params.writer.flush()?;
		for _full_reads in 0..f_sizes[0]/super::SECTION_ALIGNMENT {
			temp_reader.read_exact(&mut write_buff).unwrap_or_else(|e| { panic!("code.bin file reader error: {}", e); });
			_bytes_written += buf_writer(&mut params.writer, &write_buff,ALIGN_SECTION)?;
		}
		let mut end_vec: Vec<u8> = Vec::new();
		let code_data_end = temp_reader.read_to_end(&mut end_vec)?;
		_bytes_written += buf_writer(&mut params.writer, end_vec.as_slice(), end_vec.len())?;
		let align = ALIGN_SECTION - (code_data_end % ALIGN_SECTION);

		if align > 0 {
			write_buff = [super::CODE_PADDING; ALIGN_SECTION];
			_bytes_written += buf_writer(&mut params.writer, &write_buff, align)?;
		}
		params.writer.flush()?;
		if f_data_valid {
			let mut data_reader: BufReader<&File> = BufReader::new(params.data_fd.as_mut().unwrap());
			data_reader.rewind()?;
			write_buff = [super::CODE_PADDING; ALIGN_SECTION];
			for _full_buf in 0..f_sizes[1]/super::SECTION_ALIGNMENT {
				data_reader.read_exact(&mut write_buff).unwrap_or_else(|e| { panic!("data.bin file reader error: {}", e); });
				_bytes_written += buf_writer(&mut params.writer, &write_buff,ALIGN_SECTION)?;
			}
			params.writer.flush()?;
			let mut end_vec: Vec<u8> = Vec::new();
			let data_end = data_reader.read_to_end(&mut end_vec)?;
			_bytes_written += buf_writer(&mut params.writer, end_vec.as_slice(), end_vec.len())?;
			let align_data = ALIGN_SECTION - (data_end % ALIGN_SECTION);
			if align_data > 0 {
				write_buff = [super::CODE_PADDING; ALIGN_SECTION];
				_bytes_written += buf_writer(&mut params.writer, &write_buff, align_data)?;
			}
			params.writer.flush()?;
		}
		Ok(())
	}
	const fn get_header_size(num_prog_headers: usize, num_sections: usize) -> usize {
		use std::mem::size_of;
		size_of::<ELF32Header>() + (size_of::<ELF32ProgHeaderEnt>()*num_prog_headers)  + (size_of::<ELF32SectHeaderEnt>() * num_sections) + 24
	}

	fn calc_alignment(len: u32, align_val: u32) -> u32 {
		(match len % align_val {
			0 => len/align_val,
			_ => (len/align_val) + 1,
		})*align_val
	}

	fn generate_elf_executable_code_only(code_len: u32, entry_byte_offset: u32) -> FileHeader {
		let phead_size: usize = get_header_size(3, 0);
		let code_align = calc_alignment(phead_size as u32, super::SECTION_ALIGNMENT);

		let p_head_list: Vec<ELF32ProgHeaderEnt> = vec![ELF32ProgHeaderEnt::new(phead_size as u32),
			ELF32ProgHeaderEnt::with_params(1, 0, code_align, super::SECTION_ALIGNMENT),
			ELF32ProgHeaderEnt::with_params_flags(1, code_align, code_len+(2*entry_byte_offset), 0x05, super::SECTION_ALIGNMENT)
		];
		//let sec_head_list: Vec<ELF32SectHeaderEnt> = vec![ELF32SectHeaderEnt::default(), ELF32SectHeaderEnt::as_code(1, 0x1000, 0, code_len), ELF32SectHeaderEnt::as_shrtrtab(2, 0, )];
		FileHeader::with_defaults(ELF32Header::new(code_align+entry_byte_offset,3), p_head_list)
		//FileHeader::with_target_arch_abi(ELF32Header::new(code_align+entry_byte_offset, 3), p_head_list, 97, 40) //ARM RISC isa
	}

	fn generate_elf_executable_header(sec_sizes: Vec<u32>, num_sections: u16, entry_byte_offset: u32) -> FileHeader {
		if num_sections == 1 {
			generate_elf_executable_code_only(sec_sizes[0], entry_byte_offset)
		}
		else {
			let phead_size: usize = get_header_size(4, 0);
			let code_len = sec_sizes[0]+ 2*entry_byte_offset;
			let code_sec_offset = calc_alignment(phead_size as u32, super::SECTION_ALIGNMENT);
			let data_sec_offset = code_sec_offset + calc_alignment(code_len, super::SECTION_ALIGNMENT);
			let p_head_list: Vec<ELF32ProgHeaderEnt> = vec![ELF32ProgHeaderEnt::new(phead_size as u32),
				ELF32ProgHeaderEnt::with_params(1, 0, code_sec_offset, super::SECTION_ALIGNMENT),
				ELF32ProgHeaderEnt::with_params_flags(1, code_sec_offset, code_len, 0x05, super::SECTION_ALIGNMENT),
				ELF32ProgHeaderEnt::with_params_flags(1, data_sec_offset, sec_sizes[1], 0x06, super::SECTION_ALIGNMENT)
			];
			//let sec_head_list: Vec<ELF32SectHeaderEnt> = vec![ELF32SectHeaderEnt::default(), ELF32SectHeaderEnt::with_params_flags(1, 1, 0x06, 0x1000, code_len, 0, 0, 16, 0)];
			FileHeader::with_defaults(ELF32Header::new(code_sec_offset+entry_byte_offset,4), p_head_list)
			//FileHeader::with_target_arch_abi(ELF32Header::new(code_sec_offset+entry_byte_offset, 4), p_head_list, 97, 40) //ARM RISC isa
		}
	}
#[derive(Clone)]
#[repr(C)]
pub struct FileHeader {
	pub magic: [u8; 4], //Always [0x7F, 0x45, 0x4C, 0x46],
	pub e_class: u8, //1 = 32-bit, 2 = 64 //EI_CLASS
	pub common: [u8; 2], // Always [0x01, 0x01], EI_DATA, EI_VERSION
	pub abi: u8, // Default to 0 (Non-system specific System-V), 97=ARM OS ABI, 64=Embedded ARM OS ABI
	pub pad: [u8; 0x08 ], // 8 byte padding [ 0x00; 0x08],
	pub ofile_type: u16, //Always 0x02 = ET_EXEC
	pub machine_isa: u16, // if 32-bit use 0x03 = x86 (Intel 80386) otherwise, 0x28 (ARM) get from system
	pub e_ver: u32, //Always 1,
	pub e_dyn_elf: ELF32Header,
	pub prog_headers: Vec<ELF32ProgHeaderEnt>,
	//pub sect_headers: Vec<ELF32SectHeaderEnt>
}
#[derive(Copy, Clone)]
#[repr(C)]
pub struct ELF32Header {
	pub e_start_addr: u32,
	pub phead_off: u32, //Always 0x34
	pub sheader_off: u32,
	pub flags: u32, // Always 0
	pub eh_size: u16, //Always 52
	pub pheadent_size: u16, // Always 32
	pub ph_num: u16,
	pub shead_size: u16, // Always 40,
	pub shead_num: u16,
	pub sh_str_ndx: u16
}

#[derive(Copy, Clone)]
#[repr(C)]
pub struct ELF32ProgHeaderEnt {
	pub pt_type: u32, // PHDR = 0x06, PT_LOAD = 1,
	pub p_off: u32, //header offset
	pub p_vaddr: u32, // Virtual start address 0x34,
	pub p_paddr: u32, // physical address
	pub p_file_size: u32,
	pub p_mem_size: u32,
	pub p_flags: u32, //default is read, 0x04 = read, 1=execute, 2= write
	pub p_align: u32 //Address alignment, typically 0x04 or 0x1000
}

#[derive(Copy, Clone)]
#[repr(C)]

pub struct ELF32SectHeaderEnt {
	pub sh_name: u32, //Section Header symbol table offset
	pub sh_type: u32, //Section Header type: 0 (SHT_NULL), 1 (SHT_PROGBITS), 3 (SHT_STRTAB)
	pub sh_flags: u32, // 0 = none, 1=write, 2 = alloc, 4 =execute
	pub sh_vaddr: u32, //va_address
	pub sh_offset: u32, //offset within section
	pub sh_size: u32, //can be 0
	pub sh_link: u32, //can be zero
	pub sh_info: u32, //Extra info about section
	pub sh_addralign: u32,
	pub sh_entsize: u32 // Can be 0
}

impl FileHeader {
	#[cfg(not(target_arch = "arm"))]
	pub fn with_defaults(e_dyn_elf: ELF32Header, prog_headers: Vec<ELF32ProgHeaderEnt>) -> Self {
		FileHeader { magic: [0x7F, 0x45, 0x4C, 0x46], e_class: 1, common: [0x01, 0x01], abi: 0, pad: [0x00; 8], ofile_type: 0x02, machine_isa: 0x03, e_ver: 1, e_dyn_elf, prog_headers }
	}

	#[cfg(target_arch = "arm")]
	pub fn with_defaults(e_dyn_elf: ELF32Header, prog_headers: Vec<ELF32ProgHeaderEnt>) -> Self {
		FileHeader { magic: [0x7F, 0x45, 0x4C, 0x46], e_class: 1, common: [0x01, 0x01], abi: 0x61, pad: [0x00; 8], ofile_type: 0x02, machine_isa: 0x28, e_ver: 1, e_dyn_elf, prog_headers }
	}

	#[allow(dead_code)]
	pub fn with_target_arch_abi(e_dyn_elf: ELF32Header, prog_headers: Vec<ELF32ProgHeaderEnt>, abi: u8, isa: u16) -> Self {
		FileHeader { magic: [0x7F, 0x45, 0x4C, 0x46], e_class: 1, common: [0x01, 0x01], abi, pad: [0x00; 8], ofile_type: 0x02, machine_isa: isa, e_ver: 1, e_dyn_elf, prog_headers }
	}

}

impl ELF32Header {
	pub fn new(start_off: u32, num_ph: u16) -> Self {
		ELF32Header { e_start_addr: start_off+ELF32_VADDR_OFFSET, phead_off: 0x34, sheader_off: 0, flags: 0, eh_size: 52,
			pheadent_size: 32, ph_num: num_ph, shead_size: 0, shead_num: 0, sh_str_ndx: 0
		}
	}

	#[allow(dead_code)]
	pub fn with_sheader(start_off: u32, sh_off: u32, num_ph: u16, num_sh: u16) -> Self {
		ELF32Header { e_start_addr: start_off+ELF32_VADDR_OFFSET, phead_off: 0x34, sheader_off: sh_off, flags: 0, eh_size: 52,
			pheadent_size: 32, ph_num: num_ph, shead_size: 40, shead_num: num_sh, sh_str_ndx: num_sh-1
		}
	}

	//pub fn new(start_addr: u32, shead_off: u32, num_ph: u16, num_sh: u16, sh_str_idx: u16) -> Self {
	//	ELF32Header { e_start_addr: start_addr, phead_off: 0x34, sheader_off: shead_off, flags: 0, eh_size: 52,
	//		pheadent_size: 32, ph_num: num_ph, shead_size: 40, shead_num: num_sh, sh_str_ndx: sh_str_idx
	//	}
	//}
}

impl ELF32ProgHeaderEnt {
	pub fn new(size: u32) -> Self { //PHDR ONLY
		ELF32ProgHeaderEnt { pt_type: 0x06, p_off: 0x34,
			p_vaddr: ELF32_VADDR_OFFSET+0x34, p_paddr: ELF32_VADDR_OFFSET+0x34,
			p_file_size: size, p_mem_size: size, p_flags: 0x04, p_align: 0x4
		}
	}

	pub fn with_params(kind: u32, off: u32, size: u32, align: u32) -> Self {
		ELF32ProgHeaderEnt { pt_type: kind, p_off: off, p_vaddr:ELF32_VADDR_OFFSET+off, p_paddr: ELF32_VADDR_OFFSET+off, p_file_size: size, p_mem_size: size, p_flags: 0x4, p_align: align }
	}

	pub fn with_params_flags(kind: u32, off: u32, size: u32, flags: u32, align: u32) -> Self {
		ELF32ProgHeaderEnt { pt_type: kind, p_off: off, p_vaddr: ELF32_VADDR_OFFSET+off, p_paddr: ELF32_VADDR_OFFSET+off, p_file_size: size, p_mem_size: size, p_flags: flags, p_align: align }
	}
}

impl ELF32SectHeaderEnt {
	pub fn as_shrtrtab(sh_name: u32, sh_addr_base: u32, sh_offset: u32, sh_size: u32) -> Self {
		ELF32SectHeaderEnt { sh_name, sh_type: 3, sh_flags: 0, sh_vaddr: sh_addr_base, sh_offset, sh_size, sh_link: 0, sh_info: 0, sh_addralign: 1, sh_entsize: 0}
	}

	pub fn as_code(sh_name: u32, sh_addr_base: u32, sh_offset: u32, sh_size: u32) -> Self {
		ELF32SectHeaderEnt { sh_name, sh_type: 1, sh_flags: 6, sh_vaddr: sh_addr_base+ELF32_VADDR_OFFSET, sh_offset, sh_size, sh_link: 0, sh_info: 0, sh_addralign: 16, sh_entsize: 0}
	}

	pub fn as_data(sh_name: u32, sh_addr_base: u32, sh_offset: u32, sh_size: u32) -> Self {
		ELF32SectHeaderEnt { sh_name, sh_type: 1, sh_flags: 5, sh_vaddr: sh_addr_base+ELF32_VADDR_OFFSET, sh_offset, sh_size, sh_link: 0, sh_info: 0, sh_addralign: 4, sh_entsize: 0}
	}

}

impl Default for ELF32SectHeaderEnt {
	fn default() -> Self { //For PHDR ONLY
		ELF32SectHeaderEnt { sh_name: 0, sh_type: 0, sh_flags: 0, sh_vaddr: 0, sh_offset: 0, sh_size: 0,
			sh_link: 0, sh_info: 0, sh_addralign: 0, sh_entsize: 0
		}
	}
}

impl Serialize for FileHeader {
	fn serialize(&self) -> Vec<u8> {
		let mut ret: Vec<u8> = self.magic.to_vec();
		ret.push(self.e_class);
		ret.extend_from_slice(&self.common);
		ret.push(self.abi);
		ret.extend_from_slice(&self.pad);
		ret.extend_from_slice(&self.ofile_type.to_le_bytes());
		ret.extend_from_slice(&self.machine_isa.to_le_bytes());
		ret.extend_from_slice(&self.e_ver.to_le_bytes());
		ret.append(&mut self.e_dyn_elf.serialize());
		for p_entry in (&self.prog_headers).iter() {
			ret.append(&mut p_entry.serialize());
		}
		//for s_entry in (&self.sect_headers).iter() {
		//	ret.append(&mut s_entry.serialize());
		//}
		ret
	}
}

impl Serialize for ELF32SectHeaderEnt {
	fn serialize(&self) -> Vec<u8> {
		let mut ret: Vec<u8> = self.sh_name.to_le_bytes().to_vec();
		ret.extend_from_slice(&self.sh_type.to_le_bytes());
		ret.extend_from_slice(&self.sh_flags.to_le_bytes());
		ret.extend_from_slice(&self.sh_vaddr.to_le_bytes());
		ret.extend_from_slice(&self.sh_offset.to_le_bytes());
		ret.extend_from_slice(&self.sh_size.to_le_bytes());
		ret.extend_from_slice(&self.sh_link.to_le_bytes());
		ret.extend_from_slice(&self.sh_info.to_le_bytes());
		ret.extend_from_slice(&self.sh_addralign.to_le_bytes());
		ret.extend_from_slice(&self.sh_entsize.to_le_bytes());
		ret
	}
}

impl Serialize for ELF32ProgHeaderEnt {
	fn serialize(&self) -> Vec<u8> {
		let mut ret: Vec<u8> = self.pt_type.to_le_bytes().to_vec();
		ret.extend_from_slice(&self.p_off.to_le_bytes());
		ret.extend_from_slice(&self.p_vaddr.to_le_bytes());
		ret.extend_from_slice(&self.p_paddr.to_le_bytes());
		ret.extend_from_slice(&self.p_file_size.to_le_bytes());
		ret.extend_from_slice(&self.p_mem_size.to_le_bytes());
		ret.extend_from_slice(&self.p_flags.to_le_bytes());
		ret.extend_from_slice(&self.p_align.to_le_bytes());
		ret
	}
}

impl Serialize for ELF32Header {
	fn serialize(&self) -> Vec<u8> {
		let mut ret: Vec<u8> = self.e_start_addr.to_le_bytes().to_vec();
		ret.extend_from_slice(&self.phead_off.to_le_bytes());
		ret.extend_from_slice(&self.sheader_off.to_le_bytes());
		ret.extend_from_slice(&self.flags.to_le_bytes());
		ret.extend_from_slice(&self.eh_size.to_le_bytes());
		ret.extend_from_slice(&self.pheadent_size.to_le_bytes());
		ret.extend_from_slice(&self.ph_num.to_le_bytes());
		ret.extend_from_slice(&self.shead_size.to_le_bytes());
		ret.extend_from_slice(&self.shead_num.to_le_bytes());
		ret.extend_from_slice(&self.sh_str_ndx.to_le_bytes());
		ret
	}
}

}


