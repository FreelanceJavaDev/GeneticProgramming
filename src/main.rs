use std::env;
use std::fs::File;
use std::fs::OpenOptions;
use std::io::BufRead;
use std::io::BufReader;
use std::io::BufWriter;
use std::io::Write;
use std::path::PathBuf;

enum OutputFormat {
	Win,
	Linux,
	Unsupported,
}

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

fn main() -> std::io::Result<()> {
	let mut table: [u8; 256] = [255; 256];
	table['A' as usize] = 0;
	table['a' as usize]=0;
	table['C' as usize] = 1; 
	table['c' as usize]=1;
	table['G' as usize] = 2;
	table['g' as usize]=2;
	table['T' as usize] = 3;
	table['t' as usize]=3;

	let main_args: Vec<String> = env::args().collect();
	if main_args.len() < 2 {
		panic!("Error: Minimum args is 2. Use --help.");
	}
	else if main_args.len() == 2 {
		if main_args[1] != "--help" { panic!("Invalid argument."); }
		else { println!("Usage: <prog_path> --input=<path to fastfa file> --outputFmt=[win|linux]"); }
	}
	else if main_args.len() >= 3 { //.\test_data\GCA_003367075_2_ASM336707v2_genomic.fna
		let arg_parse: Vec<&str> = main_args[1].split('=').collect();
		let in_file = File::open(arg_parse[1]).unwrap_or_else(|e|{panic!("Failed to open input file: {}", e);});
		let mut reader = BufReader::new(in_file);
		let mut buffer= String::new();
		let curr_dir: PathBuf = env::current_dir()?;
		println!("Current directory: {}", curr_dir.display());
		let arg_parse2: Vec<&str> = main_args[2].split('=').collect();
		let file_ext: OutputFormat= match arg_parse2[1] {
			"win" => OutputFormat::Win,
			"linux" => OutputFormat::Linux,
			_ => OutputFormat::Unsupported,
		};
	
		reader.read_line(&mut buffer)?;
		println!("First Line from input: {}", buffer);
		if buffer.chars().nth(0).unwrap_or_default() == '>' {
			buffer.remove(0);
		}
		let mut end = buffer.trim_ascii_end().len();
		let end_early: Option<usize> = buffer.find(|c: char| (c == '.') || (c.is_whitespace()));
		if end_early.is_some() { end = end_early.unwrap(); }
		buffer.truncate(end);
		println!("Output filename: {}", buffer);

		let mut outfile = PathBuf::from(".");
		outfile.push("test_data");
		outfile.push(&buffer);
		match file_ext {
			OutputFormat::Win => { outfile.set_extension("exe") },
			OutputFormat::Linux => { outfile.set_extension("bin") },
			_ => { outfile.set_extension("");
			panic!("Only windows and linux executable file formats are supported for output") 
			},
		};

		println!("Output file is: {}", outfile.as_path().display());
		buffer.clear();
		let file_out: File = (match outfile.exists() {
			true => OpenOptions::new().write(true).truncate(true).read(true).open(outfile.as_path()),
			false => OpenOptions::new().write(true).create(true).read(true).open(outfile.as_path()),
		}).unwrap_or_else(|e|{ panic!("Failed to create output file: {}", e);});
		let mut writer: BufWriter<File> = BufWriter::new(file_out);
		buf_line_reader(&mut reader, &mut buffer).unwrap_or_else(|e|{panic!("Failed to read line in open input file: {}", e);});
		println!("Second line is: {}\n", buffer);
		let mut buff_2bit: [u8; 1000] = [0; 1000];
		let mut buff_2bit_temp: u8 = 0;
		let mut buff_2bit_size: usize = 0;
		let mut temp_bits_used: u8 = 0;
		let mut written_bytes: usize = 0;
		// println!("{} ->", buffer);
		while !buffer.is_empty() {
			for c in buffer.bytes() {
			// print!("{:0>2b}", table[c as usize]);
			buff_2bit_temp <<= 2;
			// Method 1: buff_2bit_temp += (c & 0x06) >> 1; //table[c as usize]; //<< temp_bits_used;
			buff_2bit_temp += table[c as usize]; // Method 2
			temp_bits_used += 2;
			if temp_bits_used == 8 {
				// print!(" ");
				// print!("{:0>8b} ", buff_2bit_temp);
				buff_2bit[buff_2bit_size] = buff_2bit_temp;
				temp_bits_used = 0;
				buff_2bit_temp = 0;
				buff_2bit_size += 1;
				if buff_2bit_size == 1000 {
					while written_bytes < buff_2bit_size-1 {
					let n = writer.write(&buff_2bit[written_bytes..999])?;
					written_bytes += n;
					}
					buff_2bit_size = 0;
					written_bytes = 0;
					
					//TODO: write the buffer to the file until empty.
				}
			}
			}
			buffer.clear();
			buf_line_reader(&mut reader, &mut buffer).unwrap_or_else(|e|{panic!("Failed to read line in open input file: {}", e);});
		}
		if temp_bits_used < 8 {
			// print!("{:b}", buff_2bit_temp);
			buff_2bit[buff_2bit_size] = buff_2bit_temp;
			temp_bits_used = 0;
			buff_2bit_temp = 0;
			buff_2bit_size += 1;
			while written_bytes < buff_2bit_size-1 {
			let n = writer.write(&buff_2bit[written_bytes..buff_2bit_size-1])?;
			written_bytes += n;
			}
			buff_2bit_size = 0;
			written_bytes = 0;
			
			//TODO: write the buffer to the file until empty.
		}
		// println!("");
		writer.flush()?;
	}
	Ok(())
}

fn buf_line_reader<R: BufRead>(reader: &mut R, buf: &mut String) -> std::io::Result<usize> {
	let mut len: usize = reader.read_line(buf)?;
	if !buf.is_empty() { 
		buf.pop();
		len -= 1;
	}
	
	if len > 0 && buf.chars().last().unwrap() == '\r' { 
		buf.pop();
		len -= 1;
	}
	Ok(len)
}
