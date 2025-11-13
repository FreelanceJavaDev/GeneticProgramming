# GeneticProgramming
 Experiment transforming genomic sequences into executables.

## The concept
This project began as several years ago as a joke about DNA being "God's change log". But the more I thought about it the more DNA felt like a program to me. I soon forgot about it other than as a science fiction plot device. I decided to give it a chance as a way to learn Rust that kept my interest and leverage what I learned at a genetic sequencing company in an entirely different concept.

When I started this project I fully expected that I would get the answer "No DNA doesn't contain computer code" I expected random data. Once I saw actual regions of assembly code that could be recognized by a disassembler coming out as a result I became intrigued to see if different methods of parsing and producing an executable that _could_ theoretically run on a system. The one thing I don't intend to do is inject any code that would alter the resulting executable's behavior. I do use single byte opcodes for section and file alignments.

## Interesting coincidences and DNA data encoding of digital data
There are some interesting correlations between DNA and binary encoding systems.
While binary is a base 2 system consisting of 0 and 1. DNA can be considered a base 4 system consisting of A, T, C and G. It only takes 2 bits to represent each nucleotide uniquely as 00 (0), 01 (1), 10(2) and 11 (3).
The remaining correlations requires some specific background knowledge about DNA sequencing, analysis and processing.

There are several categories for DNA. Broadly speaking there is encoding DNA and non-encoding DNA. Encoding DNA is DNA that is used by RNA transcribers to generate amino-acid protein chains. The way a specific amino acid is chosen is directly related to a DNA base pair triplet referred to as a codon. If individual base pairs are a letter, a codon would be analogous to a word. There are a total of 4^3^ unique codons. Each encoding DNA sequence is contained between a start codon and a stop codon. There are three stop codons: TAA TGA and TAG. The start codons are a bit more complex and vary between organisms, but the one that is universally agreed to is a single codon ATG. The primary alternates are typically used for single organisms are TTG, GTG and CTG. All codons except for the stop codons map to one of twenty amino-acids.
Non-encoding DNA is the DNA before a start codon and after a stop codon.

Now for the correlations to computers. Using the 2-bit nucleotide mapping a codon would occupy 6 bits, meaning a codon can reside within a byte with two bits (one nucleotide to spare). The number of codons is  is even nicer 4^3^ comes out to 64.

It doesn't end there, just as there is overlap in dna encoding of amino acids, so too is there overlap in computer machine op-codes. On x86 0x60-0x6F are copies of the conditional jumps typically 0x70-0x7F, 0xF0 and 0xF1 both act as a lock prefix.

## Compilation methods
All the methods center around converting a base-4 code (A, T, C & G) into a base 2 code (0 & 1). They all involve using two bits for each base.
Method 1. The ASCII masking & bit shift trick
	This is the most obvious one to software engineers.
	From the ASCII table A, C, T and G can be uniquely mapped to two bits by the bitwise expression `(Q & 0x06) >> 1`.

	This maps A => 00, C => 01, T => 10 and G => 11.
	However this only applies to capital letters.

Method 2. The Biological base pair complement map
	This method is one biologists would prefer. It also makes the most sense for a program encoded in DNA. In DNA most of the time A binds with T and C binds with G. These pairings can be seen as complemental pairings.
	This does require an ASCII table to map the sequences into the proper mapping. Though this mapping can be modified easily. The first will seem a bit like the first method.
	The 2-bit mappings that are currently implemented are:
        M2A0C1: A => 00, C => 01, G => 10 and T => 11
        M2A1C0: A => 01, C => 00, G => 11 and T => 10
        M2A2C3: A => 10, C => 11, G => 00 and T => 01
        M2A3C2: A => 11, C => 10, G => 01 and T => 00

Method 3. Hydrogen Bond Number (Partially implemented)
    This method uses hydrogen bonds between the base pairs as the differentiating factor. A and T have 2 hydrogen bonds whereas C and G have three. This is why analysis of DNA frequently refers to the CG content percentage as a DNA stability metric. It also can be used to map DNA base pairs to a single bit.

Method 4. The ASCII Table dynamic map. (unimplemented)
	This mapping is a variation on the biological base pair compliment map and the remaining mappings that aren't covered under the previous methods.


## Testing and Examining results
It is expected nobody should run a program of unknown origin or purpose their computer without examining it first.
The conventional method is via a disassembler. This is used to examine computer virus code and determine what they do. Cutter was used for examining the resulting executables. The first organism tested is in the family with the oldest form of life Cyanobacteria. Specifically a _Cylindrospermopsis raciborskii_ strain. It also has everything in a single chromosome.

Initially the sequence was turned directly into an executable and treated as code. This defaulted to 64-bit Windows assembly on AMD on the test computer. This produced an executable file, but it would never run.  This is due to the executable file header being absent.  In 64-bit interpretations there was a lot of assembly code, but there were some invalid instructions. This was due to bad prefix codes that are specific for 64-bit assembly, these went away when I specified 32-bit assembly.

An initial attempt using codons didn't lead to viable code. A second attempt will be made after the the PE and ELF header formats are implemented so data and code segments can be handled separately then merged.

Currently only the windows PE header file is fully implemented and tested for 32-bit.


## Future Improvements:
1. Linux ELF file formats for 32 bit x86 ISA
2. Linux ELF file formats for ARM 32-bit ISA (for RISC)
3. After initial raw executable generated validate the executable for problematic/invalid op codes (0x0F is a dangerous single byte opcode on x86, but 0xF0 is a lock prefix)
4. Parsing codons while handling problematic/invalid opcodes that result from the parsing.