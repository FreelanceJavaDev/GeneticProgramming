use std::env;
//TODO use rust cfg target_arch and target_os
//TODO: retrieve processor architecture
// Windows env string: PROCESSOR_ARCHITECTURE
// Possible outputs: AMD64, IA64, ARM64 or x86 - use x86 for AMD64, IA64, x86-64 or x86, ARM used for ARM64 (TODO)
//TODO: retrieve processor architecture on linux (env var and possible outputs)
/**
 * X86 opcode format: Maximum of 15 byte instructions allowed
 * Optional prefixes:
 *	Instruction Prefix (0 or 1 byte)
 *  Address-size Prefix (0 or 1 byte)
 *  Operand-size Prefix (0 or 1 byte)
 *  Segment Override (0 or 1 byte)
 *
 * General instruction format:
 * 	Opcode (1 or 2 bytes)
 *  Mode/reg/mem (0 or 1 byte)
 * 			bit: |0 1 2|3 4 5| 6 7|
 * 		Meaning: | R/M |R/Opc|Mode|
 *	Scaled Index Byte (0 or 1 byte)
 *			bit: |0 1 2|3 4 5| 6 7 |
 * 		Meaning: |Base |Index|Scale|
 *	Displacement (0, 1, 2 or 4 bytes)
 *	Immediate (0, 1, 2 or 4 bytes)
 */

#[cfg(not(target_arch="arm"))]
#[repr(C)]
pub struct OpCodeValidator {

}






//ARM 32-bit all instructions are 32 bits long, but field sizes can vary
#[cfg(target_arch="arm")]
#[repr(u8)]
pub enum CondBits {
	EQ, // ==
	NE, // !=
	CS, // If carry flag set.
    CC, // If Carry flag cleared.
    MI, // Minus, negative.
    PL, // Plus, positive.
    VS, // Overflow flag set.
    VC, // Overflow Flag cleared.
    HI, // Unsigned higher.
    LS, // Lower or Same, unsigned.
    GE, // Greater than or equal, signed.
    LT, // Less than, signed.
    GT, // Greater or equal, signed.
    LE, // Less than or equal, signed.
    AL, // Always do (no condition).
    NV  // RESERVED, used to extend opcodes in later ARMs.
}