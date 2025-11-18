use std::env;
pub enum ProcessorEnv{
	Unsupported =-1,
	WIA64,
	WAMD64,
	WARM64,
	Wx86,
}
//TODO use rust cfg target_arch and target_os
//TODO: retrieve processor architecture
// Windows env string: PROCESSOR_ARCHITECTURE
// Possible outputs: AMD64, IA64, ARM64 or x86
//TODO: retrieve processor architecture on linux (env var and possible outputs)

//x86 has two important bytes: 1st byte = is opcode or a opcode group specifier,
//2nd byte is Mod R/M byte, if group specifier then the reg field holds the opcode specifier
//



#[cfg(target_arch="arm")]
#[repr(u8)]
pub enum CondBits {
	EQ, // ==
	NE, // !=
	CS, // If carry flag set.
    CC, // If Carry flag cleared.
    MI, // Minus, negitive.
    PL, // Plus, positive.
    VS, // Overflow flag set.
    VC, // Overflow Flag cleared.
    HI, // Unsigned higher.
    LS, // Lower or Same, unsigned.
    GE, // Greater than or equal, signed.
    LT, // Less than, signed.
    GT, // Greater or equal, signed.
    LE, // Less than or equal, signed.
    AL, // Always do.
    NV  // RESERVED, used to extend opcodes in later ARMs.
}