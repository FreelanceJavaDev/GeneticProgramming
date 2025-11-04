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
