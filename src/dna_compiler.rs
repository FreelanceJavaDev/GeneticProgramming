pub fn gen_binary_table(bp_map_val: [u8; 4]) -> [u8; 256] {
	 const BP_MAP: [usize; 8] = ['A' as usize, 'a' as usize, 'C' as usize, 'c' as usize, 'G' as usize, 'g' as usize, 'T' as usize, 't' as usize,];
	let mut ret: [u8; 256] = [255; 256];
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

