//! CIGAR helpers matching ``arcsv2.engine.evidence``.

/// Parse CIGAR into (length, op) pairs. Unknown ops are kept as-is.
pub fn cigar_ops(cigar: &str) -> Vec<(u32, char)> {
    let mut out = Vec::new();
    if cigar.is_empty() || cigar == "*" {
        return out;
    }
    let mut num: u32 = 0;
    let mut seen_digit = false;
    for ch in cigar.chars() {
        if ch.is_ascii_digit() {
            num = num
                .saturating_mul(10)
                .saturating_add(ch.to_digit(10).unwrap_or(0));
            seen_digit = true;
        } else {
            if seen_digit {
                out.push((num, ch));
            }
            num = 0;
            seen_digit = false;
        }
    }
    out
}

pub fn aligned_length(cigar: &str) -> i64 {
    cigar_ops(cigar)
        .into_iter()
        .filter(|(_, op)| matches!(op, 'M' | 'D' | 'N' | '=' | 'X'))
        .map(|(n, _)| n as i64)
        .sum()
}

pub fn query_len_from_cigar(cigar: &str) -> i64 {
    cigar_ops(cigar)
        .into_iter()
        .filter(|(_, op)| matches!(op, 'M' | 'I' | 'S' | '=' | 'X'))
        .map(|(n, _)| n as i64)
        .sum()
}

/// Leading and trailing soft-clip lengths. Hard clips are skipped.
pub fn softclip_counts(cigar: &str) -> (i64, i64) {
    let ops = cigar_ops(cigar);
    if ops.is_empty() {
        return (0, 0);
    }
    let mut i = 0usize;
    while i < ops.len() && ops[i].1 == 'H' {
        i += 1;
    }
    let left = if i < ops.len() && ops[i].1 == 'S' {
        ops[i].0 as i64
    } else {
        0
    };
    let mut j = ops.len() as isize - 1;
    while j >= 0 && ops[j as usize].1 == 'H' {
        j -= 1;
    }
    let right = if j >= 0 && ops[j as usize].1 == 'S' && (j as usize) != i {
        ops[j as usize].0 as i64
    } else {
        0
    };
    (left, right)
}

pub fn cigar_has_softclip(cigar: &str) -> bool {
    !cigar.is_empty() && cigar.contains('S')
}
