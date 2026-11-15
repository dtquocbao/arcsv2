//! Block names. orig:helper.py::block_idx_to_name

const BLOCK_CHARS: usize = 52;
const LOWER_OFFSET: u8 = 6;

pub fn block_idx_to_name(block_idx: usize) -> String {
    let mult = block_idx / BLOCK_CHARS;
    let mut mod_ = (block_idx % BLOCK_CHARS) as u8;
    if mod_ > 25 {
        mod_ += LOWER_OFFSET;
    }
    let letter = (b'A' + mod_) as char;
    if mult == 0 {
        letter.to_string()
    } else {
        format!("{letter}{mult}")
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn names_wrap_the_alphabet() {
        assert_eq!(block_idx_to_name(0), "A");
        assert_eq!(block_idx_to_name(25), "Z");
        assert_eq!(block_idx_to_name(26), "a");
        assert_eq!(block_idx_to_name(51), "z");
        assert_eq!(block_idx_to_name(52), "A1");
    }
}
