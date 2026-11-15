//! SHA-256 helpers matching ``samio._lib_id`` / ``_qhash``.

use sha2::{Digest, Sha256};

/// First 8 hex digits of sha256(qname) as integer (Python int(..., 16)).
pub fn qhash(qname: &str) -> i64 {
    let digest = Sha256::digest(qname.as_bytes());
    let hex = hex_encode(&digest[..4]);
    i64::from_str_radix(&hex, 16).unwrap_or(0)
}

/// ``int(sha256(rg).hexdigest()[:4], 16) % 1000``.
pub fn lib_id(rg_or_lb: &str) -> u16 {
    let digest = Sha256::digest(rg_or_lb.as_bytes());
    let hex = hex_encode(&digest[..2]);
    (u32::from_str_radix(&hex, 16).unwrap_or(0) % 1000) as u16
}

fn hex_encode(bytes: &[u8]) -> String {
    const HEX: &[u8; 16] = b"0123456789abcdef";
    let mut s = String::with_capacity(bytes.len() * 2);
    for &b in bytes {
        s.push(HEX[(b >> 4) as usize] as char);
        s.push(HEX[(b & 0xf) as usize] as char);
    }
    s
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn qhash_matches_python_shape() {
        // Smoke: stable and in u32 range when cast from 8 hex digits.
        let h = qhash("READ1");
        assert!(h >= 0);
        assert!(h <= 0xffff_ffff);
    }
}
