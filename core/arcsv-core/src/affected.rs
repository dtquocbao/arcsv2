//! Affected length. Needleman–Wunsch, match 1000, mismatch -1, gap -1.
//! Every co-optimal move is kept. The alignment is repeated with the strings
//! swapped, and a block is affected if any alignment leaves it unmatched.

#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum Token {
    Block(u32, i8),
}

pub fn align(s1: &[Token], s2: &[Token]) -> (Vec<usize>, Vec<usize>) {
    let n = s1.len();
    let m = s2.len();
    let mut score = vec![vec![0i64; m + 1]; n + 1];
    let mut tb: Vec<Vec<Vec<u8>>> = vec![vec![Vec::new(); m + 1]; n + 1];
    for i in 1..=n {
        score[i][0] = -(i as i64);
        tb[i][0] = vec![1];
    }
    for j in 1..=m {
        score[0][j] = -(j as i64);
        tb[0][j] = vec![2];
    }
    for i in 1..=n {
        for j in 1..=m {
            let diag_score = if s1[i - 1] == s2[j - 1] { 1000 } else { -1 };
            let diag = score[i - 1][j - 1] + diag_score;
            let up = score[i - 1][j] - 1;
            let left = score[i][j - 1] - 1;
            let best = diag.max(up).max(left);
            let mut hows = Vec::new();
            if diag == best {
                hows.push(0);
            }
            if up == best {
                hows.push(1);
            }
            if left == best {
                hows.push(2);
            }
            score[i][j] = best;
            tb[i][j] = hows;
        }
    }
    let mut aff1 = vec![false; n];
    let mut aff2 = vec![false; m];
    fn walk(
        i: usize,
        j: usize,
        s1: &[Token],
        s2: &[Token],
        tb: &[Vec<Vec<u8>>],
        aff1: &mut [bool],
        aff2: &mut [bool],
    ) {
        if i == 0 && j == 0 {
            return;
        }
        for &how in &tb[i][j] {
            if how == 0 {
                if s1[i - 1] != s2[j - 1] {
                    aff1[i - 1] = true;
                    aff2[j - 1] = true;
                }
                walk(i - 1, j - 1, s1, s2, tb, aff1, aff2);
            } else if how == 1 {
                aff1[i - 1] = true;
                walk(i - 1, j, s1, s2, tb, aff1, aff2);
            } else {
                aff2[j - 1] = true;
                walk(i, j - 1, s1, s2, tb, aff1, aff2);
            }
        }
    }
    walk(n, m, s1, s2, &tb, &mut aff1, &mut aff2);
    let left: Vec<usize> = aff1
        .iter()
        .enumerate()
        .filter(|(_, on)| **on)
        .map(|(i, _)| i)
        .collect();
    let right: Vec<usize> = aff2
        .iter()
        .enumerate()
        .filter(|(_, on)| **on)
        .map(|(i, _)| i)
        .collect();
    (left, right)
}

pub fn affected_length(reference: &[Token], alt: &[Token], lengths: &[i64]) -> i64 {
    let (r1, a1) = align(reference, alt);
    let (a2, r2) = align(alt, reference);
    let mut ids = vec![false; lengths.len()];
    for i in r1.into_iter().chain(r2) {
        let Token::Block(id, _) = reference[i];
        ids[id as usize] = true;
    }
    for j in a1.into_iter().chain(a2) {
        let Token::Block(id, _) = alt[j];
        ids[id as usize] = true;
    }
    ids.iter()
        .enumerate()
        .filter(|(_, on)| **on)
        .map(|(i, _)| lengths[i])
        .sum()
}

#[cfg(test)]
mod tests {
    use super::*;

    fn block(id: u32, strand: i8) -> Token {
        Token::Block(id, strand)
    }

    #[test]
    fn inversion_and_deletion() {
        let lengths = [100, 50, 150, 25, 100];
        let reference: Vec<_> = (0..5).map(|i| block(i, 1)).collect();
        assert_eq!(affected_length(&reference, &reference, &lengths), 0);
        let inverted = [block(0, 1), block(1, 1), block(2, 1), block(3, -1), block(4, 1)];
        assert_eq!(affected_length(&reference, &inverted, &lengths), 25);
        let deletion = [block(0, 1), block(4, 1)];
        assert_eq!(affected_length(&reference, &deletion, &lengths), 225);
    }
}
