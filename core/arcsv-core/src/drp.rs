//! Discordant-pair compatibility and largest cliques (§4.3c).
//! orig:pecluster.py::cluster_handle_component

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Drp {
    pub a: i64,
    pub b: i64,
    pub index: usize,
}

pub fn compatible(p: Drp, q: Drp, max_distance: i64) -> bool {
    (p.a - q.a).abs() < max_distance
        && (p.b - q.b).abs() < max_distance
        && p.a.max(q.a) < p.b.min(q.b)
}

/// Largest cliques. Ties break toward the lexicographically smallest index set.
/// Returns whether the node cap was hit.
pub fn largest_cliques(pairs: &[Drp], max_distance: i64, node_cap: usize) -> (Vec<Vec<usize>>, bool) {
    let n = pairs.len();
    if n > node_cap {
        return (Vec::new(), true);
    }
    let mut adj = vec![vec![false; n]; n];
    for i in 0..n {
        for j in (i + 1)..n {
            if compatible(pairs[i], pairs[j], max_distance) {
                adj[i][j] = true;
                adj[j][i] = true;
            }
        }
    }
    let mut best: Vec<Vec<usize>> = Vec::new();
    let mut best_size = 0usize;
    fn rec(
        adj: &[Vec<bool>],
        r: &mut Vec<usize>,
        p: &mut Vec<usize>,
        x: &mut Vec<usize>,
        best: &mut Vec<Vec<usize>>,
        best_size: &mut usize,
    ) {
        if p.is_empty() && x.is_empty() {
            if r.len() > *best_size {
                *best_size = r.len();
                best.clear();
                let mut clique = r.clone();
                clique.sort();
                best.push(clique);
            } else if r.len() == *best_size && r.len() >= 2 {
                let mut clique = r.clone();
                clique.sort();
                best.push(clique);
            }
            return;
        }
        let snapshot: Vec<usize> = p.clone();
        for v in snapshot {
            if !p.contains(&v) {
                continue;
            }
            r.push(v);
            let p2: Vec<usize> = p.iter().copied().filter(|u| adj[v][*u]).collect();
            let x2: Vec<usize> = x.iter().copied().filter(|u| adj[v][*u]).collect();
            let mut p2 = p2;
            let mut x2 = x2;
            rec(adj, r, &mut p2, &mut x2, best, best_size);
            r.pop();
            p.retain(|u| *u != v);
            x.push(v);
        }
    }
    let mut r = Vec::new();
    let mut p: Vec<usize> = (0..n).collect();
    let mut x = Vec::new();
    rec(&adj, &mut r, &mut p, &mut x, &mut best, &mut best_size);
    best.retain(|c| c.len() >= 2);
    best.sort();
    (best, false)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn triangle_is_one_clique() {
        let pairs = [
            Drp { a: 100, b: 400, index: 2 },
            Drp { a: 110, b: 410, index: 0 },
            Drp { a: 90, b: 390, index: 1 },
        ];
        let (cliques, hit) = largest_cliques(&pairs, 50, 40);
        assert!(!hit);
        assert_eq!(cliques, vec![vec![0, 1, 2]]);
    }

    #[test]
    fn order_of_input_does_not_change_indexes() {
        let a = [
            Drp { a: 0, b: 100, index: 0 },
            Drp { a: 10, b: 110, index: 1 },
        ];
        let b = [a[1], a[0]];
        let (c1, _) = largest_cliques(&a, 50, 40);
        let (c2, _) = largest_cliques(&b, 50, 40);
        assert_eq!(c1, c2);
    }
}
