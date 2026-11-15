//! Iterative path enumeration. orig:sv_inference.py::get_paths_iterative
//!
//! Neighbor ids are sorted so the set of paths does not depend on insertion order.
//! A locus with more than `max_paths` traversals is skipped.

#[derive(Clone, Debug)]
pub struct Graph {
    pub n_blocks: usize,
    /// Undirected adjacency edges, not the implicit block edges.
    pub neighbors: Vec<Vec<usize>>,
}

impl Graph {
    pub fn new(n_blocks: usize) -> Self {
        Self {
            n_blocks,
            neighbors: vec![Vec::new(); n_blocks * 2],
        }
    }

    pub fn add_edge(&mut self, u: usize, v: usize) {
        if !self.neighbors[u].contains(&v) {
            self.neighbors[u].push(v);
            self.neighbors[u].sort();
        }
        if !self.neighbors[v].contains(&u) {
            self.neighbors[v].push(u);
            self.neighbors[v].sort();
        }
    }
}

fn flip(node: usize) -> usize {
    if (node & 1) == 0 {
        node + 1
    } else {
        node - 1
    }
}

pub fn enumerate_paths(graph: &Graph, max_paths: usize) -> Result<Vec<Vec<usize>>, &'static str> {
    let end = graph.n_blocks * 2 - 1;
    let mut paths = Vec::new();
    let mut path = Vec::new();
    let mut cur = 0usize;
    let mut stacks: Vec<Vec<usize>> = Vec::new();
    let mut guard = 0usize;
    let guard_max = 20 * max_paths.max(1) + 100;
    loop {
        guard += 1;
        if guard > guard_max {
            return Err("skipped_backtrack");
        }
        let nxt = flip(cur);
        path.push(cur);
        path.push(nxt);
        if nxt == end {
            paths.push(path.clone());
            if paths.len() > max_paths {
                return Err("skipped_npaths");
            }
        }
        let mut neigh = graph.neighbors.get(nxt).cloned().unwrap_or_default();
        neigh.retain(|n| *n < graph.n_blocks * 2);
        if neigh.is_empty() {
            let mut last_stack = Vec::new();
            let mut nback = 0usize;
            while last_stack.is_empty() && !stacks.is_empty() {
                last_stack = stacks.pop().unwrap();
                nback += 1;
            }
            if last_stack.is_empty() {
                return Ok(paths);
            }
            path.truncate(path.len() - 2 * nback);
            cur = last_stack.pop().unwrap();
            stacks.push(last_stack);
        } else {
            stacks.push(neigh);
            let frame = stacks.last_mut().unwrap();
            cur = frame.pop().unwrap();
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn reference_and_deletion_paths() {
        let mut g = Graph::new(3);
        g.add_edge(1, 2);
        g.add_edge(3, 4);
        g.add_edge(1, 4);
        let paths = enumerate_paths(&g, 100).unwrap();
        assert!(paths.iter().any(|p| p == &vec![0, 1, 2, 3, 4, 5]));
        assert!(paths.iter().any(|p| p == &vec![0, 1, 4, 5]));
    }

    #[test]
    fn too_many_paths_are_skipped() {
        let mut g = Graph::new(3);
        g.add_edge(1, 2);
        g.add_edge(3, 4);
        g.add_edge(1, 4);
        let err = enumerate_paths(&g, 0).unwrap_err();
        assert_eq!(err, "skipped_npaths");
    }
}
