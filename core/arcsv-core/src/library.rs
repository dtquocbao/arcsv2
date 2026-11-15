//! Insert pmf, CDF sums, and the O(1) rectangle normalizer.
//! orig:bamparser_streaming.py::pmf_kernel_smooth
//! orig:sv_inference_insertions.py::double_int

/// Gaussian KDE on the integer grid. Bandwidth is Silverman's rule as used by
/// ARC-SV: `max(1, 0.785 * IQR / n^(1/5))`.
pub fn pmf_kernel_smooth(inserts: &[i64], xmin: i64, xmax: i64) -> Vec<f64> {
    let truncated: Vec<i64> = inserts
        .iter()
        .copied()
        .filter(|x| *x >= xmin && *x <= xmax)
        .collect();
    assert!(!truncated.is_empty(), "no insert lengths");
    let mut sorted = truncated.clone();
    sorted.sort();
    let q1 = percentile(&sorted, 0.25);
    let q3 = percentile(&sorted, 0.75);
    let iqr = q3 - q1;
    let n = truncated.len() as f64;
    let bw = (0.785 * iqr / n.powf(0.2)).max(1.0);
    let width = (xmax - xmin + 1) as usize;
    let mut hist = vec![0usize; width];
    for x in &truncated {
        hist[(x - xmin) as usize] += 1;
    }
    let radius = (8.0 * bw).ceil().max(1.0) as i64;
    let inv = 1.0 / (bw * (2.0 * std::f64::consts::PI).sqrt() * n);
    let mut kernel = Vec::with_capacity((2 * radius + 1) as usize);
    for k in -radius..=radius {
        let z = k as f64 / bw;
        kernel.push((-0.5 * z * z).exp());
    }
    let mut pmf = vec![0.0; width];
    for (i, count) in hist.iter().enumerate() {
        if *count == 0 {
            continue;
        }
        let weight = *count as f64 * inv;
        let lo = (i as i64 - radius).max(0) as usize;
        let hi = ((i as i64 + radius + 1) as usize).min(width);
        for (j, slot) in pmf.iter_mut().enumerate().take(hi).skip(lo) {
            let k = (j as i64 - i as i64 + radius) as usize;
            *slot += weight * kernel[k];
        }
    }
    let total: f64 = pmf.iter().sum();
    pmf.iter().map(|p| p / total).collect()
}

pub fn percentile(sorted: &[i64], q: f64) -> f64 {
    if sorted.len() == 1 {
        return sorted[0] as f64;
    }
    let rank = (sorted.len() - 1) as f64 * q;
    let lo = rank.floor() as usize;
    let hi = rank.ceil() as usize;
    if lo == hi {
        sorted[lo] as f64
    } else {
        let w = rank - lo as f64;
        sorted[lo] as f64 * (1.0 - w) + sorted[hi] as f64 * w
    }
}

pub fn cdf_sum(pmf: &[f64]) -> Vec<f64> {
    let mut cdf = Vec::with_capacity(pmf.len());
    let mut acc = 0.0;
    for p in pmf {
        acc += p;
        cdf.push(acc);
    }
    let mut out = Vec::with_capacity(cdf.len());
    acc = 0.0;
    for f in cdf {
        acc += f;
        out.push(acc);
    }
    out
}

pub fn eval_cdf_sum(cdf_sum: &[f64], t: i64) -> f64 {
    if t < 0 {
        0.0
    } else if t as usize >= cdf_sum.len() {
        let n = cdf_sum.len() as i64;
        cdf_sum[cdf_sum.len() - 1] + (t - n + 1) as f64
    } else {
        cdf_sum[t as usize]
    }
}

/// orig:sv_inference_insertions.py::double_int
pub fn double_int(cdf_sum: &[f64], a1: i64, b1: i64, a2: i64, b2: i64, xmax: i64, ymax: i64) -> f64 {
    if a1 > xmax || a2 > ymax || b1 < 0 || b2 < 0 {
        return 0.0;
    }
    let a1 = a1.max(0);
    let a2 = a2.max(0);
    let b1 = b1.min(xmax);
    let b2 = b2.min(ymax);
    if a1 > b1 || a2 > b2 {
        return 0.0;
    }
    eval_cdf_sum(cdf_sum, b2 - a1 + 1)
        + eval_cdf_sum(cdf_sum, a2 - b1 - 1)
        - eval_cdf_sum(cdf_sum, b2 - b1)
        - eval_cdf_sum(cdf_sum, a2 - a1)
}

/// orig:call_sv.py::compute_pi_robust
pub fn pi_robust(pmf: &[f64], p: f64) -> f64 {
    let mut ordered = pmf.to_vec();
    ordered.sort_by(|a, b| b.partial_cmp(a).unwrap());
    let mut acc = 0.0;
    for value in ordered.iter().copied() {
        acc += value;
        if acc >= 1.0 - p {
            return value;
        }
    }
    ordered[ordered.len() - 1]
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn pmf_is_normalized() {
        let inserts = vec![400i64; 80];
        let pmf = pmf_kernel_smooth(&inserts, 0, 1200);
        let sum: f64 = pmf.iter().sum();
        assert!((sum - 1.0).abs() < 1e-12);
        assert!(pi_robust(&pmf, 1e-4) > 0.0);
    }

    #[test]
    fn double_int_matches_direct_sum() {
        let pmf = vec![0.0, 0.2, 0.5, 0.3];
        let cs = cdf_sum(&pmf);
        let mut brute = 0.0;
        for x in 0..=5 {
            for y in 0..=5 {
                let d = y - x + 1;
                if (0..pmf.len() as i64).contains(&d) {
                    brute += pmf[d as usize];
                }
            }
        }
        let got = double_int(&cs, 0, 5, 0, 5, 5, 5);
        assert!((got - brute).abs() < 1e-12);
    }
}
