use anyhow;

pub(crate) fn likely_valid_permit_list(
    num_unmatched: usize,
    total_mapped: usize,
    thresh: f64,
) -> anyhow::Result<f64> {
    if total_mapped > 0 {
        let unmatched_frac = (num_unmatched as f64) / (total_mapped as f64);
        if unmatched_frac < thresh {
            anyhow::Ok(unmatched_frac)
        } else {
            anyhow::bail!(
                "Percentage of mapped reads not matching a known barcode exaclty ({}%) is > the suggested fraction ({}%)",
                unmatched_frac * 100.0f64,
                thresh * 100.0f64
            )
        }
    } else {
        anyhow::bail!("Cannot determine (likely) valid permit list if not reads are mapped")
    }
}
