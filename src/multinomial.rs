use anyhow;
use rand::Rng;
use rand_distr::Distribution;
use rand_distr::uniform::{SampleBorrow, SampleUniform};
use rand_distr::weighted::{Weight, WeightedIndex};

// TODO: check how many samples are required before WeightedAliasIndex outperforms WeightedIndex

pub struct Multinomial<X: SampleUniform + PartialOrd + std::fmt::Debug> {
    alias_idx: WeightedIndex<X>,
    num_classes: usize,
    nsamp: usize,
}

impl<X: SampleUniform + PartialOrd + std::fmt::Debug> Multinomial<X> {
    pub fn new<I>(w: I, n: usize) -> anyhow::Result<Self>
    where
        I: IntoIterator,
        <I as IntoIterator>::Item: SampleBorrow<X>,
        X: Weight,
    {
        let weights: Vec<I::Item> = w.into_iter().collect();
        let num_classes = weights.len();
        let idx = WeightedIndex::new(weights)?;
        Ok(Self {
            alias_idx: idx,
            num_classes,
            nsamp: n,
        })
    }

    pub fn sample_u64<R: Rng + ?Sized>(&mut self, rnd: &mut R) -> Vec<u64> {
        let mut v = vec![0u64; self.num_classes];
        for _ in 0..self.nsamp {
            let idx = self.alias_idx.sample(rnd);
            v[idx] += 1;
        }
        v
    }

    pub fn sample_u32<R: Rng + ?Sized>(&mut self, rnd: &mut R) -> Vec<u32> {
        let mut v = vec![0u32; self.num_classes];
        for _ in 0..self.nsamp {
            let idx = self.alias_idx.sample(rnd);
            v[idx] += 1;
        }
        v
    }
}
