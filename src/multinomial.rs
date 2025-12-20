use anyhow;
use rand_distr::weighted::{Weight, WeightedIndex};
use rand_distr::Distribution;
use rand_distr::uniform::{SampleUniform, SampleBorrow};
use rand::Rng;
use nalgebra::base::OVector;
use nalgebra::base::dimension::Dyn;

// TODO: check how many samples are required before WeightedAliasIndex outperforms WeightedIndex

pub struct Multinomial<X: SampleUniform + PartialOrd + std::fmt::Debug> {
    alias_idx: WeightedIndex<X>,
    nsamp: usize
}

impl<X: SampleUniform + PartialOrd + std::fmt::Debug> Multinomial<X> {
    pub fn new<I>(w: I, n: usize) -> anyhow::Result<Self> 
    where
        I: IntoIterator,
        <I as IntoIterator>::Item: SampleBorrow<X>,
        X: Weight,
    {
        let idx = WeightedIndex::new(w)?;
        Ok(Self {
            alias_idx: idx,
            nsamp:n
        })
    }

    pub fn sample_u64<R: Rng + ?Sized>(&mut self, rnd: &mut R) -> OVector<u64, Dyn> {
        let mut v: Vec<u64> = Vec::with_capacity(self.nsamp);
        for _ in 0..self.nsamp {
            let idx = self.alias_idx.sample(rnd);
            v[idx] += 1;
        }
        OVector::<u64, Dyn>::from_vec(v)
    }

    pub fn sample_u32<R: Rng + ?Sized>(&mut self, rnd: &mut R) -> OVector<u32, Dyn> {
        let mut v: Vec<u32> = Vec::with_capacity(self.nsamp);
        for _ in 0..self.nsamp {
            let idx = self.alias_idx.sample(rnd);
            v[idx] += 1;
        }
        OVector::<u32, Dyn>::from_vec(v)
    }
}

