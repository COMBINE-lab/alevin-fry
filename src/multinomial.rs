use anyhow;
use rand_distr::weighted::{AliasableWeight, WeightedAliasIndex};
use rand_distr::Distribution;
use rand::Rng;
use nalgebra::base::OVector;
use nalgebra::base::dimension::Dyn;

pub struct Multinomial<W: AliasableWeight> {
    alias_idx: WeightedAliasIndex<W>,
    nsamp: usize
}

impl<W: AliasableWeight + std::fmt::Debug> Multinomial<W> {
    pub fn new(w: Vec<W>, n: usize) -> anyhow::Result<Self> {
        let idx = WeightedAliasIndex::new(w)?;
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

