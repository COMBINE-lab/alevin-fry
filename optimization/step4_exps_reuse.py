"""Step 4: stop rebuilding the truncation model's lookup tables from scratch on every EM iteration.

`Exps::new` built, per iteration of a cell's EM, two tables of one entry per base up to the longest
transcript (the end-kernel density and its running sum) plus the length-law tables, each a fresh
allocation: a few megabytes allocated, zeroed and filled per iteration, for hundreds of iterations
per cell. The tables are now refreshed in place into buffers that live for the whole fit, with the
same arithmetic in the same order, so every value is bit-identical.
"""
import os
P = "/scratch1/zahra/source-code/alevin-fry/src/truncation.rs"
s = open(P).read()

def rep(old, new, count=1):
    global s
    n = s.count(old); assert n == count, f"expected {count} of {old[:70]!r}, found {n}"
    s = s.replace(old, new)

# GapLaw: fill the tables into existing buffers
rep("""    /// Per-nt density `dens[d]` and distribution function `cdf[x] = P(gap <= x)` for gaps up
    /// to `max_len`.
    fn tables(&self, max_len: u32) -> (Vec<f64>, Vec<f64>) {
        let dens: Vec<f64> = (0..=max_len)
            .map(|d| {
                let b = Self::bin(d);
                self.h.get(b).map_or(0.0, |p| p / Self::width(b))
            })
            .collect();
        let mut acc = 0.0;
        let cdf = dens
            .iter()
            .map(|x| {
                acc += x;
                acc
            })
            .collect();
        (dens, cdf)
    }
""",
"""    /// Per-nt density `dens[d]` and distribution function `cdf[x] = P(gap <= x)` for gaps up
    /// to `max_len`.
    fn tables(&self, max_len: u32) -> (Vec<f64>, Vec<f64>) {
        let (mut dens, mut cdf) = (Vec::new(), Vec::new());
        self.tables_into(max_len, &mut dens, &mut cdf);
        (dens, cdf)
    }

    /// [`GapLaw::tables`] written into existing buffers (the same values, no new allocation
    /// once the buffers have their capacity).
    fn tables_into(&self, max_len: u32, dens: &mut Vec<f64>, cdf: &mut Vec<f64>) {
        dens.clear();
        dens.extend((0..=max_len).map(|d| {
            let b = Self::bin(d);
            self.h.get(b).map_or(0.0, |p| p / Self::width(b))
        }));
        let mut acc = 0.0;
        cdf.clear();
        cdf.extend(dens.iter().map(|x| {
            acc += x;
            acc
        }));
    }
""")

# Exps: a refresh that reuses the buffers
rep("""impl Exps {
    fn new(par: &Params, max_len: u32) -> Self {
        let tabs = |laws: &[LenLaw]| laws.iter().map(LenTab::new).collect();
        let (g3d, g3c) = par.k3.abs.tables(max_len);
        let (g5d, g5c) = par.k5.abs.tables(max_len);
        Exps {
            l3: tabs(&par.len3),
            l5: tabs(&par.len5),
            lu: tabs(&par.lenu),
            g3d,
            g3c,
            g5d,
            g5c,
        }
    }
}
""",
"""impl Exps {
    fn new(par: &Params, max_len: u32) -> Self {
        let tabs = |laws: &[LenLaw]| laws.iter().map(LenTab::new).collect();
        let (g3d, g3c) = par.k3.abs.tables(max_len);
        let (g5d, g5c) = par.k5.abs.tables(max_len);
        Exps {
            l3: tabs(&par.len3),
            l5: tabs(&par.len5),
            lu: tabs(&par.lenu),
            g3d,
            g3c,
            g5d,
            g5c,
        }
    }

    /// The tables for new parameters, rebuilt in place: the end-kernel tables keep their
    /// buffers (one entry per base up to the longest transcript, rebuilt every iteration of a
    /// fit), and the small per-bin length tables are replaced.
    fn refresh(&mut self, par: &Params, max_len: u32) {
        let tabs = |laws: &[LenLaw], out: &mut Vec<LenTab>| {
            out.clear();
            out.extend(laws.iter().map(LenTab::new));
        };
        tabs(&par.len3, &mut self.l3);
        tabs(&par.len5, &mut self.l5);
        tabs(&par.lenu, &mut self.lu);
        par.k3.abs.tables_into(max_len, &mut self.g3d, &mut self.g3c);
        par.k5.abs.tables_into(max_len, &mut self.g5d, &mut self.g5c);
    }
}
""")

# the per-cell fit: one Exps for the whole fit
rep("""    let mut has_positional_evidence =
        global.is_some() && cell.alns.iter().any(|a| a.span.is_some());

    while niter < MAX_EM_ITER {""",
"""    let mut has_positional_evidence =
        global.is_some() && cell.alns.iter().any(|a| a.span.is_some());
    // the likelihood tables, refreshed in place whenever the parameters change
    let mut exps = Exps::new(&par, max_len);

    while niter < MAX_EM_ITER {""")
rep("""            None => {
                let e = Exps::new(&par, max_len);
                let st = step_with(&cell, &par, &e, &prev, &mut curr);
                steps += 1;
                has_positional_evidence |= st.n > 0.0;
                let ll_gain = st.ll - last_ll;
                last_ll = st.ll;
                let next = st.params(&par);
                let reason = freeze.update(moved(&next, &par), ll_gain);
                par = next;
                if reason.is_some() {
                    let e = Exps::new(&par, max_len);
                    frozen = Some(Fixed::new(&cell, &par, &e));""",
"""            None => {
                exps.refresh(&par, max_len);
                let st = step_with(&cell, &par, &exps, &prev, &mut curr);
                steps += 1;
                has_positional_evidence |= st.n > 0.0;
                let ll_gain = st.ll - last_ll;
                last_ll = st.ll;
                let next = st.params(&par);
                let reason = freeze.update(moved(&next, &par), ll_gain);
                par = next;
                if reason.is_some() {
                    exps.refresh(&par, max_len);
                    frozen = Some(Fixed::new(&cell, &par, &exps));""")

# the global fits: same per-iteration pattern
rep("""    while iterations < MAX_EM_ITER as usize {
        let e = Exps::new(&par, max_len);
        let (par_ref, e_ref) = (&par, &e);""",
"""    let mut exps = Exps::new(&par, max_len);
    while iterations < MAX_EM_ITER as usize {
        exps.refresh(&par, max_len);
        let e = &exps;
        let (par_ref, e_ref) = (&par, e);""")
rep("""    while iterations < MAX_EM_ITER as usize {
        let e = Exps::new(&par, max_len);
        let (st, next_alpha) = pooled_step(cells, &par, &e, &alpha, n_threads);""",
"""    let mut exps = Exps::new(&par, max_len);
    while iterations < MAX_EM_ITER as usize {
        exps.refresh(&par, max_len);
        let (st, next_alpha) = pooled_step(cells, &par, &exps, &alpha, n_threads);""")
open(P, "w").write(s)
print("step 4: truncation tables reused across iterations")
