use std::sync::{Arc, Condvar, Mutex};

#[derive(Clone)]
/// Coordinates bounded admission of sequence buffers across workers.
pub struct MemoryBudget {
    inner: Arc<Inner>,
}

struct Inner {
    limit: u64,
    state: Mutex<State>,
    changed: Condvar,
}

#[derive(Default)]
struct State {
    used: u64,
    peak: u64,
    oversized: u64,
}

/// Releases an admitted allocation when dropped.
pub struct MemoryPermit {
    budget: MemoryBudget,
    bytes: u64,
}

#[derive(Clone, Copy)]
/// Summary of tracked memory use during an audit.
pub struct MemoryTelemetry {
    pub limit_bytes: u64,
    pub peak_tracked_bytes: u64,
    pub oversized_contigs: u64,
}

impl MemoryBudget {
    /// Creates a shared budget with the given byte limit.
    pub fn new(limit_bytes: u64) -> Self {
        Self {
            inner: Arc::new(Inner {
                limit: limit_bytes,
                state: Mutex::new(State::default()),
                changed: Condvar::new(),
            }),
        }
    }

    /// Waits for capacity, admitting an oversized contig only when alone.
    pub fn acquire(&self, bytes: u64) -> Result<MemoryPermit, String> {
        self.acquire_inner(bytes, true)?
            .ok_or("memory budget rejected allocation".into())
    }

    /// Waits for capacity without admitting an allocation larger than the limit.
    pub fn acquire_within_limit(&self, bytes: u64) -> Result<Option<MemoryPermit>, String> {
        self.acquire_inner(bytes, false)
    }

    fn acquire_inner(
        &self,
        bytes: u64,
        allow_oversized: bool,
    ) -> Result<Option<MemoryPermit>, String> {
        if bytes > self.inner.limit && !allow_oversized {
            return Ok(None);
        }
        let mut state = self
            .inner
            .state
            .lock()
            .map_err(|_| "memory budget poisoned")?;
        if bytes > self.inner.limit {
            while state.used != 0 {
                state = self
                    .inner
                    .changed
                    .wait(state)
                    .map_err(|_| "memory budget poisoned")?;
            }
            state.oversized += 1;
        } else {
            while state.used.saturating_add(bytes) > self.inner.limit {
                state = self
                    .inner
                    .changed
                    .wait(state)
                    .map_err(|_| "memory budget poisoned")?;
            }
        }
        state.used += bytes;
        state.peak = state.peak.max(state.used);
        Ok(Some(MemoryPermit {
            budget: self.clone(),
            bytes,
        }))
    }

    /// Returns the final limit, peak, and oversized-allocation counters.
    pub fn telemetry(&self) -> Result<MemoryTelemetry, String> {
        let state = self
            .inner
            .state
            .lock()
            .map_err(|_| "memory budget poisoned")?;
        Ok(MemoryTelemetry {
            limit_bytes: self.inner.limit,
            peak_tracked_bytes: state.peak,
            oversized_contigs: state.oversized,
        })
    }
}

impl Drop for MemoryPermit {
    fn drop(&mut self) {
        if let Ok(mut state) = self.budget.inner.state.lock() {
            state.used -= self.bytes;
            self.budget.inner.changed.notify_all();
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn tracks_peak_and_oversized_admissions() {
        let budget = MemoryBudget::new(10);
        let first = budget.acquire(6).unwrap();
        let second = budget.acquire(4).unwrap();
        drop((first, second));
        let oversized = budget.acquire(12).unwrap();
        drop(oversized);
        let telemetry = budget.telemetry().unwrap();
        assert_eq!(telemetry.peak_tracked_bytes, 12);
        assert_eq!(telemetry.oversized_contigs, 1);
    }

    #[test]
    fn rejects_oversized_alignment() {
        let budget = MemoryBudget::new(10);
        assert!(budget.acquire_within_limit(11).unwrap().is_none());
    }
}
