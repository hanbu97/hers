use rand_core::{RngCore, SeedableRng};

// 保持原有的 FixedRandom 不变

#[derive(Debug, Clone)]
pub struct FixedRandom1 {
    initial_value: u64,
    current_value: u64,
    index: usize,
}

impl FixedRandom1 {
    pub fn new(value: u64) -> Self {
        FixedRandom1 {
            initial_value: value,
            current_value: value,
            index: 0,
        }
    }
}

impl RngCore for FixedRandom1 {
    fn next_u32(&mut self) -> u32 {
        self.next_u64() as u32
    }

    fn next_u64(&mut self) -> u64 {
        if self.index == 0 {
            self.index = 1;
            self.initial_value
        } else {
            self.current_value = self.current_value.wrapping_add(1);
            self.current_value
        }
    }

    fn fill_bytes(&mut self, dest: &mut [u8]) {
        for chunk in dest.chunks_mut(8) {
            let value = self.next_u64().to_le_bytes();
            chunk.copy_from_slice(&value[..chunk.len()]);
        }
    }

    fn try_fill_bytes(&mut self, dest: &mut [u8]) -> Result<(), rand_core::Error> {
        self.fill_bytes(dest);
        Ok(())
    }
}

impl SeedableRng for FixedRandom1 {
    type Seed = [u8; 8];

    fn from_seed(seed: Self::Seed) -> Self {
        let value = u64::from_le_bytes(seed);
        FixedRandom1::new(value)
    }
}
