use rand_core::{RngCore, SeedableRng};

#[derive(Debug, Clone)]
pub struct FixedRandom {
    values: Vec<u64>,
    index: usize,
    pub count: usize,
}

impl FixedRandom {
    pub fn new(values: Vec<u64>) -> Self {
        FixedRandom {
            values,
            index: 0,
            count: 0,
        }
    }
}

impl RngCore for FixedRandom {
    fn next_u32(&mut self) -> u32 {
        self.next_u64() as u32
    }

    fn next_u64(&mut self) -> u64 {
        let value = self.values[self.index];
        self.index = (self.index + 1) % self.values.len();
        value
    }

    fn fill_bytes(&mut self, dest: &mut [u8]) {
        self.count += 1;
        for chunk in dest.chunks_mut(8) {
            let value = self.next_u64().to_le_bytes();
            chunk.copy_from_slice(&value[..chunk.len()]);
        }

        println!(
            "------------------------- count: {} ------------------------- ",
            self.count
        );
    }

    fn try_fill_bytes(&mut self, dest: &mut [u8]) -> Result<(), rand_core::Error> {
        self.fill_bytes(dest);
        Ok(())
    }
}

impl SeedableRng for FixedRandom {
    type Seed = [u8; 0];

    fn from_seed(_seed: Self::Seed) -> Self {
        FixedRandom {
            values: Vec::new(),
            index: 0,
            count: 0,
        }
    }
}
