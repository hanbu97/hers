/// CRed reduce returns a mod q where a is between 0 and 2*q-1.
pub fn c_red(a: u64, q: u64) -> u64 {
    if a >= q {
        a - q
    } else {
        a
    }
}
