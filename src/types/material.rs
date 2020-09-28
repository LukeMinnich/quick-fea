#[allow(non_snake_case)]
pub struct IsotropicMaterial {
    pub E: f64,
    pub G: f64,
    pub nu: f64,
}

#[allow(non_snake_case)]
impl IsotropicMaterial {
    pub fn new(E: f64, nu: f64) -> IsotropicMaterial {
        let G = E / (2. * (1. + nu));
        IsotropicMaterial { E, G, nu }
    }
}
