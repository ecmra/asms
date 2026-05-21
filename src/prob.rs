/// log of mu^x(1-mu)^(1-x)
pub fn log_bernoulli_term(x: f32, mu: f32) -> f32 {
    if (0.0 == x && 0.0 == mu) || (1.0 == x && 1.0 == mu) {
        return f32::ln(1.0);
    } else {
        return x * f32::ln(mu) + (1.0 - x) * f32::ln_1p(-mu);
    }
}

/// likelihood of a string of bits given a multivariate mean
pub fn logmvbern(x: &[f32], mu: &[f32]) -> f32 {
    let d = x.len();
    let mut logp: f32 = 0.0;
    for j in 0..d {
        if -1.0 == x[j] {
            continue;
        }
        logp += log_bernoulli_term(x[j], mu[j]);
        //println!("logp:{} x={} mu={}", log_bernoulli_term(x[j], mu[j]) , x[j], mu[j]);
    }
    //println!("logmvbern:{} x={:?} mu={:?}", logp, x, mu);
    return logp;
}

/// logsumexp(x1,x2,x3) = log(exp(x1)+exp(x2)+exp(x3))
pub fn logsumexp(v: &[f32]) -> f32 {
    let mut m = f32::MIN;
    for i in 0..v.len() {
        if v[i] > m {
            m = v[i];
        }
    }
    let mut s = 0.0;
    for i in 0..v.len() {
        s += f32::exp(v[i] - m);
    }
    return m + f32::ln(s);
}


#[test]
fn logsumexp_test() {
    assert!(f32::abs(logsumexp(&[1., 2., 3.]) - 3.40760596444438) < 1e-6);
    /*
    In [15]: np.log(np.exp(1)+np.exp(2)+np.exp(3))
    Out[15]: 3.40760596444438
     */
}
#[test]
fn log_bernoulli_term_test() {
    assert!(f32::abs(log_bernoulli_term(1.0, 0.2) - -1.6094379124341003) < 1e-6);
}
#[test]
fn logmvbern_test() {
    let x = &[0.0, 1.0];
    let mu = &[0.2, 0.8];
    assert!(f32::abs(logmvbern(x, mu) - -0.446287) < 1e-6);
}
