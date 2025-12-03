
test_that("sample_size_nbinom adjusts for event_gap", {
  # No gap
  res_no_gap <- sample_size_nbinom(lambda1 = 0.5, lambda2 = 0.5, dispersion = 0.1, power = 0.8, exposure = 1)
  # With gap, effective rate is lower. 
  # Effective lambda = lambda / (1 + lambda * gap)
  # If gap is huge, effective rate -> 1/gap.
  
  # Let's check mu calculation indirectly via sample size or events.
  # For same sample size, power should drop if we add gap (less information).
  # Or for fixed power, sample size should increase.
  
  # Same parameters, gap = 0.5 (huge gap)
  # effective lambda = 0.5 / (1 + 0.5*0.5) = 0.5 / 1.25 = 0.4
  res_gap <- sample_size_nbinom(lambda1 = 0.5, lambda2 = 0.5, dispersion = 0.1, power = 0.8, exposure = 1, event_gap = 0.5)
  
  # Since lambdas are same, this is just a size check? 
  # Wait, if lambdas are same, N is infinite/undefined or just based on error?
  # sample_size_nbinom requires different lambdas for power?
  # The function checks inputs.
  
  # Let's use different lambdas
  l1 <- 0.5
  l2 <- 0.3
  gap <- 0.1
  
  res1 <- sample_size_nbinom(lambda1 = l1, lambda2 = l2, dispersion = 0, power = 0.8, exposure = 1)
  res2 <- sample_size_nbinom(lambda1 = l1, lambda2 = l2, dispersion = 0, power = 0.8, exposure = 1, event_gap = gap)
  
  expect_true(res2$n_total > res1$n_total)
  
  # specific check
  # effective l1 = 0.5 / (1 + 0.05) = 0.5/1.05 = 0.476
  # effective l2 = 0.3 / (1 + 0.03) = 0.3/1.03 = 0.291
  
  # res3 calculates sample size assuming the RATES themselves are lower (eff_l1, eff_l2).
  # This implies a smaller effect size (log(eff_l1/eff_l2) < log(l1/l2)), 
  # which would result in an even larger sample size than res2.
  # res2 maintains the original effect size (log(l1/l2)) because the analysis 
  # is assumed to correctly recover the rates using net-time offsets, 
  # but accounts for higher variance due to fewer events.
  
  eff_l1 <- l1 / (1 + l1 * gap)
  eff_l2 <- l2 / (1 + l2 * gap)
  
  res3 <- sample_size_nbinom(lambda1 = eff_l1, lambda2 = eff_l2, dispersion = 0, power = 0.8, exposure = 1)
  
  # res2 should be smaller than res3 because res2 has a stronger signal (larger effect size denominator)
  expect_true(res2$n_total < res3$n_total)
  
  # However, the number of events in res2 should be consistent with the effective rates
  # Variance calculation involves 1/mu. mu in res2 is calculated using eff_l.
  # So the "Variance" term in the formula (V) is identical for res2 and res3.
  # The difference is only in the Denominator (effect size).
  
  # Let's verify V is the same? We can't access V directly but we can check variance output
  # variance = V / n_approx. 
  # Since n differs, variance output differs.
  
  # We can check that the calculated events match manually
  mu1_eff <- eff_l1 * 1
  expect_equal(res2$events_n1, res2$n1 * mu1_eff)
})

