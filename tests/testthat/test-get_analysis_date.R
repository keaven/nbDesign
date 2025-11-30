test_that("get_analysis_date returns correct time", {
  # Simple simulation
  enroll_rate <- data.frame(rate = 10, duration = 1)
  fail_rate <- data.frame(treatment = "A", rate = 5)
  sim <- nb_sim(enroll_rate, fail_rate, max_followup = 5, n = 10)

  # Check if we have events
  total_events <- sum(sim$event == 1)
  skip_if(total_events < 5, "Not enough events generated for test")

  # Calculate manually
  events_df <- sim[sim$event == 1, ]
  events_df <- events_df[order(events_df$calendar_time), ]
  target_time <- events_df$calendar_time[5]

  res <- get_analysis_date(sim, planned_events = 5)
  expect_equal(res, target_time)
})

test_that("get_analysis_date handles insufficient events", {
  enroll_rate <- data.frame(rate = 10, duration = 1)
  fail_rate <- data.frame(treatment = "A", rate = 0.1) # Low rate
  sim <- nb_sim(enroll_rate, fail_rate, max_followup = 0.1, n = 10)

  # Likely very few events
  total_events <- sum(sim$event == 1)
  target <- total_events + 10

  expect_message(
    res <- get_analysis_date(sim, planned_events = target),
    sprintf("Only %d events in trial", total_events)
  )

  expect_equal(res, max(sim$calendar_time))
})
