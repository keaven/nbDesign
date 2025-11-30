#' Find calendar date for target event count
#'
#' Finds the calendar time (since start of randomization) at which a specified
#' total number of events is reached in the simulated dataset.
#'
#' @param data A data frame of simulated data, typically from [nb_sim()].
#' @param planned_events Integer. The target number of events.
#'
#' @return Numeric. The calendar date when `planned_events` is achieved.
#'   If the dataset contains fewer than `planned_events`, returns the maximum
#'   calendar time in the dataset and prints a message.
#'
#' @examples
#' enroll_rate <- data.frame(rate = 20 / (5 / 12), duration = 5 / 12)
#' fail_rate <- data.frame(treatment = c("Control", "Experimental"), rate = c(0.5, 0.3))
#' dropout_rate <- data.frame(
#'   treatment = c("Control", "Experimental"),
#'   rate = c(0.1, 0.05), duration = c(100, 100)
#' )
#' sim <- nb_sim(enroll_rate, fail_rate, dropout_rate, max_followup = 2, n = 40)
#' get_analysis_date(sim, planned_events = 15)
#' @export
get_analysis_date <- function(data, planned_events) {
  # Filter for actual events (event == 1)
  events_df <- data[data$event == 1, ]

  total_events <- nrow(events_df)

  if (total_events < planned_events) {
    message(sprintf("Only %d events in trial", total_events))
    # Return max calendar time in the entire dataset (including censoring times)
    return(max(data$calendar_time))
  }

  # Sort events by calendar time
  events_df <- events_df[order(events_df$calendar_time), ]

  # The calendar time of the planned_events-th event
  return(events_df$calendar_time[planned_events])
}
