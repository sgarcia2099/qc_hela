#!/usr/bin/env Rscript

suppressPackageStartupMessages({
	library(dplyr)
	library(readr)
	library(stringr)
	library(tidyr)
	library(purrr)
})

parse_cli_args <- function(args) {
	defaults <- list(
		input_dir = ".",
		output_csv = NULL,
		pattern = "MSMSSpectrumInfo\\.txt$",
		recursive = FALSE
	)

	if (length(args) > 0) {
		kv <- args[str_detect(args, "^--")]
		for (item in kv) {
			key_val <- str_split_fixed(str_remove(item, "^--"), "=", 2)
			key <- key_val[, 1]
			val <- key_val[, 2]

			if (key == "input-dir" && nzchar(val)) defaults$input_dir <- val
			if (key == "output-csv" && nzchar(val)) defaults$output_csv <- val
			if (key == "pattern" && nzchar(val)) defaults$pattern <- val
			if (key == "recursive" && nzchar(val)) defaults$recursive <- tolower(val) %in% c("1", "true", "yes")
		}
	}

	if (is.null(defaults$output_csv) || !nzchar(defaults$output_csv)) {
		defaults$output_csv <- file.path(defaults$input_dir, "msms_duty_cycle_metrics.csv")
	}

	defaults
}

extract_sample_id <- function(path) {
	file <- basename(path)
	base <- str_remove(file, "\\.txt$")
	base <- str_remove(base, "_MSMSSpectrumInfo$")

	# Normalize common pre/post labels while preserving sample prefix.
	base <- base |>
		str_replace("_preRun$", "_preRun") |>
		str_replace("_postRun$", "_postRun")

	base
}

safe_numeric <- function(x) {
	suppressWarnings(as.numeric(x))
}

compute_file_metrics <- function(path) {
	df <- read_tsv(path, show_col_types = FALSE, progress = FALSE)

	required <- c("First Scan", "RT in min", "Ion Inject Time in ms")
	missing <- setdiff(required, names(df))
	if (length(missing) > 0) {
		return(tibble(
			sample_id = extract_sample_id(path),
			file = basename(path),
			file_path = path,
			status = "missing_required_columns",
			error = paste(missing, collapse = ", "),
			total_spectra = NA_real_,
			acquisition_window_s = NA_real_,
			scan_rate_hz = NA_real_,
			total_ion_inject_ms = NA_real_,
			mean_ion_inject_ms = NA_real_,
			median_ion_inject_ms = NA_real_,
			mean_cycle_time_ms = NA_real_,
			median_cycle_time_ms = NA_real_,
			duty_cycle_pct = NA_real_
		))
	}

	df2 <- df |>
		transmute(
			scan = suppressWarnings(as.integer(`First Scan`)),
			rt_min = safe_numeric(`RT in min`),
			ion_inject_ms = safe_numeric(`Ion Inject Time in ms`)
		) |>
		filter(!is.na(rt_min), rt_min >= 0, !is.na(scan)) |>
		arrange(scan)

	if (nrow(df2) == 0) {
		return(tibble(
			sample_id = extract_sample_id(path),
			file = basename(path),
			file_path = path,
			status = "no_valid_rows",
			error = "No rows with valid First Scan and RT in min",
			total_spectra = 0,
			acquisition_window_s = NA_real_,
			scan_rate_hz = NA_real_,
			total_ion_inject_ms = NA_real_,
			mean_ion_inject_ms = NA_real_,
			median_ion_inject_ms = NA_real_,
			mean_cycle_time_ms = NA_real_,
			median_cycle_time_ms = NA_real_,
			duty_cycle_pct = NA_real_
		))
	}

	df3 <- df2 |>
		mutate(cycle_time_ms = c(diff(rt_min * 60) * 1000, NA_real_))

	total_spectra <- nrow(df3)
	min_rt <- min(df3$rt_min, na.rm = TRUE)
	max_rt <- max(df3$rt_min, na.rm = TRUE)
	acq_s <- (max_rt - min_rt) * 60

	total_ion_inject_ms <- sum(df3$ion_inject_ms, na.rm = TRUE)
	duty_cycle_pct <- if (is.finite(acq_s) && acq_s > 0) total_ion_inject_ms / (acq_s * 1000) * 100 else NA_real_
	scan_rate_hz <- if (is.finite(acq_s) && acq_s > 0) total_spectra / acq_s else NA_real_

	tibble(
		sample_id = extract_sample_id(path),
		file = basename(path),
		file_path = path,
		status = "ok",
		error = NA_character_,
		total_spectra = total_spectra,
		acquisition_window_s = round(acq_s, 3),
		scan_rate_hz = round(scan_rate_hz, 3),
		total_ion_inject_ms = round(total_ion_inject_ms, 3),
		mean_ion_inject_ms = round(mean(df3$ion_inject_ms, na.rm = TRUE), 3),
		median_ion_inject_ms = round(median(df3$ion_inject_ms, na.rm = TRUE), 3),
		mean_cycle_time_ms = round(mean(df3$cycle_time_ms, na.rm = TRUE), 3),
		median_cycle_time_ms = round(median(df3$cycle_time_ms, na.rm = TRUE), 3),
		duty_cycle_pct = round(duty_cycle_pct, 3)
	)
}

main <- function() {
	args <- parse_cli_args(commandArgs(trailingOnly = TRUE))

	if (!dir.exists(args$input_dir)) {
		stop(sprintf("Input directory does not exist: %s", args$input_dir))
	}

	files <- list.files(
		path = args$input_dir,
		pattern = args$pattern,
		full.names = TRUE,
		recursive = isTRUE(args$recursive)
	)

	if (length(files) == 0) {
		stop(sprintf("No files matched pattern '%s' in: %s", args$pattern, args$input_dir))
	}

	message(sprintf("Found %d MSMSSpectrumInfo file(s).", length(files)))

	# map_dfr is memory-efficient for many files and creates a single output table.
	metrics <- map_dfr(files, compute_file_metrics) |>
		arrange(sample_id, file)

	out_dir <- dirname(args$output_csv)
	dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
	write_csv(metrics, args$output_csv, na = "")

	ok_n <- sum(metrics$status == "ok", na.rm = TRUE)
	bad_n <- nrow(metrics) - ok_n
	message(sprintf("Wrote %d row(s) to: %s", nrow(metrics), args$output_csv))
	message(sprintf("Status summary: ok=%d, issues=%d", ok_n, bad_n))
}

if (identical(environment(), globalenv())) {
	main()
}
