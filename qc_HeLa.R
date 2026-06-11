##### Script for QC Analysis of Pre/Post HeLa Runs #####

# Date: 2026-06-09
# Author(s): Sarah Garcia, GenAI

suppressPackageStartupMessages({
	library(dplyr)
	library(ggplot2)
	library(readr)
	library(stringr)
	library(tidyr)
	library(purrr)
	library(tibble)
})

RUN_COLORS <- c(preRun = "#1f77b4", postRun = "#d62728")

qc_theme <- function() {
	theme_minimal(base_size = 12) +
		theme(
			plot.title = element_text(face = "bold"),
			legend.title = element_blank(),
			panel.grid.minor = element_blank()
		)
}

qc_fill_scale <- function() {
	scale_fill_manual(
		values = RUN_COLORS,
		breaks = c("preRun", "postRun"),
		labels = c("Pre Run", "Post Run")
	)
}

parse_cli_args <- function(args) {
	# Defaults are intentionally data-centric so outputs sit with the dataset.
	defaults <- list(
		input_dir = "/home/jkg/OneDrive/Documents/All Files - personal/School_UT/Work/SG-01-008_BlackSoldierFlyGutMicrobiome/qc_hela/data",
		output_dir = NULL,
		render_report = TRUE,
		write_preprocess_script = TRUE
	)

	if (length(args) > 0) {
		kv <- args[str_detect(args, "^--")]
		for (item in kv) {
			key_val <- str_split_fixed(str_remove(item, "^--"), "=", 2)
			key <- key_val[, 1]
			val <- key_val[, 2]

			if (key == "input-dir" && nzchar(val)) defaults$input_dir <- val
			if (key == "output-dir" && nzchar(val)) defaults$output_dir <- val
			if (key == "render-report" && nzchar(val)) defaults$render_report <- tolower(val) %in% c("1", "true", "yes")
			if (key == "write-preprocess-script" && nzchar(val)) defaults$write_preprocess_script <- tolower(val) %in% c("1", "true", "yes")
		}
	}

	if (is.null(defaults$output_dir) || !nzchar(defaults$output_dir)) {
		# By default, write output next to the input data directory.
		defaults$output_dir <- file.path(dirname(defaults$input_dir), "output")
	}

	defaults
}

stop_if_missing_columns <- function(df, cols, label) {
	missing <- setdiff(cols, names(df))
	if (length(missing) > 0) {
		stop(sprintf("%s missing required columns: %s", label, paste(missing, collapse = ", ")))
	}
}

sanitize_metric_name <- function(x) {
	x |>
		str_replace_all("[^A-Za-z0-9]+", "_") |>
		str_replace_all("_+", "_") |>
		str_replace_all("^_|_$", "") |>
		tolower()
}

normalize_sample_metric <- function(metric) {
	metric |>
		str_replace_all("_f[0-9]+_", "_sample_") |>
		str_replace_all("_f[0-9]+$", "_sample")
}

read_txt_table <- function(path) {
	read_tsv(path, show_col_types = FALSE, progress = FALSE)
}

discover_input_files <- function(input_dir) {
	if (!dir.exists(input_dir)) {
		stop(sprintf("Input directory does not exist: %s", input_dir))
	}

	files <- list.files(input_dir, pattern = "\\.txt$", full.names = TRUE)
	if (length(files) == 0) {
		stop(sprintf("No .txt files found in input directory: %s", input_dir))
	}

	expected_types <- c(
		"MSMSSpectrumInfo", "PSMs", "PeptideGroups",
		"ProteinGroups", "ResultStatistics", "StatisticsInsights"
	)

	index <- tibble(path = files) |>
		mutate(
			file = basename(path),
			sample = case_when(
				str_detect(file, "_preRun_") | str_detect(file, "_preRun") ~ "preRun",
				str_detect(file, "_postRun_") | str_detect(file, "_postRun") ~ "postRun",
				TRUE ~ NA_character_
			),
			file_type = str_match(file, "_(MSMSSpectrumInfo|PSMs|PeptideGroups|ProteinGroups|ResultStatistics|StatisticsInsights)\\.txt$")[, 2]
		) |>
		filter(!is.na(sample), !is.na(file_type))

	if (nrow(index) == 0) {
		stop("No preRun/postRun files matched expected naming patterns.")
	}

	expected_grid <- expand_grid(sample = c("preRun", "postRun"), file_type = expected_types)
	missing_pairs <- anti_join(expected_grid, index |> select(sample, file_type), by = c("sample", "file_type"))
	if (nrow(missing_pairs) > 0) {
		msg <- missing_pairs |>
			mutate(label = paste(sample, file_type, sep = "::")) |>
			pull(label) |>
			paste(collapse = ", ")
		stop(sprintf("Missing required file pair(s): %s", msg))
	}

	duplicate_pairs <- index |>
		count(sample, file_type, name = "n") |>
		filter(n > 1)
	if (nrow(duplicate_pairs) > 0) {
		msg <- duplicate_pairs |>
			mutate(label = paste(sample, file_type, sep = "::")) |>
			pull(label) |>
			paste(collapse = ", ")
		stop(sprintf("Duplicate files found for pair(s): %s", msg))
	}

	index
}

load_tables <- function(index) {
	split(index, interaction(index$sample, index$file_type, drop = TRUE)) |>
		imap(function(df, key) {
			list(
				sample = df$sample[[1]],
				file_type = df$file_type[[1]],
				path = df$path[[1]],
				data = read_txt_table(df$path[[1]])
			)
		})
}

get_table <- function(tbls, sample, file_type) {
	key <- paste(sample, file_type, sep = ".")
	tbls[[key]]$data
}

extract_statistics_insights <- function(tbls) {
	bind_rows(
		get_table(tbls, "preRun", "StatisticsInsights") |> mutate(run = "preRun"),
		get_table(tbls, "postRun", "StatisticsInsights") |> mutate(run = "postRun")
	) |>
		mutate(
			metric = sanitize_metric_name(Name),
			value = suppressWarnings(as.numeric(Value))
		) |>
		select(run, metric, value, Name, Units)
}

extract_result_statistics <- function(tbls) {
	keep_patterns <- c(
		"^MS/MS Spectrum Info -",
		"^PSMs - Precursor Abundance$",
		"^PSMs - q-Value$",
		"^PSMs - PEP$",
		"^PSMs - SVM Score$",
		"^PSMs - Charge$",
		"^PSMs - Isolation Interference"
	)

	pattern <- paste(keep_patterns, collapse = "|")

	bind_rows(
		get_table(tbls, "preRun", "ResultStatistics") |> mutate(run = "preRun"),
		get_table(tbls, "postRun", "ResultStatistics") |> mutate(run = "postRun")
	) |>
		filter(str_detect(Name, pattern), !str_detect(Name, "RT")) |>
		mutate(
			metric = sanitize_metric_name(Name),
			arithmetic_mean = suppressWarnings(as.numeric(`Arithmetic Mean`)),
			median = suppressWarnings(as.numeric(Median)),
			count = suppressWarnings(as.numeric(Count))
		) |>
		select(run, metric, arithmetic_mean, median, count, Name)
}

build_summary_metrics <- function(insights, result_stats) {
	insight_metrics <- insights |>
		group_by(run, metric) |>
		summarise(value = first(value), .groups = "drop") |>
		mutate(metric = normalize_sample_metric(metric))

	rs_metrics <- result_stats |>
		transmute(run, metric = normalize_sample_metric(paste0(metric, "_mean")), value = arithmetic_mean)

	merged <- bind_rows(insight_metrics, rs_metrics) |>
		filter(!is.na(value)) |>
		group_by(run, metric) |>
		summarise(value = mean(value, na.rm = TRUE), .groups = "drop")

	wide <- merged |>
		pivot_wider(names_from = run, values_from = value) |>
		mutate(
			delta_abs = postRun - preRun,
			delta_pct = if_else(abs(preRun) > 0, (postRun - preRun) / preRun * 100, NA_real_),
			decline_flag = case_when(
				is.na(delta_pct) ~ NA,
				delta_pct <= -20 ~ "major_decline",
				delta_pct <= -5 ~ "minor_decline",
				delta_pct < 5 ~ "stable",
				TRUE ~ "increase"
			)
		) |>
		arrange(metric)

	wide
}

safe_log2fc <- function(post, pre) {
	log2((post + 1) / (pre + 1))
}

analyze_shared_peptides <- function(tbls) {
	pre <- get_table(tbls, "preRun", "PeptideGroups")
	post <- get_table(tbls, "postRun", "PeptideGroups")

	required <- c("Annotated Sequence", "Modifications")
	stop_if_missing_columns(pre, required, "PeptideGroups preRun")
	stop_if_missing_columns(post, required, "PeptideGroups postRun")

	pre_abundance_col <- names(pre)[str_detect(names(pre), "^Abundance ")][1]
	post_abundance_col <- names(post)[str_detect(names(post), "^Abundance ")][1]

	if (is.na(pre_abundance_col) || is.na(post_abundance_col)) {
		stop("Could not locate abundance columns in PeptideGroups files.")
	}

	pre_df <- pre |>
		transmute(
			feature_key = paste(`Annotated Sequence`, Modifications, sep = "||"),
			abundance_pre = suppressWarnings(as.numeric(.data[[pre_abundance_col]]))
		) |>
		group_by(feature_key) |>
		summarise(
			abundance_pre = sum(abundance_pre, na.rm = TRUE),
			.groups = "drop"
		)

	post_df <- post |>
		transmute(
			feature_key = paste(`Annotated Sequence`, Modifications, sep = "||"),
			abundance_post = suppressWarnings(as.numeric(.data[[post_abundance_col]]))
		) |>
		group_by(feature_key) |>
		summarise(
			abundance_post = sum(abundance_post, na.rm = TRUE),
			.groups = "drop"
		)

	shared <- inner_join(pre_df, post_df, by = "feature_key") |>
		mutate(
			log2_fc_post_pre = safe_log2fc(abundance_post, abundance_pre),
			mean_log10_abundance = log10((abundance_pre + abundance_post) / 2 + 1)
		)

	overlap_summary <- tibble(
		metric = c(
			"peptides_pre_total",
			"peptides_post_total",
			"peptides_shared_total",
			"peptides_shared_vs_pre_pct",
			"peptides_shared_vs_post_pct",
			"peptides_shared_vs_union_pct"
		),
		value = c(
			nrow(pre_df),
			nrow(post_df),
			nrow(shared),
			ifelse(nrow(pre_df) > 0, nrow(shared) / nrow(pre_df) * 100, NA_real_),
			ifelse(nrow(post_df) > 0, nrow(shared) / nrow(post_df) * 100, NA_real_),
			ifelse((nrow(pre_df) + nrow(post_df) - nrow(shared)) > 0,
				nrow(shared) / (nrow(pre_df) + nrow(post_df) - nrow(shared)) * 100,
				NA_real_)
		)
	)

	list(shared = shared, overlap_summary = overlap_summary)
}

analyze_protein_groups <- function(tbls) {
	pre <- get_table(tbls, "preRun", "ProteinGroups")
	post <- get_table(tbls, "postRun", "ProteinGroups")

	key_col_pre <- names(pre)[str_detect(names(pre), "Protein Group ID|Group Description")][1]
	key_col_post <- names(post)[str_detect(names(post), "Protein Group ID|Group Description")][1]

	psm_col_pre <- names(pre)[str_detect(names(pre), "Number of PSMs")][1]
	psm_col_post <- names(post)[str_detect(names(post), "Number of PSMs")][1]

	if (any(is.na(c(key_col_pre, key_col_post, psm_col_pre, psm_col_post)))) {
		return(tibble())
	}

	pre_df <- pre |>
		transmute(
			protein_group_key = as.character(.data[[key_col_pre]]),
			psm_pre = suppressWarnings(as.numeric(.data[[psm_col_pre]]))
		) |>
		group_by(protein_group_key) |>
		summarise(psm_pre = sum(psm_pre, na.rm = TRUE), .groups = "drop")

	post_df <- post |>
		transmute(
			protein_group_key = as.character(.data[[key_col_post]]),
			psm_post = suppressWarnings(as.numeric(.data[[psm_col_post]]))
		) |>
		group_by(protein_group_key) |>
		summarise(psm_post = sum(psm_post, na.rm = TRUE), .groups = "drop")

	inner_join(pre_df, post_df, by = "protein_group_key") |>
		mutate(
			log2_fc_psm_post_pre = safe_log2fc(psm_post, psm_pre)
		)
}

build_psm_quality_table <- function(tbls) {
	# PSM table drives abundance/q-value trend summaries.
	pre <- get_table(tbls, "preRun", "PSMs")
	post <- get_table(tbls, "postRun", "PSMs")

	fetch_col <- function(nms, pattern) {
		hit <- nms[str_detect(nms, pattern)]
		if (length(hit) == 0) NA_character_ else hit[[1]]
	}

	map_psm <- function(df, run_label) {
		abundance_col <- fetch_col(names(df), "Precursor Abundance")
		q_col <- fetch_col(names(df), "^q-Value$")
		pep_col <- fetch_col(names(df), "^PEP$")

		tibble(
			run = run_label,
			precursor_abundance = if (!is.na(abundance_col)) suppressWarnings(as.numeric(df[[abundance_col]])) else NA_real_,
			q_value = if (!is.na(q_col)) suppressWarnings(as.numeric(df[[q_col]])) else NA_real_,
			pep = if (!is.na(pep_col)) suppressWarnings(as.numeric(df[[pep_col]])) else NA_real_
		)
	}

	bind_rows(map_psm(pre, "preRun"), map_psm(post, "postRun"))
}

build_msms_quality_table <- function(tbls) {
	# Isolation interference is reliably available in MS/MS spectrum-level exports.
	pre <- get_table(tbls, "preRun", "MSMSSpectrumInfo")
	post <- get_table(tbls, "postRun", "MSMSSpectrumInfo")

	fetch_col <- function(nms, pattern) {
		hit <- nms[str_detect(nms, pattern)]
		if (length(hit) == 0) NA_character_ else hit[[1]]
	}

	map_msms <- function(df, run_label) {
		interference_col <- fetch_col(names(df), "Isolation Interference")

		tibble(
			run = run_label,
			isolation_interference = if (!is.na(interference_col)) suppressWarnings(as.numeric(df[[interference_col]])) else NA_real_
		)
	}

	bind_rows(map_msms(pre, "preRun"), map_msms(post, "postRun"))
}

compute_duty_cycle <- function(tbls) {
	# Derive scan-rate (Hz) and cycle-time (ms) from MS/MS spectrum-level timestamps.
	process_run <- function(df, run_label) {
		required_cols <- c("RT in min", "First Scan", "Ion Inject Time in ms")
		missing_cols <- setdiff(required_cols, names(df))
		if (length(missing_cols) > 0) {
			stop(sprintf("MSMSSpectrumInfo (%s) missing columns: %s", run_label, paste(missing_cols, collapse = ", ")))
		}

		df |>
			transmute(
				run = run_label,
				scan = suppressWarnings(as.integer(`First Scan`)),
				rt_min = suppressWarnings(as.numeric(`RT in min`)),
				ion_inject_ms = suppressWarnings(as.numeric(`Ion Inject Time in ms`))
			) |>
			filter(!is.na(rt_min), !is.na(scan), rt_min >= 0) |>
			arrange(scan) |>
			mutate(
				# Cycle time: elapsed seconds to the next scan, expressed in ms
				cycle_time_ms = c(diff(rt_min * 60) * 1000, NA_real_)
			)
	}

	pre_scans  <- process_run(get_table(tbls, "preRun",  "MSMSSpectrumInfo"), "preRun")
	post_scans <- process_run(get_table(tbls, "postRun", "MSMSSpectrumInfo"), "postRun")
	scan_data  <- bind_rows(pre_scans, post_scans)

	# Run-level summary: total spectra, scan rate, mean/median cycle times
	summary <- scan_data |>
		group_by(run) |>
		summarise(
			total_spectra          = n(),
			acquisition_window_s   = (max(rt_min, na.rm = TRUE) - min(rt_min, na.rm = TRUE)) * 60,
			scan_rate_hz           = round(total_spectra / acquisition_window_s, 3),
			mean_ion_inject_ms     = round(mean(ion_inject_ms,  na.rm = TRUE), 2),
			median_ion_inject_ms   = round(median(ion_inject_ms, na.rm = TRUE), 2),
			mean_cycle_time_ms     = round(mean(cycle_time_ms,  na.rm = TRUE), 2),
			median_cycle_time_ms   = round(median(cycle_time_ms, na.rm = TRUE), 2),
			.groups = "drop"
		)

	# 1-minute RT bins for the scans-over-time line chart
	binned <- scan_data |>
		mutate(rt_bin_min = floor(rt_min)) |>
		group_by(run, rt_bin_min) |>
		summarise(
			scans_per_min      = n(),
			hz                 = round(scans_per_min / 60, 3),
			mean_ion_inject_ms = round(mean(ion_inject_ms, na.rm = TRUE), 2),
			.groups = "drop"
		)

	list(summary = summary, scan_data = scan_data, binned = binned)
}

compute_bias_metrics <- function(shared_peptides) {
	peptide_bias <- shared_peptides |>
		summarise(
			shared_features = n(),
			median_log2_fc = median(log2_fc_post_pre, na.rm = TRUE),
			frac_negative_fc = mean(log2_fc_post_pre < 0, na.rm = TRUE),
			frac_abs_fc_gt1 = mean(abs(log2_fc_post_pre) > 1, na.rm = TRUE),
			wilcox_p = tryCatch(wilcox.test(log2_fc_post_pre, mu = 0)$p.value, error = function(e) NA_real_)
		)

	list(peptide_bias = peptide_bias)
}

build_peptide_intensity_table <- function(tbls) {
	pre <- get_table(tbls, "preRun", "PeptideGroups")
	post <- get_table(tbls, "postRun", "PeptideGroups")

	fetch_abundance_col <- function(df) {
		names(df)[str_detect(names(df), "^Abundance ")][1]
	}

	pre_col <- fetch_abundance_col(pre)
	post_col <- fetch_abundance_col(post)

	if (is.na(pre_col) || is.na(post_col)) {
		stop("Could not locate peptide abundance columns for intensity histogram.")
	}

	bind_rows(
		tibble(run = "preRun", peptide_abundance = suppressWarnings(as.numeric(pre[[pre_col]]))),
		tibble(run = "postRun", peptide_abundance = suppressWarnings(as.numeric(post[[post_col]])))
	) |>
		filter(!is.na(peptide_abundance), peptide_abundance > 0) |>
		mutate(log10_peptide_abundance = log10(peptide_abundance + 1))
}

write_preprocess_bash_script <- function(output_dir, input_dir) {
	script_path <- file.path(output_dir, "preprocess_tables.sh")
	script <- c(
		"#!/usr/bin/env bash",
		"set -euo pipefail",
		sprintf("INPUT_DIR=%s", shQuote(input_dir)),
		"echo 'Compressing large .txt tables to .txt.gz for archival (originals kept)...'",
		"find \"$INPUT_DIR\" -maxdepth 1 -type f -name '*.txt' -size +10M -print0 | while IFS= read -r -d '' f; do",
		"  gzip -c \"$f\" > \"$f.gz\"",
		"done",
		"echo 'Done.'"
	)
	writeLines(script, script_path)
	Sys.chmod(script_path, mode = "0755")
	script_path
}

save_plots <- function(out_plot_dir, insights, shared_peptides, peptide_intensity, psm_quality, msms_quality, duty_cycle) {
	dir.create(out_plot_dir, recursive = TRUE, showWarnings = FALSE)
	# Remove prior plot files so stale legacy plots do not leak into the HTML report.
	unlink(list.files(out_plot_dir, pattern = "\\.png$", full.names = TRUE), force = TRUE)

	counts_plot <- insights |>
		filter(str_detect(metric, "high_confident")) |>
		mutate(metric_label = str_replace_all(metric, "_", " ")) |>
		ggplot(aes(x = metric_label, y = value, fill = run)) +
		geom_col(position = position_dodge(width = 0.7), width = 0.65) +
		coord_flip() +
		qc_fill_scale() +
		qc_theme() +
		labs(title = "High-Confidence IDs: Pre vs Post", x = NULL, y = "Count")
	ggsave(file.path(out_plot_dir, "01_high_confidence_counts.png"), counts_plot, width = 9, height = 5, dpi = 140)

	scatter <- ggplot(shared_peptides, aes(x = log10(abundance_pre + 1), y = log10(abundance_post + 1))) +
		geom_point(alpha = 0.35, size = 1.1, color = "#1f77b4") +
		geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "#d62728") +
		qc_theme() +
		labs(
			title = "Shared Peptides: Pre vs Post Abundance",
			x = "log10 Pre abundance",
			y = "log10 Post abundance"
		)
	ggsave(file.path(out_plot_dir, "02_shared_peptide_scatter.png"), scatter, width = 7, height = 6, dpi = 140)

	ma <- ggplot(shared_peptides, aes(x = mean_log10_abundance, y = log2_fc_post_pre)) +
		geom_point(alpha = 0.3, size = 1, color = "#2ca02c") +
		geom_hline(yintercept = 0, linetype = "dashed", color = "#d62728") +
		qc_theme() +
		labs(title = "MA Plot (Shared Peptides)", x = "Mean log10 abundance", y = "log2(Post/Pre)")
	ggsave(file.path(out_plot_dir, "03_ma_plot.png"), ma, width = 7, height = 6, dpi = 140)

	hist_fc <- ggplot(shared_peptides, aes(x = log2_fc_post_pre)) +
		geom_histogram(bins = 70, fill = "#4c78a8", color = "white", alpha = 0.9) +
		geom_vline(xintercept = 0, linetype = "dashed", color = "#d62728") +
		qc_theme() +
		labs(title = "Distribution of Shared-Peptide log2 Fold-Change", x = "log2(Post/Pre)", y = "Features")
	ggsave(file.path(out_plot_dir, "04_log2fc_histogram.png"), hist_fc, width = 8, height = 5.5, dpi = 140)

	peptide_intensity_hist <- peptide_intensity |>
		ggplot(aes(x = log10_peptide_abundance, fill = run)) +
		geom_histogram(alpha = 0.45, bins = 70, position = "identity") +
		qc_fill_scale() +
		qc_theme() +
		labs(
			title = "Peptide Intensity Distribution: Pre vs Post",
			x = "log10 peptide abundance",
			y = "Features"
		)
	ggsave(file.path(out_plot_dir, "05_peptide_intensity_histogram.png"), peptide_intensity_hist, width = 8, height = 5.5, dpi = 140)

	psm_abundance <- psm_quality |>
		filter(!is.na(precursor_abundance), precursor_abundance > 0) |>
		mutate(log10_precursor_abundance = log10(precursor_abundance + 1)) |>
		ggplot(aes(x = run, y = log10_precursor_abundance, fill = run)) +
		geom_boxplot(outlier.alpha = 0.15) +
		qc_fill_scale() +
		qc_theme() +
		labs(title = "PSM Precursor Abundance by Run", x = NULL, y = "log10 abundance")
	ggsave(file.path(out_plot_dir, "06_psm_precursor_abundance_boxplot.png"), psm_abundance, width = 7, height = 5, dpi = 140)

	interference <- msms_quality |>
		filter(!is.na(isolation_interference)) |>
		ggplot(aes(x = run, y = isolation_interference, fill = run)) +
		geom_violin(trim = TRUE, alpha = 0.85) +
		qc_fill_scale() +
		qc_theme() +
		labs(title = "Isolation Interference by Run", x = NULL, y = "Isolation Interference (%)")
	ggsave(file.path(out_plot_dir, "07_isolation_interference_violin.png"), interference, width = 7, height = 5, dpi = 140)

	# Build per-run mean Hz label for annotation
	hz_labels <- duty_cycle$summary |>
		transmute(
			run,
			label = sprintf("%s\nMean: %.2f Hz", run, scan_rate_hz),
			x_pos = max(duty_cycle$binned$rt_bin_min, na.rm = TRUE) * 0.98,
			y_pos = scan_rate_hz
		)

	# Plot 08: scans/sec over acquisition time (1-min bins), colored by run
	hz_over_time <- duty_cycle$binned |>
		ggplot(aes(x = rt_bin_min, y = hz, color = run)) +
		geom_line(linewidth = 0.9, alpha = 0.85) +
		geom_hline(
			data = duty_cycle$summary,
			aes(yintercept = scan_rate_hz, color = run),
			linetype = "dashed", linewidth = 0.55, alpha = 0.7
		) +
		geom_label(
			data = hz_labels,
			aes(x = x_pos, y = y_pos, label = label, color = run),
			hjust = 1, vjust = 0.5, size = 3.2, show.legend = FALSE,
			fill = "white", linewidth = 0.2
		) +
		scale_color_manual(values = RUN_COLORS, breaks = c("preRun", "postRun"),
						   labels = c("Pre Run", "Post Run")) +
		qc_theme() +
		labs(
			title  = "MS2 Scan Rate Over Acquisition Time",
			subtitle = "Dashed lines show run mean; labels show average Hz",
			x = "Retention time (min)",
			y = "Scan rate (Hz)",
			color  = NULL
		)
	ggsave(file.path(out_plot_dir, "08_duty_cycle_hz_over_time.png"), hz_over_time, width = 10, height = 5, dpi = 140)

	# Plot 09: side-by-side bars for scan rate and mean cycle time
	dc_long <- duty_cycle$summary |>
		select(run, scan_rate_hz, mean_cycle_time_ms) |>
		pivot_longer(c(scan_rate_hz, mean_cycle_time_ms),
					 names_to = "metric", values_to = "value") |>
		mutate(metric_label = if_else(
			metric == "scan_rate_hz", "Scan Rate (Hz)", "Mean Cycle Time (ms)"
		))

	dc_comparison <- dc_long |>
		ggplot(aes(x = run, y = value, fill = run)) +
		geom_col(width = 0.6) +
		geom_text(aes(label = round(value, 2)), vjust = -0.4, size = 3.5) +
		facet_wrap(~ metric_label, scales = "free_y") +
		qc_fill_scale() +
		qc_theme() +
		theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
		labs(
			title = "Duty Cycle: Scan Rate and Mean Cycle Time",
			x = NULL, y = NULL
		)
	ggsave(file.path(out_plot_dir, "09_duty_cycle_comparison.png"), dc_comparison, width = 8, height = 5, dpi = 140)
}

write_html_report <- function(output_dir, summary_metrics, peptide_bias, overlap_summary, duty_cycle_summary, plot_dir) {
	report_path <- file.path(output_dir, "qc_report.html")
	key <- peptide_bias |>
		mutate(across(everything(), ~ format(., digits = 4)))

	overlap_vals <- overlap_summary |>
		select(metric, value) |>
		deframe()

	pre_total <- ifelse(is.na(overlap_vals[["peptides_pre_total"]]), 0, overlap_vals[["peptides_pre_total"]])
	post_total <- ifelse(is.na(overlap_vals[["peptides_post_total"]]), 0, overlap_vals[["peptides_post_total"]])
	shared_total <- ifelse(is.na(overlap_vals[["peptides_shared_total"]]), 0, overlap_vals[["peptides_shared_total"]])
	shared_vs_pre <- overlap_vals[["peptides_shared_vs_pre_pct"]]
	shared_vs_post <- overlap_vals[["peptides_shared_vs_post_pct"]]
	shared_vs_union <- overlap_vals[["peptides_shared_vs_union_pct"]]

	table_rows <- summary_metrics |>
		mutate(
			flag_priority = case_when(
				decline_flag == "major_decline" ~ 1,
				decline_flag == "minor_decline" ~ 2,
				decline_flag == "stable" ~ 3,
				decline_flag == "increase" ~ 4,
				TRUE ~ 5
			)
		) |>
		filter(!is.na(decline_flag), decline_flag != "stable") |>
		arrange(flag_priority, desc(abs(delta_pct))) |>
		head(30) |>
		mutate(across(c(preRun, postRun, delta_abs, delta_pct), ~ format(., digits = 4))) |>
		transmute(row = sprintf(
			"<tr><td>%s</td><td>%s</td><td>%s</td><td>%s</td><td>%s</td><td>%s</td></tr>",
			metric, preRun, postRun, delta_abs, delta_pct, decline_flag
		)) |>
		pull(row)

	images <- list.files(plot_dir, pattern = "\\.png$", full.names = FALSE) |>
		sort()

	image_html <- paste0(
		"<div style='display:flex;flex-wrap:wrap;gap:18px;'>",
		paste0(
			sprintf(
				"<div style='width:47%%'><p><b>%s</b></p><img src='plots/%s' style='max-width:100%%;border:1px solid #ccc;'/></div>",
				images, images
			),
			collapse = ""
		),
		"</div>"
	)

	html <- c(
		"<!doctype html>",
		"<html><head><meta charset='utf-8'/><title>HeLa QC Report</title></head><body style='font-family:Helvetica,Arial,sans-serif;'>",
		"<h1>Pre/Post HeLa Run QC Report</h1>",
		sprintf("<p><b>Generated:</b> %s</p>", as.character(Sys.time())),
		"<h2>Run-order bias summary</h2>",
		sprintf("<p>Shared features: %s</p>", key$shared_features),
		sprintf("<p>Shared proportion vs preRun: %s / %s (%s%%)</p>",
			format(shared_total, big.mark = ",", scientific = FALSE),
			format(pre_total, big.mark = ",", scientific = FALSE),
			format(round(shared_vs_pre, 2), nsmall = 2)
		),
		sprintf("<p>Shared proportion vs postRun: %s / %s (%s%%)</p>",
			format(shared_total, big.mark = ",", scientific = FALSE),
			format(post_total, big.mark = ",", scientific = FALSE),
			format(round(shared_vs_post, 2), nsmall = 2)
		),
		sprintf("<p>Shared proportion vs union (Jaccard): %s%%</p>",
			format(round(shared_vs_union, 2), nsmall = 2)
		),
		sprintf("<p>Median log2(Post/Pre): %s</p>", key$median_log2_fc),
		sprintf("<p>Fraction with negative fold-change: %s</p>", key$frac_negative_fc),
		sprintf("<p>Wilcoxon p-value (log2 FC vs 0): %s</p>", key$wilcox_p),
		"<h2>Duty cycle</h2>",
		"<table border='1' cellpadding='6' cellspacing='0'>",
		"<tr><th>Run</th><th>Total Spectra</th><th>Scan Rate (Hz)</th><th>Mean Cycle Time (ms)</th><th>Median Cycle Time (ms)</th><th>Mean Ion Inject Time (ms)</th></tr>",
		paste0(apply(duty_cycle_summary, 1, function(r) sprintf(
			"<tr><td>%s</td><td>%s</td><td>%s</td><td>%s</td><td>%s</td><td>%s</td></tr>",
			r["run"], r["total_spectra"], r["scan_rate_hz"],
			r["mean_cycle_time_ms"], r["median_cycle_time_ms"], r["mean_ion_inject_ms"]
		)), collapse = ""),
		"</table>",
		"<h2>Top metric changes</h2>",
		"<table border='1' cellpadding='6' cellspacing='0'>",
		"<tr><th>Metric</th><th>Pre</th><th>Post</th><th>Delta</th><th>Delta %</th><th>Flag</th></tr>",
		table_rows,
		"</table>",
		"<h2>QC plots</h2>",
		image_html,
		"</body></html>"
	)

	writeLines(html, report_path)
	report_path
}

main <- function() {
	# Main orchestrator for discovery, quantification, plotting, and reporting.
	args <- parse_cli_args(commandArgs(trailingOnly = TRUE))

	if (is.null(args$output_dir) || !nzchar(args$output_dir)) {
		args$output_dir <- file.path(dirname(args$input_dir), "output")
	}

	dir.create(args$output_dir, recursive = TRUE, showWarnings = FALSE)
	plot_dir <- file.path(args$output_dir, "plots")
	dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
	# Remove legacy RT-bias artifact if present from older workflow versions.
	unlink(file.path(args$output_dir, "qc_psm_rt_bias.csv"), force = TRUE)

	message("[1/6] Discovering input files")
	index <- discover_input_files(args$input_dir)

	message("[2/6] Loading tabular data")
	tbls <- load_tables(index)

	message("[3/6] Extracting run-level metrics")
	insights <- extract_statistics_insights(tbls)
	result_stats <- extract_result_statistics(tbls)
	summary_metrics <- build_summary_metrics(insights, result_stats)

	message("[4/6] Computing shared-feature and bias analyses")
	peptide_results <- analyze_shared_peptides(tbls)
	peptide_intensity <- build_peptide_intensity_table(tbls)
	protein_shared <- analyze_protein_groups(tbls)
	psm_quality <- build_psm_quality_table(tbls)
	msms_quality <- build_msms_quality_table(tbls)
	bias <- compute_bias_metrics(peptide_results$shared)
	duty_cycle <- compute_duty_cycle(tbls)

	message("[5/6] Writing tabular outputs")
	write_csv(summary_metrics, file.path(args$output_dir, "qc_metrics_summary.csv"))
	write_csv(peptide_results$shared, file.path(args$output_dir, "qc_shared_peptides.csv"))
	write_csv(peptide_results$overlap_summary, file.path(args$output_dir, "qc_overlap_summary.csv"))
	write_csv(protein_shared, file.path(args$output_dir, "qc_shared_protein_groups.csv"))
	write_csv(bias$peptide_bias, file.path(args$output_dir, "qc_peptide_bias_summary.csv"))
	write_csv(duty_cycle$summary, file.path(args$output_dir, "qc_duty_cycle.csv"))

	message("[6/6] Generating plots and report")
	save_plots(plot_dir, insights, peptide_results$shared, peptide_intensity, psm_quality, msms_quality, duty_cycle)

	if (isTRUE(args$render_report)) {
		report_path <- write_html_report(args$output_dir, summary_metrics, bias$peptide_bias, peptide_results$overlap_summary, duty_cycle$summary, plot_dir)
		message(sprintf("Report written: %s", report_path))
	}

	if (isTRUE(args$write_preprocess_script)) {
		preprocess_script <- write_preprocess_bash_script(args$output_dir, args$input_dir)
		message(sprintf("Preprocess helper script written: %s", preprocess_script))
	}

	message("QC workflow complete.")
}

if (identical(environment(), globalenv())) {
	main()
}
