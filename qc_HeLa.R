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

parse_cli_args <- function(args) {
	# Defaults are intentionally data-centric so outputs sit with the dataset.
	defaults <- list(
		input_dir = "/home/jkg/OneDrive/Documents/All Files - personal/School_UT/Work/SG-01-008_BlackSoldierFlyGutMicrobiome/qc_hela/data",
		output_dir = NULL,
		render_report = TRUE,
		write_preprocess_script = TRUE
	)

	if (length(args) == 0) {
		return(defaults)
	}

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
		"^Consensus Features -",
		"^LCMS Features -",
		"^MS/MS Spectrum Info -",
		"^PSMs - Precursor Abundance$",
		"^PSMs - q-Value$",
		"^PSMs - PEP$",
		"^PSMs - SVM Score$",
		"^PSMs - Charge$",
		"^PSMs - RT \\[min\\]$",
		"^PSMs - Isolation Interference"
	)

	pattern <- paste(keep_patterns, collapse = "|")

	bind_rows(
		get_table(tbls, "preRun", "ResultStatistics") |> mutate(run = "preRun"),
		get_table(tbls, "postRun", "ResultStatistics") |> mutate(run = "postRun")
	) |>
		filter(str_detect(Name, pattern)) |>
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

	pre_rt_col <- names(pre)[str_detect(names(pre), "Top Apex RT")][1]
	post_rt_col <- names(post)[str_detect(names(post), "Top Apex RT")][1]

	pre_df <- pre |>
		transmute(
			feature_key = paste(`Annotated Sequence`, Modifications, sep = "||"),
			abundance_pre = suppressWarnings(as.numeric(.data[[pre_abundance_col]])),
			rt_pre = if (!is.na(pre_rt_col)) suppressWarnings(as.numeric(.data[[pre_rt_col]])) else NA_real_
		) |>
		group_by(feature_key) |>
		summarise(
			abundance_pre = sum(abundance_pre, na.rm = TRUE),
			rt_pre = median(rt_pre, na.rm = TRUE),
			.groups = "drop"
		)

	post_df <- post |>
		transmute(
			feature_key = paste(`Annotated Sequence`, Modifications, sep = "||"),
			abundance_post = suppressWarnings(as.numeric(.data[[post_abundance_col]])),
			rt_post = if (!is.na(post_rt_col)) suppressWarnings(as.numeric(.data[[post_rt_col]])) else NA_real_
		) |>
		group_by(feature_key) |>
		summarise(
			abundance_post = sum(abundance_post, na.rm = TRUE),
			rt_post = median(rt_post, na.rm = TRUE),
			.groups = "drop"
		)

	shared <- inner_join(pre_df, post_df, by = "feature_key") |>
		mutate(
			log2_fc_post_pre = safe_log2fc(abundance_post, abundance_pre),
			mean_log10_abundance = log10((abundance_pre + abundance_post) / 2 + 1),
			rt_mean = rowMeans(cbind(rt_pre, rt_post), na.rm = TRUE)
		)

	overlap_summary <- tibble(
		metric = c("peptides_pre_total", "peptides_post_total", "peptides_shared_total", "peptides_overlap_pct"),
		value = c(nrow(pre_df), nrow(post_df), nrow(shared), ifelse(nrow(pre_df) > 0, nrow(shared) / nrow(pre_df) * 100, NA_real_))
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
		rt_col <- fetch_col(names(df), "RT in min")
		q_col <- fetch_col(names(df), "^q-Value$")
		pep_col <- fetch_col(names(df), "^PEP$")
		interference_col <- fetch_col(names(df), "Isolation Interference")

		tibble(
			run = run_label,
			precursor_abundance = if (!is.na(abundance_col)) suppressWarnings(as.numeric(df[[abundance_col]])) else NA_real_,
			rt_min = if (!is.na(rt_col)) suppressWarnings(as.numeric(df[[rt_col]])) else NA_real_,
			q_value = if (!is.na(q_col)) suppressWarnings(as.numeric(df[[q_col]])) else NA_real_,
			pep = if (!is.na(pep_col)) suppressWarnings(as.numeric(df[[pep_col]])) else NA_real_,
			isolation_interference = if (!is.na(interference_col)) suppressWarnings(as.numeric(df[[interference_col]])) else NA_real_
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
		rt_col <- fetch_col(names(df), "RT in min")
		interference_col <- fetch_col(names(df), "Isolation Interference")

		tibble(
			run = run_label,
			rt_min = if (!is.na(rt_col)) suppressWarnings(as.numeric(df[[rt_col]])) else NA_real_,
			isolation_interference = if (!is.na(interference_col)) suppressWarnings(as.numeric(df[[interference_col]])) else NA_real_
		)
	}

	bind_rows(map_msms(pre, "preRun"), map_msms(post, "postRun"))
}

compute_bias_metrics <- function(shared_peptides, psm_quality) {
	peptide_bias <- shared_peptides |>
		summarise(
			shared_features = n(),
			median_log2_fc = median(log2_fc_post_pre, na.rm = TRUE),
			frac_negative_fc = mean(log2_fc_post_pre < 0, na.rm = TRUE),
			frac_abs_fc_gt1 = mean(abs(log2_fc_post_pre) > 1, na.rm = TRUE),
			wilcox_p = tryCatch(wilcox.test(log2_fc_post_pre, mu = 0)$p.value, error = function(e) NA_real_)
		)

	psm_bias <- psm_quality |>
		filter(!is.na(precursor_abundance), precursor_abundance > 0) |>
		mutate(rt_bin = ntile(rt_min, 10)) |>
		group_by(run, rt_bin) |>
		summarise(median_log10_abundance = median(log10(precursor_abundance + 1), na.rm = TRUE), .groups = "drop") |>
		pivot_wider(names_from = run, values_from = median_log10_abundance) |>
		mutate(post_minus_pre = postRun - preRun)

	list(peptide_bias = peptide_bias, psm_bias_by_rt = psm_bias)
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

save_plots <- function(out_plot_dir, insights, shared_peptides, psm_quality, msms_quality, psm_bias) {
	dir.create(out_plot_dir, recursive = TRUE, showWarnings = FALSE)

	counts_plot <- insights |>
		filter(str_detect(metric, "high_confident")) |>
		mutate(metric_label = str_replace_all(metric, "_", " ")) |>
		ggplot(aes(x = metric_label, y = value, fill = run)) +
		geom_col(position = position_dodge(width = 0.7), width = 0.65) +
		coord_flip() +
		theme_minimal(base_size = 12) +
		labs(title = "High-Confidence IDs: Pre vs Post", x = NULL, y = "Count")
	ggsave(file.path(out_plot_dir, "01_high_confidence_counts.png"), counts_plot, width = 9, height = 5, dpi = 140)

	scatter <- ggplot(shared_peptides, aes(x = log10(abundance_pre + 1), y = log10(abundance_post + 1))) +
		geom_point(alpha = 0.35, size = 1.1, color = "#1f77b4") +
		geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "#d62728") +
		theme_minimal(base_size = 12) +
		labs(
			title = "Shared Peptides: Pre vs Post Abundance",
			x = "log10 Pre abundance",
			y = "log10 Post abundance"
		)
	ggsave(file.path(out_plot_dir, "02_shared_peptide_scatter.png"), scatter, width = 7, height = 6, dpi = 140)

	ma <- ggplot(shared_peptides, aes(x = mean_log10_abundance, y = log2_fc_post_pre)) +
		geom_point(alpha = 0.3, size = 1, color = "#2ca02c") +
		geom_hline(yintercept = 0, linetype = "dashed", color = "#d62728") +
		theme_minimal(base_size = 12) +
		labs(title = "MA Plot (Shared Peptides)", x = "Mean log10 abundance", y = "log2(Post/Pre)")
	ggsave(file.path(out_plot_dir, "03_ma_plot.png"), ma, width = 7, height = 6, dpi = 140)

	hist_fc <- ggplot(shared_peptides, aes(x = log2_fc_post_pre)) +
		geom_histogram(bins = 70, fill = "#4c78a8", color = "white", alpha = 0.9) +
		geom_vline(xintercept = 0, linetype = "dashed", color = "#d62728") +
		theme_minimal(base_size = 12) +
		labs(title = "Distribution of Shared-Peptide log2 Fold-Change", x = "log2(Post/Pre)", y = "Features")
	ggsave(file.path(out_plot_dir, "04_log2fc_histogram.png"), hist_fc, width = 8, height = 5.5, dpi = 140)

	rt_trend <- shared_peptides |>
		filter(!is.na(rt_mean)) |>
		mutate(rt_bin = ntile(rt_mean, 10)) |>
		group_by(rt_bin) |>
		summarise(median_log2_fc = median(log2_fc_post_pre, na.rm = TRUE), .groups = "drop") |>
		ggplot(aes(x = rt_bin, y = median_log2_fc)) +
		geom_line(color = "#ff7f0e", linewidth = 1) +
		geom_point(color = "#ff7f0e", size = 2) +
		geom_hline(yintercept = 0, linetype = "dashed", color = "#d62728") +
		theme_minimal(base_size = 12) +
		labs(title = "Median log2(Post/Pre) by RT Decile", x = "RT decile", y = "Median log2 FC")
	ggsave(file.path(out_plot_dir, "05_rt_binned_log2fc_trend.png"), rt_trend, width = 8, height = 5, dpi = 140)

	psm_abundance <- psm_quality |>
		filter(!is.na(precursor_abundance), precursor_abundance > 0) |>
		mutate(log10_precursor_abundance = log10(precursor_abundance + 1)) |>
		ggplot(aes(x = run, y = log10_precursor_abundance, fill = run)) +
		geom_boxplot(outlier.alpha = 0.15) +
		theme_minimal(base_size = 12) +
		labs(title = "PSM Precursor Abundance by Run", x = NULL, y = "log10 abundance")
	ggsave(file.path(out_plot_dir, "06_psm_precursor_abundance_boxplot.png"), psm_abundance, width = 7, height = 5, dpi = 140)

	interference <- msms_quality |>
		filter(!is.na(isolation_interference)) |>
		ggplot(aes(x = run, y = isolation_interference, fill = run)) +
		geom_violin(trim = TRUE, alpha = 0.85) +
		theme_minimal(base_size = 12) +
		labs(title = "Isolation Interference by Run", x = NULL, y = "Isolation Interference (%)")
	ggsave(file.path(out_plot_dir, "07_isolation_interference_violin.png"), interference, width = 7, height = 5, dpi = 140)

	psm_rt_trend <- psm_bias |>
		ggplot(aes(x = rt_bin, y = post_minus_pre)) +
		geom_col(fill = "#9467bd") +
		geom_hline(yintercept = 0, linetype = "dashed", color = "#d62728") +
		theme_minimal(base_size = 12) +
		labs(title = "Post-Pre Median log10 Precursor Abundance by RT Decile", x = "RT decile", y = "Post - Pre")
	ggsave(file.path(out_plot_dir, "08_psm_rt_bias_bar.png"), psm_rt_trend, width = 8, height = 5, dpi = 140)
}

write_html_report <- function(output_dir, summary_metrics, peptide_bias, plot_dir) {
	report_path <- file.path(output_dir, "qc_report.html")
	key <- peptide_bias |>
		mutate(across(everything(), ~ format(., digits = 4)))

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
		sprintf("<p>Median log2(Post/Pre): %s</p>", key$median_log2_fc),
		sprintf("<p>Fraction with negative fold-change: %s</p>", key$frac_negative_fc),
		sprintf("<p>Wilcoxon p-value (log2 FC vs 0): %s</p>", key$wilcox_p),
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

	dir.create(args$output_dir, recursive = TRUE, showWarnings = FALSE)
	plot_dir <- file.path(args$output_dir, "plots")
	dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

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
	protein_shared <- analyze_protein_groups(tbls)
	psm_quality <- build_psm_quality_table(tbls)
	msms_quality <- build_msms_quality_table(tbls)
	bias <- compute_bias_metrics(peptide_results$shared, psm_quality)

	message("[5/6] Writing tabular outputs")
	write_csv(summary_metrics, file.path(args$output_dir, "qc_metrics_summary.csv"))
	write_csv(peptide_results$shared, file.path(args$output_dir, "qc_shared_peptides.csv"))
	write_csv(peptide_results$overlap_summary, file.path(args$output_dir, "qc_overlap_summary.csv"))
	write_csv(protein_shared, file.path(args$output_dir, "qc_shared_protein_groups.csv"))
	write_csv(bias$psm_bias_by_rt, file.path(args$output_dir, "qc_psm_rt_bias.csv"))
	write_csv(bias$peptide_bias, file.path(args$output_dir, "qc_peptide_bias_summary.csv"))

	message("[6/6] Generating plots and report")
	save_plots(plot_dir, insights, peptide_results$shared, psm_quality, msms_quality, bias$psm_bias_by_rt)

	if (isTRUE(args$render_report)) {
		report_path <- write_html_report(args$output_dir, summary_metrics, bias$peptide_bias, plot_dir)
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
