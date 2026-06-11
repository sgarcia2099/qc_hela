#!/usr/bin/env bash
set -euo pipefail

input_dir="."
output_csv=""
name_glob="*MSMSSpectrumInfo*.txt"
recursive="false"

usage() {
  cat <<'EOF'
Usage:
  msms_duty_cycle_batch_fast.sh [options]

Options:
  --input-dir=PATH      Directory containing MSMSSpectrumInfo files (default: .)
  --output-csv=PATH     Output CSV path (default: <input-dir>/msms_duty_cycle_metrics.csv)
  --name-glob=GLOB      File glob for matching MSMSSpectrumInfo files (default: *MSMSSpectrumInfo*.txt)
  --recursive=true|false  Recurse into subdirectories (default: false)
  -h, --help            Show this help

Output columns:
  sample_id,file,file_path,status,error,total_spectra,acquisition_window_s,
  scan_rate_hz,total_ion_inject_ms,mean_ion_inject_ms,median_ion_inject_ms,
  mean_cycle_time_ms,median_cycle_time_ms,duty_cycle_pct
EOF
}

for arg in "$@"; do
  case "$arg" in
    --input-dir=*) input_dir="${arg#*=}" ;;
    --output-csv=*) output_csv="${arg#*=}" ;;
    --name-glob=*) name_glob="${arg#*=}" ;;
    --recursive=*) recursive="${arg#*=}" ;;
    -h|--help)
      usage
      exit 0
      ;;
    *)
      echo "Unknown argument: $arg" >&2
      usage >&2
      exit 1
      ;;
  esac
done

if [[ ! -d "$input_dir" ]]; then
  echo "Input directory does not exist: $input_dir" >&2
  exit 1
fi

if [[ -z "$output_csv" ]]; then
  output_csv="$input_dir/msms_duty_cycle_metrics.csv"
fi

find_args=("$input_dir")
if [[ "${recursive,,}" == "true" || "${recursive,,}" == "1" || "${recursive,,}" == "yes" ]]; then
  find_args+=( -type f -name "$name_glob" -print0 )
else
  find_args+=( -maxdepth 1 -type f -name "$name_glob" -print0 )
fi

mapfile -d '' files < <(find "${find_args[@]}" | sort -z)

if [[ "${#files[@]}" -eq 0 ]]; then
  echo "No files matched '$name_glob' in: $input_dir" >&2
  exit 1
fi

mkdir -p "$(dirname "$output_csv")"

echo "Found ${#files[@]} MSMSSpectrumInfo file(s)."

awk -v total_files="${#files[@]}" '
BEGIN {
  OFS = ","
  processed_files = 0
  progress_step = (total_files >= 200 ? 5 : 1)
  print "sample_id,file,file_path,status,error,total_spectra,acquisition_window_s,scan_rate_hz,total_ion_inject_ms,mean_ion_inject_ms,median_ion_inject_ms,mean_cycle_time_ms,median_cycle_time_ms,duty_cycle_pct"
}

function reset_state(    k) {
  missing_cols = ""
  sample_id = ""
  file_name = ""
  status = ""
  err = ""

  scan_idx = 0
  rt_idx = 0
  ion_idx = 0

  scan_count = 0
  ion_n = 0
  cycle_n = 0
  total_ion = 0
  min_rt = 0
  max_rt = 0

  delete seen_scan
  delete rt_by_scan
  delete ion_by_scan
  delete ion_vals
  delete cycle_vals
  delete sorted_scans
}

function q(s,   t) {
  t = s
  gsub(/"/, "\"\"", t)
  return "\"" t "\""
}

function clean(s,   t) {
  t = s
  gsub(/\r/, "", t)
  gsub(/^"/, "", t)
  gsub(/"$/, "", t)
  return t
}

function parse_line(line, arr, sep,   n, i) {
  n = split(line, arr, sep)
  for (i = 1; i <= n; i++) arr[i] = clean(arr[i])
  return n
}

function fmt(x) {
  if (x == "" || x != x) return ""
  return sprintf("%.3f", x)
}

function median(arr, n,   i, mid, tmp) {
  if (n < 1) return ""
  delete tmp
  for (i = 1; i <= n; i++) tmp[i] = arr[i]
  asort(tmp)
  mid = int((n + 1) / 2)
  if (n % 2 == 1) return tmp[mid]
  return (tmp[mid] + tmp[mid + 1]) / 2
}

function print_na_row() {
  print q(sample_id), q(file_name), q(FILENAME), q(status), q(err), "", "", "", "", "", "", "", "", ""
}

FNR == 1 {
  reset_state()

  file_name = FILENAME
  sub(/^.*\//, "", file_name)
  sample_id = file_name
  sub(/\.txt$/, "", sample_id)
  sub(/_MSMSSpectrumInfo$/, "", sample_id)

  sep = (index($0, "\t") > 0 ? "\t" : ",")
  header_n = parse_line($0, header, sep)

  for (i = 1; i <= header_n; i++) {
    if (header[i] == "First Scan") scan_idx = i
    else if (header[i] == "RT in min") rt_idx = i
    else if (header[i] == "Ion Inject Time in ms") ion_idx = i
  }

  if (scan_idx == 0) missing_cols = missing_cols (missing_cols ? ";" : "") "First Scan"
  if (rt_idx == 0) missing_cols = missing_cols (missing_cols ? ";" : "") "RT in min"
  if (ion_idx == 0) missing_cols = missing_cols (missing_cols ? ";" : "") "Ion Inject Time in ms"

  next
}

{
  if (missing_cols != "") next

  row_n = parse_line($0, row, sep)

  if (row_n < scan_idx || row_n < rt_idx || row_n < ion_idx) next

  scan_raw = row[scan_idx]
  rt_raw = row[rt_idx]
  ion_raw = row[ion_idx]

  if (scan_raw == "" || rt_raw == "") next

  scan = int(scan_raw + 0)
  rt = rt_raw + 0
  if (rt < 0) next

  ion = (ion_raw == "" ? 0 : ion_raw + 0)

  # Keep the first occurrence of each scan to avoid duplicate-induced distortion.
  if (!(scan in seen_scan)) {
    seen_scan[scan] = 1
    rt_by_scan[scan] = rt
    ion_by_scan[scan] = ion

    ion_vals[++ion_n] = ion
    total_ion += ion

    if (scan_count == 0 || rt < min_rt) min_rt = rt
    if (scan_count == 0 || rt > max_rt) max_rt = rt
    scan_count++
  }
}

ENDFILE {
  processed_files++

  if (missing_cols != "") {
    status = "missing_required_columns"
    err = missing_cols
    print_na_row()
  } else if (scan_count == 0) {
    status = "no_valid_rows"
    err = "No rows with valid First Scan and RT in min"
    print q(sample_id), q(file_name), q(FILENAME), q(status), q(err), 0, "", "", "", "", "", "", "", ""
  } else {
    acq_s = (max_rt - min_rt) * 60
    scan_rate_hz = (acq_s > 0 ? scan_count / acq_s : "")
    duty_cycle_pct = (acq_s > 0 ? (total_ion / (acq_s * 1000)) * 100 : "")

    mean_ion = (ion_n > 0 ? total_ion / ion_n : "")
    median_ion = median(ion_vals, ion_n)

    n_sorted = asorti(seen_scan, sorted_scans, "@ind_num_asc")
    cycle_sum = 0
    for (i = 2; i <= n_sorted; i++) {
      prev_scan = sorted_scans[i - 1]
      curr_scan = sorted_scans[i]
      cycle = (rt_by_scan[curr_scan] - rt_by_scan[prev_scan]) * 60 * 1000
      cycle_vals[++cycle_n] = cycle
      cycle_sum += cycle
    }

    mean_cycle = (cycle_n > 0 ? cycle_sum / cycle_n : "")
    median_cycle = median(cycle_vals, cycle_n)

    status = "ok"
    err = ""

    print q(sample_id), q(file_name), q(FILENAME), q(status), q(err),
          scan_count, fmt(acq_s), fmt(scan_rate_hz), fmt(total_ion), fmt(mean_ion), fmt(median_ion),
          fmt(mean_cycle), fmt(median_cycle), fmt(duty_cycle_pct)
  }

  if (processed_files % progress_step == 0 || processed_files == total_files) {
    printf("\rProgress: %d/%d files (%.1f%%)", processed_files, total_files, (processed_files * 100.0) / total_files) > "/dev/stderr"
  }
}

END {
  if (processed_files > 0) {
    printf("\n") > "/dev/stderr"
  }
}
' "${files[@]}" > "$output_csv"

ok_count=$(awk -F, 'NR>1 && $4=="""ok""" {n++} END{print n+0}' "$output_csv")
issue_count=$(( ${#files[@]} - ok_count ))

echo "Wrote ${#files[@]} row(s) to: $output_csv"
echo "Status summary: ok=$ok_count, issues=$issue_count"
