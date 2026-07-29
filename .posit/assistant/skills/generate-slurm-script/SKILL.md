---
name: generate-slurm-script
description: Generate a paired SLURM .sh submission script for an existing .R script. Use when the user wants to create a shell script for submitting an R script as a batch job via sbatch, or mentions creating a .sh file for an .R script.
---

# Generate SLURM Shell Script

When asked to generate a `.sh` file for a given `.R` script, follow these steps:

1. **Read the R script** to determine:
   - Whether it accepts command-line arguments (check for `commandArgs()`)
   - If arguments come from a list of items (samples, cell types) → **array job**
   - If arguments represent parameter combinations → **array job with factorial design**
   - If no arguments are taken → **single job**

2. **Determine resource requirements** based on the R script's content:
   - Heavy object loading (SCE, Seurat, ArchR) → 80-100GB memory, `katun` partition
   - Lightweight analysis (plotting, stats) → 20-50GB memory, `shared` partition
   - ChromVAR, large matrix ops → 100GB, 10 cores, long time limit
   - Default: 50GB, 1 core, 1 day, `shared` partition

3. **Generate the shell script** using the appropriate template below.

## Naming Convention

The `.sh` file must have the same base name as the `.R` file and live in the same directory. For example, `code/17_species_diverg/01_hashikawa_markers.R` → `code/17_species_diverg/01_hashikawa_markers.sh`.

## Template A: Single Job (no array)

Use when the R script takes no command-line arguments.

```bash
#!/bin/bash
#SBATCH -p {partition}
#SBATCH --mem={memory}
#SBATCH --job-name={script_basename_no_ext}
#SBATCH -c {cpus}
#SBATCH -t {time_limit}
#SBATCH -o ../../processed-data/{subdir}/logs/{script_basename_no_ext}.txt
#SBATCH -e ../../processed-data/{subdir}/logs/{script_basename_no_ext}.txt

set -e

echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOB_ID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Node name: ${HOSTNAME}"

module load conda_R/4.5

## List current modules for reproducibility
module list

Rscript {script_basename}.R

echo "**** Job ends ****"
date
```

## Template B: Array Job with Item List (bash array of names)

Use when the R script iterates over a set of discrete items (cell types, samples, conditions) passed as a single `--flag value` argument.

```bash
#!/bin/bash
#SBATCH -p {partition}
#SBATCH --mem={memory}
#SBATCH --job-name={script_basename_no_ext}
#SBATCH -c {cpus}
#SBATCH -t {time_limit}
#SBATCH -o /dev/null
#SBATCH -e /dev/null
#SBATCH --array=0-{N-1}%{max_concurrent}

set -eo pipefail

# Bash array of items
items=(
  "item1" "item2" "item3"
)

# Allow local testing; SLURM sets this in the array
i=${SLURM_ARRAY_TASK_ID:-0}
m=${#items[@]}

# guard
if (( i < 0 || i >= m )); then
  echo "Invalid task index: $i (must be 0..$((m-1)))" >&2
  exit 1
fi

res="${items[$i]}"

mkdir -p logs
log_path="logs/{script_basename_no_ext}_${res}.log"

{
  echo "**** Job starts ****"
  date

  echo "**** JHPCE info ****"
  echo "User: ${USER}"
  echo "Job id: ${SLURM_JOB_ID}"
  echo "Job name: ${SLURM_JOB_NAME}"
  echo "Node name: ${HOSTNAME}"
  echo "Task id: ${SLURM_ARRAY_TASK_ID}"
  echo "Selected item: ${res}"

  ## Load the R module
  module load conda_R/4.5

  ## List current modules for reproducibility
  module list

  Rscript {script_basename}.R --{flag_name} "${res}"

  echo "**** Job ends ****"
  date

} > "$log_path" 2>&1
```

## Template C: Array Job with Factorial Parameter Combinations

Use when the R script accepts multiple parameters and you want to sweep over all combinations.

```bash
#!/bin/bash
#SBATCH -p {partition}
#SBATCH --mem={memory}
#SBATCH --job-name={script_basename_no_ext}
#SBATCH -c {cpus}
#SBATCH -t {time_limit}
#SBATCH -o /dev/null
#SBATCH -e /dev/null
#SBATCH --array=0-{total-1}%{max_concurrent}

set -eo pipefail

# Parameter arrays
param1_values=({values})   # m={count}
param2_values=({values})   # n={count}
# total = m * n = {total}

# Task mapping:
# 0: param1[0], param2[0]
# 1: param1[0], param2[1]
# ...

i=${SLURM_ARRAY_TASK_ID}
m=${#param1_values[@]}
n=${#param2_values[@]}
total=$(( m * n ))

if (( i < 0 || i >= total )); then
  echo "Invalid task index: $i (total=$total)"; exit 1
fi

param1_idx=$(( i / n ))
param2_idx=$(( i % n ))

p1="${param1_values[$param1_idx]}"
p2="${param2_values[$param2_idx]}"

mkdir -p logs
log_path="logs/{script_basename_no_ext}_${p1}_${p2}_task_${i}.log"

{
  echo "**** Job starts ****"
  date

  echo "**** JHPCE info ****"
  echo "User: ${USER}"
  echo "Job id: ${SLURM_JOB_ID}"
  echo "Job name: ${SLURM_JOB_NAME}"
  echo "Node name: ${HOSTNAME}"
  echo "Task id: ${SLURM_ARRAY_TASK_ID}"

  ## Load the R module
  module load conda_R/4.5

  ## List current modules for reproducibility
  module list

  Rscript {script_basename}.R --{flag1} "${p1}" --{flag2} "${p2}"

  echo "**** Job ends ****"
  date

} > "$log_path" 2>&1
```

## Template D: Array Job with External File (sample names from .txt file)

Use when sample IDs are read from a text file (one per line).

```bash
#!/bin/bash
#SBATCH -p {partition}
#SBATCH --mem={memory}
#SBATCH --job-name={script_basename_no_ext}
#SBATCH -c {cpus}
#SBATCH -t {time_limit}
#SBATCH -o /dev/null
#SBATCH -e /dev/null
#SBATCH --array=1-{N}

set -eo pipefail

id=$(sed -n ${SLURM_ARRAY_TASK_ID}p {targets_file}.txt)

mkdir -p logs
log_path="logs/{script_basename_no_ext}_${id}.log"

{
  echo "**** Job starts ****"
  date

  echo "**** JHPCE info ****"
  echo "User: ${USER}"
  echo "Job id: ${SLURM_JOB_ID}"
  echo "Job name: ${SLURM_JOB_NAME}"
  echo "Node name: ${HOSTNAME}"
  echo "Task id: ${SLURM_ARRAY_TASK_ID}"
  echo "Sample: ${id}"

  ## Load the R module
  module load conda_R/4.5

  ## List current modules for reproducibility
  module list

  Rscript {script_basename}.R $id

  echo "**** Job ends ****"
  date

} > "$log_path" 2>&1
```

## Key Decisions to Make

When generating a script, determine:

1. **Single vs. array**: Does the R script use `commandArgs()`? If yes → array. If no → single.
2. **Array type**: Are arguments from a list of items, a factorial combination, or an external file?
3. **Partition**: `katun` for high-memory jobs (≥80GB), `shared` for lighter work.
4. **Memory**: Based on data objects loaded. Default 50GB if uncertain; ask the user.
5. **Cores**: Match to parallelism in R script (e.g., `BiocParallel`, `mclapply` → multiple cores).
6. **Time**: Default 1 day; use 2-10 days for very large computations.
7. **Log path**: Use the `../../processed-data/{subdir}/logs/` convention for single jobs, or local `logs/` directory with brace-redirected output for array jobs.
8. **R module**: Default to `conda_R/4.5`. Check if script uses older packages that need `conda_R/4.3.x`.

## Important Notes

- Always use `set -e` (single jobs) or `set -eo pipefail` (array jobs with brace redirection).
- For array jobs, use `/dev/null` for SBATCH `-o` and `-e` and redirect output manually to task-specific log files.
- Include the `module list` call for reproducibility.
- The `%N` suffix on `--array` throttles concurrent tasks (e.g., `--array=0-17%10` runs max 10 at once). Default to the total number of tasks or 20, whichever is smaller.
- Ask the user if resource requirements are unclear from the R script.
