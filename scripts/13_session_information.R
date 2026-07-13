############################################################
##
## S3011 Perennial Sorghum GWAS Pipeline
##
## Script 13: Reproducibility and Session Information
##
## Purpose:
##   - Record R, operating-system, and package information
##   - Record versions of packages used across the pipeline
##   - Record Git commit information when available
##   - Record MD5 checksums for scripts and key input files
##   - Produce machine-readable and human-readable
##     reproducibility summaries
##
## Required preceding script:
##   00_setup.R
##
## Outputs:
##   results/reproducibility/
##     sessionInfo.txt
##     reproducibility_metadata.csv
##     package_versions.csv
##     script_checksums.csv
##     input_file_checksums.csv
##     git_information.txt
##
############################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tibble)
})

############################################################
## 1. Verify required pipeline objects
############################################################

required_objects <- c(
  "root_s3011"
)

missing_objects <- required_objects[
  !vapply(
    required_objects,
    exists,
    logical(1),
    inherits = TRUE
  )
]

if (length(missing_objects) > 0L) {
  stop(
    "Required pipeline objects are missing: ",
    paste(missing_objects, collapse = ", "),
    ". Run 00_setup.R before running this script."
  )
}

############################################################
## 2. Define output and project directories
############################################################

reproducibility_dir <- file.path(
  root_s3011,
  "results",
  "reproducibility"
)

scripts_dir <- file.path(
  root_s3011,
  "scripts"
)

data_dir <- file.path(
  root_s3011,
  "data"
)

dir.create(
  reproducibility_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

############################################################
## 3. Save complete session information
############################################################

session_information <- capture.output(
  sessionInfo()
)

writeLines(
  session_information,
  file.path(
    reproducibility_dir,
    "sessionInfo.txt"
  )
)

############################################################
## 4. Record general reproducibility metadata
############################################################

system_information <- Sys.info()

metadata <- tibble(
  Field = c(
    "Run timestamp",
    "Project root",
    "R version",
    "R platform",
    "R architecture",
    "Operating system",
    "Node name",
    "User",
    "Locale",
    "Working directory"
  ),

  Value = c(
    format(
      Sys.time(),
      "%Y-%m-%d %H:%M:%S %Z"
    ),

    normalizePath(
      root_s3011,
      winslash = "/",
      mustWork = FALSE
    ),

    R.version.string,

    R.version$platform,

    R.version$arch,

    paste(
      na.omit(
        c(
          system_information[
            "sysname"
          ],
          system_information[
            "release"
          ],
          system_information[
            "version"
          ],
          system_information[
            "machine"
          ]
        )
      ),
      collapse = " | "
    ),

    unname(
      system_information[
        "nodename"
      ]
    ),

    unname(
      system_information[
        "user"
      ]
    ),

    Sys.getlocale(),

    normalizePath(
      getwd(),
      winslash = "/",
      mustWork = FALSE
    )
  )
)

write_csv(
  metadata,
  file.path(
    reproducibility_dir,
    "reproducibility_metadata.csv"
  )
)

############################################################
## 5. Record package versions
############################################################

pipeline_packages <- c(
  "data.table",
  "dplyr",
  "tidyr",
  "readr",
  "stringr",
  "purrr",
  "tibble",
  "ggplot2",
  "lme4",
  "qqman"
)

package_versions <- bind_rows(
  lapply(
    pipeline_packages,
    function(package_name) {

      installed <- requireNamespace(
        package_name,
        quietly = TRUE
      )

      tibble(
        Package =
          package_name,

        Installed =
          installed,

        Version =
          if (
            installed
          ) {
            as.character(
              packageVersion(
                package_name
              )
            )
          } else {
            NA_character_
          }
      )
    }
  )
)

write_csv(
  package_versions,
  file.path(
    reproducibility_dir,
    "package_versions.csv"
  )
)

############################################################
## 6. Helper to calculate file checksums
############################################################

build_checksum_table <- function(
  files,
  project_root
) {

  if (length(files) == 0L) {
    return(
      tibble(
        Relative_Path =
          character(),
        Size_Bytes =
          numeric(),
        Modified_Time =
          character(),
        MD5 =
          character()
      )
    )
  }

  files <- files[
    file.exists(
      files
    )
  ]

  if (length(files) == 0L) {
    return(
      tibble(
        Relative_Path =
          character(),
        Size_Bytes =
          numeric(),
        Modified_Time =
          character(),
        MD5 =
          character()
      )
    )
  }

  file_information <- file.info(
    files
  )

  relative_paths <- sub(
    paste0(
      "^",
      gsub(
        "([.\\+*?\\[\\](){}^$|])",
        "\\\\\\1",
        normalizePath(
          project_root,
          winslash = "/",
          mustWork = FALSE
        )
      ),
      "/?"
    ),
    "",
    normalizePath(
      files,
      winslash = "/",
      mustWork = FALSE
    )
  )

  tibble(
    Relative_Path =
      relative_paths,

    Size_Bytes =
      file_information$size,

    Modified_Time =
      format(
        file_information$mtime,
        "%Y-%m-%d %H:%M:%S %Z"
      ),

    MD5 =
      unname(
        tools::md5sum(
          files
        )
      )
  ) %>%
    arrange(
      Relative_Path
    )
}

############################################################
## 7. Record script checksums
############################################################

script_files <- if (
  dir.exists(
    scripts_dir
  )
) {
  list.files(
    scripts_dir,
    pattern = "\\.[Rr]$",
    recursive = TRUE,
    full.names = TRUE
  )
} else {
  character()
}

script_checksums <- build_checksum_table(
  files =
    script_files,
  project_root =
    root_s3011
)

write_csv(
  script_checksums,
  file.path(
    reproducibility_dir,
    "script_checksums.csv"
  )
)

############################################################
## 8. Record input-file checksums
############################################################

## Large generated result files are intentionally excluded.
## This section records files stored under the project data
## directory. Add other immutable inputs here when necessary.

input_files <- if (
  dir.exists(
    data_dir
  )
) {
  list.files(
    data_dir,
    recursive = TRUE,
    full.names = TRUE
  )
} else {
  character()
}

input_file_checksums <- build_checksum_table(
  files =
    input_files,
  project_root =
    root_s3011
)

write_csv(
  input_file_checksums,
  file.path(
    reproducibility_dir,
    "input_file_checksums.csv"
  )
)

############################################################
## 9. Record Git information when available
############################################################

run_git_command <- function(arguments) {

  result <- tryCatch(
    system2(
      command = "git",
      args = arguments,
      stdout = TRUE,
      stderr = TRUE
    ),
    error = function(e) {
      paste0(
        "Git command failed: ",
        conditionMessage(
          e
        )
      )
    }
  )

  paste(
    result,
    collapse = "\n"
  )
}

git_available <- nzchar(
  Sys.which(
    "git"
  )
)

git_output <- c(
  paste0(
    "Git available: ",
    git_available
  )
)

if (git_available) {

  repository_check <- run_git_command(
    c(
      "-C",
      shQuote(
        root_s3011
      ),
      "rev-parse",
      "--is-inside-work-tree"
    )
  )

  git_output <- c(
    git_output,
    paste0(
      "Inside Git repository: ",
      repository_check
    )
  )

  if (
    identical(
      trimws(
        repository_check
      ),
      "true"
    )
  ) {

    git_output <- c(
      git_output,
      "",
      "Commit:",
      run_git_command(
        c(
          "-C",
          shQuote(
            root_s3011
          ),
          "rev-parse",
          "HEAD"
        )
      ),
      "",
      "Branch:",
      run_git_command(
        c(
          "-C",
          shQuote(
            root_s3011
          ),
          "rev-parse",
          "--abbrev-ref",
          "HEAD"
        )
      ),
      "",
      "Repository status:",
      run_git_command(
        c(
          "-C",
          shQuote(
            root_s3011
          ),
          "status",
          "--short"
        )
      ),
      "",
      "Most recent commit:",
      run_git_command(
        c(
          "-C",
          shQuote(
            root_s3011
          ),
          "log",
          "-1",
          "--pretty=format:%H%n%an%n%ad%n%s",
          "--date=iso"
        )
      )
    )
  }

} else {

  git_output <- c(
    git_output,
    "Git executable was not found on this system."
  )
}

writeLines(
  git_output,
  file.path(
    reproducibility_dir,
    "git_information.txt"
  )
)

############################################################
## 10. Completion messages
############################################################

message(
  "13_session_information.R completed successfully."
)

message(
  "Reproducibility outputs saved in: ",
  reproducibility_dir
)

message(
  "Scripts recorded: ",
  nrow(
    script_checksums
  )
)

message(
  "Input files recorded: ",
  nrow(
    input_file_checksums
  )
)

message(
  "Packages checked: ",
  nrow(
    package_versions
  )
)
