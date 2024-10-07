#' Install Python and necessary packages using Miniconda
#'
#' This function ensures that the required Python environment and packages are installed.
#'
#' @param conda Path to the Conda executable or "auto" to automatically detect it.
#' @param python.version Version of Python to install (default is "3.8").
#' @param install.conda Logical indicating whether to install Miniconda if not detected (default is FALSE).
#' @param miniconda.path Path to the Miniconda installation directory (default is NULL, which uses the default path).
#'
#' @return None
#' @export
#'
#' @examples
#' installPython()
installPython <- function(
    conda = "auto",
    python.version = "3.8",
    install.conda = FALSE,
    miniconda.path = NULL
) {
  suppressMessages(require(reticulate))

  # Check if Conda is available
  if (!.isConda()) {
    if (!install.conda) {
      stop("No miniconda detected, but 'install.conda' is FALSE. Please, set ",
           "'install.conda = TRUE' to install miniconda." )
    }
    message("=== No miniconda detected, installing through the reticulate R package")
    if (is.null(miniconda.path)) {
      miniconda.path <- reticulate::miniconda_path()
    }
    status1 <- tryCatch(
      reticulate::install_miniconda(path = miniconda.path),
      error = function(e) {
        return(TRUE)
      }
    )
    if (isTRUE(status1)) {
      stop(
        "Error during the installation. Please see the website of the ",
        "package and/or the vignettes for more details",
        call. = FALSE
      )
    }
  }

  # Get the path to the Conda binary
  dirConda <- reticulate::conda_binary("auto")
  message("\n=== Creating SSLAR-env environment")


  # Custom Python version
  # if (python.version != "3.8") {
  #   warning(
  #     "Please, be sure the selected Python versions are ",
  #     "compatible. Otherwise, miniconda will raise an error",
  #     call. = FALSE, immediate. = TRUE
  #   )
  # }

  # Create the Conda environment
  status2 <- tryCatch(
    reticulate::conda_create(
      envname = "SSLAR-env",
      packages = paste0("python==", python.version)
    ),
    error = function(e) {
      return(TRUE)
    }
  )

  if (isTRUE(status2)) {
    stop(
      "Error during the installation. Please see the website of the ",
      "package and/or the vignettes for more details",
      call. = FALSE
    )
  }
  # Install TensorFlow in the environment

  message("\n=== Installing ssNMF in SSLAR-env environment")
  status3 <- tryCatch(
    reticulate::py_install("ssNMF",
                           conda = dirConda,
                           envname = "SSLAR-env",
                           pip = TRUE),
    error = function(e) {
      return(TRUE)
    }
  )
  if (isTRUE(status3)) {
    stop(
      "Error during the installation. Please see the website of the ssNMF",
      "package and/or the vignettes for more details",
      call. = FALSE
    )
  }

  message("\n=== Installing torch in SSLAR-env environment")
  status4 <- tryCatch(
    reticulate::py_install("torch",
                           conda = dirConda,
                           envname = "SSLAR-env",
                           pip = TRUE),
    error = function(e) {
      return(TRUE)
    }
  )
  if (isTRUE(status4)) {
    stop(
      "Error during the installation. Please see the website of the torch",
      "package and/or the vignettes for more details",
      call. = FALSE
    )
  }

  message("Installation complete!")
  message("Restart R and load the necessary packages.")
}



# Check if Conda is available
.isConda <- function() {
  conda <- tryCatch(
    reticulate::conda_binary("auto"), error = function(e) NULL
  )
  !is.null(conda)
}

# Check if Python is available
.isPython <- function() {
  tryCatch(
    expr = reticulate::py_available(initialize = TRUE),
    error = function(e) FALSE
  )
}

.onLoad <- function(libname, pkgname) {
  if (.isConda()) {
    tryCatch(
      expr = reticulate::use_condaenv("SSLAR-env", required = TRUE),
      error = function(e) NULL
    )
  }
}
