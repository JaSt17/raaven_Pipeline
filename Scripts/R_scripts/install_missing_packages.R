.libPaths(c('~/MyRextensions', .libPaths()))

#set a CRan mirror
options(repos = c(CRAN = "https://ftp.acc.umu.se/mirror/CRAN/"))

# Function to check if a package is installed, and install it if not
install_if_missing <- function(pkg) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
    library(pkg, character.only = TRUE)
  }
}

# List of packages used in this script
packages <- c("ShortRead", "parallel", "doParallel", "data.table", "devtools", "Hmisc", "kableExtra")

# Install missing packages
for (pkg in packages) {
  install_if_missing(pkg)
}