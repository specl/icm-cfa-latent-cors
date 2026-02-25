pkgs_file = file.path("pkgs.txt")
pkgs = character(0)

if (file.exists(pkgs_file)) {
  lines = trimws(readLines(pkgs_file, warn = FALSE))
  pkgs = unique(lines[nzchar(lines) & !grepl("^#", lines)])
}

.silent = function(expr) {
  suppressWarnings(suppressMessages({
    capture.output(capture.output(force(expr), type = "message"))
  }))
  invisible(NULL)
}

if (length(pkgs) > 0) {
  .silent(sapply(pkgs, function(p) {
    ok = require(p, character.only = TRUE, quietly = TRUE)
    if (!ok) {
      install.packages(p)
      ok = require(p, character.only = TRUE, quietly = TRUE)
    }
    ok
  }))
}



