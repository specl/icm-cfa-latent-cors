out_dir = "_out"
cache_dir = "_cache"

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)

rmarkdown::render(
  input = "p.Rmd",
  output_dir = out_dir,
  intermediates_dir = cache_dir,
  clean = TRUE
)
