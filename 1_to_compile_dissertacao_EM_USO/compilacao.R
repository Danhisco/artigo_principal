rmarkdown::render(
  "1_to_compile_dissertacao_EM_USO/CongSAD.Rmd",
  output_file = "CongSAD.docx",
  envir = globalenv(),
  knit_root_dir = getwd(),
  clean = TRUE
)
