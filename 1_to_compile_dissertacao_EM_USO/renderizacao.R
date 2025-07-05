bookdown::render_book("compilacao_integral.Rmd", "bookdown::pdf_document2")
# bookdown::render_book("compilacao_integral.Rmd", "bookdown::word_document2")
library(pdftools)
pdf_combine(c("docs/_main.pdf",
              "~/Documentos/mestrado_Ecologia/disciplinas/BIE5716/relatorio_final.pdf"), 
            output = "docs/merged.pdf")