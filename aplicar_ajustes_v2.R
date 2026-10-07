# aplicar_ajustes_v2.R ------------------------------------------------------
# Gera compilacao_integral_v2.Rmd a partir de compilacao_integral.Rmd com:
#   (1) vlabel corrigido e usado no texto (Results) e na legenda da Fig. 2
#   (2) referencias cruzadas manuais -> \@ref() do bookdown
#   (3) figuras no tamanho final (largura util do docx) + chunk de exportacao
#   (4) subscritos: U_pristine -> U~pristine~, Ui/Uj -> U~i~/U~j~ (texto e legendas)
# O original NAO e alterado.
#
# Uso:
#   setwd("~/Documentos/artigo_principal")   # mesma raiz usada pelos readRDS()
#   source("1_to_compile_dissertacao_EM_USO/aplicar_ajustes_v2.R")
# (ajuste f_in / f_out se o Rmd estiver em outro lugar)

f_in  <- "1_to_compile_dissertacao_EM_USO/compilacao_integral.Rmd"
f_out <- "1_to_compile_dissertacao_EM_USO/compilacao_integral_v2.Rmd"

# Referencia cruzada da tabela: o caption do flextable NAO e alvo nativo de
# \@ref(tab:...). Deixe FALSE (mantem set_caption(), que ja numera "Table 1").
# Mude para TRUE apenas para testar tab.cap/tab.id (ver README no fim).
TAB_BOOKDOWN <- FALSE

x <- readLines(f_in, encoding = "UTF-8", warn = FALSE)

# ---- helpers ---------------------------------------------------------------
chunk_range <- function(x, label) {
  s <- grep(paste0("^```\\{r ", label, "[,}]"), x)
  stopifnot("chunk nao encontrado ou duplicado" = length(s) == 1)
  e <- s + which(trimws(x[(s + 1):length(x)]) == "```")[1]
  c(s, e)
}
replace_chunk <- function(x, label, new) {
  r <- chunk_range(x, label)
  c(x[seq_len(r[1] - 1)], new, x[-seq_len(r[2])])
}
insert_after_chunk <- function(x, label, new) {
  r <- chunk_range(x, label)
  c(x[seq_len(r[2])], "", new, x[-seq_len(r[2])])
}
lines_of <- function(s) strsplit(s, "\n", fixed = TRUE)[[1]]

# Indices das linhas de texto corrido (fora de chunks R e fora do cabecalho YAML).
# Um chunk comeca em "```{" e termina em "```": a soma acumulada de aberturas
# menos fechamentos e > 0 enquanto estamos DENTRO de um chunk.
idx_texto <- function(x) {
  starts   <- grepl("^```\\{", x)
  ends     <- trimws(x) == "```"
  in_chunk <- (cumsum(starts) - cumsum(ends)) > 0
  yaml_end <- which(trimws(x) == "---")[2]
  which(!in_chunk & seq_along(x) > yaml_end)
}

# ---- 1. chunk novo logo apos o setup ----------------------------------------
new_setup <- lines_of(r"--(```{r tamanho-figuras-e-valores,include=FALSE}
# Largura util do docx de referencia (polegadas). 6.7 in = 17 cm.
# Ajuste para (largura da pagina - margens) de resources/<ref>.docx.
W_FIG    <- 6.7
H_FIG1   <- 2.9    # Fig. 1 (3 paineis lado a lado)
H_FIG2   <- 6.2    # Fig. 2 (grade 3 x 5)
DPI_FIG  <- 300

# Valores citados no texto (inline) e na legenda da Fig. 2.
# Ordem criada em f_obs_predito_bysite(): min, max (amplitude); q05, q95.
dfrange <- readRDS("1_to_compile_dissertacao_EM_USO/00_Resultados/figuras/dfrange.rds")
stopifnot(nrow(dfrange) == 4)
vals   <- unname(dfrange$value)
vlabel <- setNames(sprintf("%.3f", vals), c("min_amp", "max_amp", "q05", "q95"))
stopifnot(vals[1] < vals[2], vals[3] < vals[4])   # protege contra ordem trocada

# Ajuste de tamanho aplicado ao objeto ggplot salvo (sem refazer os RDS)
f_fig_final <- function(p, base, ...) {
  p + ggplot2::theme(text = ggplot2::element_text(size = base), ...)
}
f_fig1 <- function(p) {
  f_fig_final(p, base = 8.5,
              strip.text = ggplot2::element_text(size = 9, face = "bold"))
}
f_fig2 <- function(p) {
  f_fig_final(p, base = 7,
              legend.position = "bottom",
              legend.key.size = grid::unit(0.35, "cm"),
              strip.text = ggplot2::element_text(size = 6, colour = "black")) +
    ggplot2::scale_x_continuous(breaks = c(0.5, 0.7, 0.9)) +
    ggplot2::labs(y = expression(log(U[i] / U[j]))) +
    ggplot2::guides(linetype = ggplot2::guide_legend(
      title.theme = ggplot2::element_text(size = 7, face = "bold")))
}
```

```{r exporta-figuras-submissao,eval=FALSE,include=FALSE}
# Rodar manualmente: arquivos separados (TIFF + PDF) para o sistema de submissao
dir.create("figuras/submissao", showWarnings = FALSE, recursive = TRUE)
save_fig <- function(p, nm, w, h) {
  ggplot2::ggsave(sprintf("figuras/submissao/%s.tiff", nm), p, width = w, height = h,
                  units = "in", dpi = DPI_FIG, compression = "lzw", bg = "white")
  ggplot2::ggsave(sprintf("figuras/submissao/%s.pdf", nm), p, width = w, height = h,
                  units = "in", device = grDevices::cairo_pdf)
}
save_fig(f_fig1(readRDS("1_to_compile_dissertacao_EM_USO/09_SI/RDS/p_interpretacao_cong.rds")),
         "Fig1", W_FIG, H_FIG1)
save_fig(f_fig2(readRDS("1_to_compile_dissertacao_EM_USO/00_Resultados/figuras/p_efeitos_porclasseCF.rds")),
         "Fig2", W_FIG, H_FIG2)
# Fig. 3: no chunk fig3display, descomente o ggsave() e use width = W_FIG, height = W_FIG
```)--")
x <- insert_after_chunk(x, "setup resultados", new_setup)

# ---- 2. vlabel: legenda da Fig. 2 e texto de Results ------------------------
r   <- chunk_range(x, "vcap p-efeitos")
blk <- x[r[1]:r[2]]
blk <- blk[!grepl("^dfrange <- readRDS|^vlabel <- sapply|^names\\(vlabel\\)", blk)]
blk <- sub("maximum = 0\\.075; minimum = -0\\.074",
           "maximum = {vlabel[['max_amp']]}; minimum = {vlabel[['min_amp']]}", blk)
blk <- sub("95th quantile = 0\\.04; 5th quantile = -0\\.038 ?\\)",
           "95th quantile = {vlabel[['q95']]}; 5th quantile = {vlabel[['q05']]})", blk)
stopifnot(any(grepl("vlabel\\[\\['max_amp", blk)), any(grepl("vlabel\\[\\['q05", blk)))
x <- c(x[seq_len(r[1] - 1)], blk, x[-seq_len(r[2])])

# Results (texto corrido, fora de chunk) -> inline R
n_before <- sum(grepl("maximum = 0\\.075", x))
stopifnot(n_before == 1)
x <- sub("maximum = 0\\.075; minimum = -0\\.074",
         "maximum = `r vlabel[['max_amp']]`; minimum = `r vlabel[['min_amp']]`", x)
x <- sub("95th quantile = 0\\.04; 5th quantile = -0\\.038\\)",
         "95th quantile = `r vlabel[['q95']]`; 5th quantile = `r vlabel[['q05']]`)", x)

# ---- 3. tamanho das figuras --------------------------------------------------
x <- replace_chunk(x, "figCong", lines_of(r"--(```{r figCong,fig.cap=vcap,fig.width=W_FIG,fig.height=H_FIG1,dpi=DPI_FIG,eval=TRUE,include=TRUE}
p <- readRDS(file = "1_to_compile_dissertacao_EM_USO/09_SI/RDS/p_interpretacao_cong.rds")
f_fig1(p)
```)--"))

x <- replace_chunk(x, "p-efeitos", lines_of(r"--(```{r p-efeitos,fig.cap=vcap,fig.width=W_FIG,fig.height=H_FIG2,dpi=DPI_FIG}
p_efeitos <- readRDS("1_to_compile_dissertacao_EM_USO/00_Resultados/figuras/p_efeitos_porclasseCF.rds")
f_fig2(p_efeitos)
```)--"))

i <- grep("^```\\{r fig3display,", x)
stopifnot(length(i) == 1, grepl("fig.width=6.8,fig.height=6.8", x[i], fixed = TRUE))
x[i] <- sub("fig.width=6.8,fig.height=6.8",
            "fig.width=W_FIG,fig.height=W_FIG,dpi=DPI_FIG", x[i], fixed = TRUE)

# ---- 4. referencias cruzadas -------------------------------------------------
# \@ref(fig:x) devolve so o numero ("1"); o prefixo "Fig." continua no texto.
# "Fig. SI1", "Fig. S4.1", "Figs. SI2.2" nao casam (exigem digito 1-3 + fronteira).
# Aplicada so ao texto corrido: linhas dentro de chunks R ficam intactas.
xref <- function(x, num, label, kind = "fig", prefix = "[Ff]ig") {
  alvo <- idx_texto(x)
  pat <- sprintf("\\b(%s)\\. %d\\b", prefix, num)
  rep <- sprintf("\\1. \\\\@ref(%s:%s)", kind, label)
  x[alvo] <- gsub(pat, rep, x[alvo], perl = TRUE)
  x
}
x <- xref(x, 1, "figCong")
x <- xref(x, 2, "p-efeitos")
x <- xref(x, 3, "fig3display")

# ---- 5. tabela (opcional) ----------------------------------------------------
if (TAB_BOOKDOWN) {
  x <- sub('^vcap <- "Model comparison', 'vcap_tab <- "Model comparison', x)
  x <- x[!grepl("set_caption(caption = vcap)", x, fixed = TRUE)]
  i <- grep("^```\\{r tab-selecaoCong,", x)
  stopifnot(length(i) == 1)
  x[i] <- sub("\\}$", ",tab.id='tab-selecaoCong',tab.cap=vcap_tab}", x[i])
  x <- xref(x, 1, "tab-selecaoCong", kind = "tab", prefix = "Tab")
}

# ---- 6. subscritos -----------------------------------------------------------
# Sintaxe de subscrito do pandoc: U~pristine~ (sem espacos entre os ~ ~).
# Vale para o docx e para o pdf. Aplicado ao texto corrido e as duas legendas
# que contem U_...: Fig. 2 e Fig. 3 (strings R dentro de chunks "vcap ...").
# Outros chunks (codigo) NAO sao tocados.
alvo <- idx_texto(x)
for (lab in c("vcap fig3", "vcap p-efeitos")) {
  r <- chunk_range(x, lab)
  alvo <- c(alvo, (r[1] + 1):(r[2] - 1))
}
alvo <- sort(unique(alvo))

subscr <- function(l) {
  l <- gsub("\\b([UM])_(pristine|fragmented|clumped|numerator|denominator|combined|friction|edge)\\b",
            "\\1~\\2~", l, perl = TRUE)
  gsub("\\bUi/Uj\\b", "U~i~/U~j~", l, perl = TRUE)
}
tmp <- subscr(x[alvo])

# "U_matrix landscape" tem espaco: no pandoc o espaco precisa de barra ("U~matrix\ landscape~").
# Como a legenda e uma string R, o arquivo deve conter DUAS barras ("\\ ").
# regmatches<- grava o texto literalmente (sem interpretar barras como o gsub faria).
m <- gregexpr("U_matrix landscape", tmp, fixed = TRUE)
regmatches(tmp, m) <- lapply(regmatches(tmp, m),
                             function(v) rep("U~matrix\\\\ landscape~", length(v)))
x[alvo] <- tmp

# Rotulo do eixo y no chunk (eval=FALSE) que recria a Fig. 2, para manter consistencia
# se os RDS forem regenerados (o render atual ja corrige via f_fig2()).
x <- sub('y="logU/U"', 'y=expression(log(U[i]/U[j]))', x, fixed = TRUE)

writeLines(x, f_out, useBytes = TRUE)

# ---- resumo ------------------------------------------------------------------
cat("Escrito:", f_out, "\n")
cat("Refs \\@ref inseridas:", sum(lengths(regmatches(x, gregexpr("\\\\@ref\\(", x)))), "\n")
cat("Subscritos: linhas com 'U_xxx' restantes em texto/legendas (esperado: nenhuma):\n")
print(grep("\\bU_[a-z]", x[alvo], value = TRUE))
cat("Linhas com 'Fig. N' manual restantes (conferir):\n")
print(grep("\\b[Ff]ig\\. [123]\\b", x, value = TRUE))

# README -----------------------------------------------------------------------
# Teste da tabela (TAB_BOOKDOWN = TRUE): se a legenda sair sem numero ou como
# "Table 1" duplicado, volte a FALSE e mantenha "Tab. 1" digitado a mao.