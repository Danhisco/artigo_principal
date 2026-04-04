# README
- O corpo principal do texto está em 1_tocomp/compilacao_integral.Rmd. Esse arquivo contempla:
-- Title, authors, abstract, introduction, methodology, results, discussion, conclussion, references 
- O apêndice 1, que trata de detalhes metodologicos, está em 1_tocomp/SI.Rmd, que contempla:
-- The features from the pre-selected forest plots; Geoprocessing; Dispersal kernel
- O apêndice 2, sobre a estimativa da escala suficiente para cada cenário de dispersão, está em 1_tocomp/SoE.Rmd, que contempla:
-- Context; Largest landscape scale approximates infinite-landscape estimates; Spatial-scale dependence of U estimates across dispersal scenarios; Conclusion

Nesses dois arquivos .Rmd (SoE e SI) há os códigos necessários para realizar as tarefas descritas neles. 
Isso pq o processo inicial de iniciação do artigo é: 
i) pre selecionar os inventários florestais no TreeCo que tenham correspondência (ano do dado) na coleção 6 do mapbiomas
ii) estimar os cenários de dispersão em função da densidade de indivíduos na parcela florestal
iii) estabelecer correspondência entre cenários de limitaçaõ de dispersão e escalas espaciais da paisagem ao redor

Uma vez completada essas etapas, então é necessário preparar os mapas de cobertura florestal:
- Em 0_apendices/mapas_p_MNEE/mapas_p_MNEE.Rmd há a preparação das paisagens e parte do trabalho apresentado no apêndice 1 sobre geoprocessamento
Ese apêndice se inicia cortando os mapas de cobertura florestal da coleção 6 do mapiomas em equivalentes de lado 16km.
Então é feito o ajuste da densidade de pixels para ser equivalente à densidade de indivíduos na parcela amostrada.
E é criada a parcela amostrada, pressuposto forma quadrada, no centro da paisagem.
Uma vez que a maior escala espacial da paisagem foi de lado 4 km, as paisagens salvas em dados/simulacao/*.txt são dessa extensão.

Na sequência, as simulações são executadas em '1_torep/simulacoes.R', que:
1) realiza as simulações (de U e SADs para cada combinação de parcela, k e land) [f_simMNEE]
2) compara as SADs simuladas e observadas
3) sumariza e armazena nas pastas: 
    a) "dados/csv_SoE/df_KSrep.csv" # comparação bruta das SADs
    b) "dados/csv_SoE/taxaU/df_U.csv" # valor bruto da taxa U
    c) "dados/csv_SoE/taxaU/df_contrastes.csv" # valor dos efeitos da paisagem FLF, FPS e LFC 
Nesse arquivo é necessário hardcoding as decisões da escala adequada por cenário de dispersão,
pois não decorrem diretamente da abordagem heurística e sim da interpretação desses resultados.

No apêndice 3 01_tocomp/CongSAD.Rmd há as análises finais para obter a figura 1 do texto principal.
Nesse apêndice é feito toda a análise a partir de dados/csv_SoE/df_KSrep.csv.
Aqui é feito a análise de autocorrelação espacial global (I de Moran),
a seleção de modelos e predição,
as implicações e a informação de suporte do apêndice.



