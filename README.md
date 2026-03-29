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
ii) estabelecer correspondência entre cenários de limitaçaõ de dispersão e escalas espaciais da paisagem ao redor

Uma vez completada essas etapas, então é necessário preparar os mapas de cobertura florestal:
- Em 0_apendices/mapas_p_MNEE/mapas_p_MNEE.Rmd há a preparação das paisagens e parte do trabalho apresentado no apêndice 1 sobre geoprocessamento
Ese apêndice se inicia cortando os mapas de cobertura florestal da coleção 6 do mapiomas em equivalentes de lado 16km.
Então é feito o ajuste da densidade de pixels para ser equivalente à densidade de indivíduos na parcela amostrada.
E é criada a parcela amostrada, pressuposto forma quadrada, no centro da paisagem.
Uma vez que a maior escala espacial da paisagem foi de lado 4 km, as paisagens salvas em dados/simulacao/*.txt são dessa extensão.

Na sequência, as simulações são executadas em:
- 0_apendices/simulacao/Simulacoes.Rmd