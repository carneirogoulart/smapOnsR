test_that("Testa a funcao de calibracao", {
  nome2 <- "pimentalt"
  modelo <- new_modelo_smap_ons(parametros[Nome == nome2], postos_plu[nome == nome2])
  kt_max <- max(which(modelo$kt[3:1] > 0)) - 1
  kt_min <- max(which(modelo$kt[3:63] > 0)) - 1

  EbInic <- 2300
  TuInic <- 0.50
  Supin <- 2000

  normal_climatologica <- historico_etp_NC[nome == nome2]
  precipitacao <- historico_precipitacao[posto %in% postos_plu[nome == nome2, posto]]
  precipitacao_ponderada <- ponderacao_espacial(precipitacao, postos_plu[nome == nome2])
  data_inicio_objetivo <- as.Date("2011-01-01")
  data_fim_objetivo <- as.Date("2021-12-31") - kt_max
  evapotranspiracao <- transforma_NC_serie(precipitacao_ponderada[data >= min(data) + kt_min & data <= data_fim_objetivo], normal_climatologica)
  vazao <- historico_vazao[data >= data_inicio_objetivo & data <= data_fim_objetivo & posto == nome2]

  area <- attributes(modelo)$area
  vetor_modelo <- unlist(modelo)
  vetor_modelo <- c(vetor_modelo[1:11], vetor_modelo[75:77])
  numero_postos_plu <- postos_plu[nome == nome2, length(unique(posto))]
  vetor_modelo[15:(14 + numero_postos_plu * 2)] <- 0.01

  limite_inferior <- vetor_modelo * 0.001
  limite_superior <- vetor_modelo * 100
  limite_inferior[8] <- 0.9999999
  limite_superior[8] <- 1.0000001
  limite_inferior[12] <- 0.8
  limite_superior[12] <- 1.2
  limite_inferior[13] <- 0.8
  limite_superior[13] <- 1.2
  limite_inferior[14] <- 0.8
  limite_superior[14] <- 1.2
  limite_inferior[15:(14 + numero_postos_plu * 2)] <- 0.000000001
  limite_superior[15:(14 + numero_postos_plu * 2)] <- 50

  fo <- funcao_objetivo_calibracao(vetor_modelo, kt_min, kt_max, area, EbInic, TuInic, Supin, precipitacao,
      evapotranspiracao, vazao, data_inicio_objetivo, data_fim_objetivo, postos_plu[nome == nome2])

  par <- calibracao(vetor_modelo,  kt_min, kt_max, area, EbInic, TuInic, Supin, precipitacao,
      evapotranspiracao, vazao, data_inicio_objetivo, data_fim_objetivo,
      limite_inferior, limite_superior, postos_plu[nome == nome2])

  expect_equal((par$value < fo), TRUE)
})
