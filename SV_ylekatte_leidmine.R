
library(arrow)
library(tidyverse)
library(data.table)

###Laeb sisse SV ja geeni andmed

geenid <- read_csv("gencode_v49.protein_coding_genes_1M_proximity.csv", col_names = FALSE)
colnames(geenid) = c("kromosoom", "algus", "l6pp", "ymbrus_algus", "ymbrus_l6pp", "id", "nimi", "DNA_ahel", "tyyp")
geenid = geenid %>% filter(!(kromosoom %in% c("chrX", "chrY", "chrM" )))

save(geenid, file="geenid.Rda") ## andmete visualiseerimiseks hiljem

SVd = read.table(file = 'SV_filtreeritud_uus.tsv', sep = '\t', header = F)
colnames(SVd) = c("kromosoom", "algus", "l6pp", "tyyp", "tyyp2", "meetod", "pikkus", "yld_sagedus", "eu_sagedus")
SVd = SVd %>% filter(!(kromosoom %in% c("chrX", "chrY", "chrM" )))



## leiab ülekatte
setDT(SVd)
setDT(geenid)

setkey(SVd, kromosoom, algus, l6pp)
setkey(geenid, kromosoom, algus, l6pp)


hits <- foverlaps(
  SVd,
  geenid,
  type = "any",
  nomatch = NULL,
  which = TRUE
)

SVd[, geenis := FALSE]
SVd[unique(hits$xid), geenis := TRUE] ### paneb kirja geeniga ülekatte





process <- function(eqtl_file, SVd, pikkus = 100) { ## leiab eqtl fails olevas olud SV-dega, +-100 eqtl kohale
  
  
  eqtl <- read_parquet(eqtl_file)
  
  eqtl <- eqtl %>% ## lisab +-100 (või valitud pikkused)
    mutate(
      tmp = str_split(variant_id, "_", simplify = TRUE),
      kromosoom = tmp[,1],
      pos = as.numeric(tmp[,2]),
      algus = pos - pikkus,
      l6pp = pos + pikkus
    ) %>%
    select(-tmp)
  

  setDT(eqtl)
  setDT(SVd)
  
  setkey(eqtl, kromosoom, algus, l6pp)
  setkey(SVd, kromosoom, algus, l6pp)
  
  hits <- foverlaps(
    SVd,
    eqtl,
    type = "any",
    nomatch = NULL,
    which = TRUE
  )
  
  out <- rep(FALSE, nrow(SVd)) ## paneb kirja ülekatte
  out[hits$xid] <- TRUE
  
  return(out)
}

eqtl_files <- list.files( ## loeb kõik failid kaustast
  "dataEQTL",
  pattern = "\\.v10\\.eQTLs\\.signif_pairs\\.parquet$",
  full.names = TRUE
)



for (f in eqtl_files) { ## leiab iga koe ülekatted
  
  koed <- sub(
    "\\.v10\\.eQTLs\\.signif_pairs\\.parquet$", "" ,basename(f))
  
  SVd[[paste0("eQTL_", koed)]] <- process(f, SVd) ## paneb kirja, kas antud koes oli SV-l eQTL ülekatet
}




SVd_valmis <- SVd %>%
  mutate(
    eQTL_any = if_any(
      starts_with("eQTL_"),
      ~ . == TRUE
    )
  ) %>%
  select(
    kromosoom, algus, l6pp, tyyp, pikkus,
    yld_sagedus, eu_sagedus, geenis,
    eQTL_any, 
    tyyp2, 
    meetod
  ) ### valib vahalikud tunnused, eQTL loeb kõik, millel on vähemalt üks ülekate.


save(SVd_valmis, file="andmed.Rda") ## salvestab faili


