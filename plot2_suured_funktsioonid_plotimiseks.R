

plot_mudel2 = function(data, mudel, facet = FALSE, ylim = c(0,1)){ ## plotib mudeli ennustused koos andmetest nähtava osakaaluga
  
  bp <- log10(5000)
  
  newdat <- expand.grid(
    log10_pikkus = seq(min(data$log10_pikkus),
                       max(data$log10_pikkus),
                       length.out = 200),
    geenis = levels(data$geenis),
    eQTL  = levels(data$eQTL)
  ) %>%
    mutate(
      after_bp   = as.numeric(log10_pikkus > bp),
      combo = interaction(geenis, eQTL, sep = " × ") ## loob kõik võimalikud väärtused, millele saab ennustada
    )
  
  
  pred_link <- predict( ## ennustab mudelist
    mudel,
    newdata = newdat,
    type = "link",
    se.fit = TRUE
  )
  
  newdat <- newdat %>%
    mutate(
      fit_link = pred_link$fit,
      se_link  = pred_link$se.fit,
      fit      = plogis(fit_link),
      lower    = plogis(fit_link - 1.96 * se_link),
      upper    = plogis(fit_link + 1.96 * se_link),
      eQTL  = factor(ifelse(as.character(eQTL) == "TRUE", "JAH", "EI    ")), ## uus sõnastus joonisel
      geenis  = factor(ifelse(as.character(geenis) == "TRUE", "JAH", "EI    ")), 
      combo = interaction(geenis, eQTL, sep = " + ")                       
    ) ##seob andmed
  
  
  binned <- data %>%  ### paneb juurde andmetest nähtud osakaalud
    mutate(
      log10_bin = round(log10_pikkus, 1),
      eQTL  = factor(ifelse(as.character(eQTL) == "TRUE", "JAH", "EI    ")), 
      geenis  = factor(ifelse(as.character(geenis) == "TRUE", "JAH", "EI    ")), 
      combo = interaction(geenis, eQTL, sep = " + ")
    ) %>%
    group_by(log10_bin, geenis, eQTL, combo) %>%
    summarise(
      mean_haruldane = mean(haruldane, na.rm = TRUE),
      n = n(),
      se = sqrt(mean_haruldane * (1 - mean_haruldane) / n),
      lower = pmax(mean_haruldane - 1.96 * se, 0),
      upper = pmin(mean_haruldane + 1.96 * se, 1),
      .groups = "drop"
    ) 
  
  
  if (!facet) { ## plotib kui pole facet
    
    p <- ggplot(newdat,
                aes(x = 10^log10_pikkus,
                    color = combo,
                    linetype = combo,
                    fill = combo)) +
      
      geom_ribbon(aes(ymin = lower, ymax = upper),
                  alpha = 0.2,
                  color = NA) +
      
      geom_line(aes(y = fit), size = 1.1) +
      
      geom_errorbar(
        data = binned,
        aes(x = 10^log10_bin,
            ymin = lower,
            ymax = upper,
            color = combo),
        inherit.aes = FALSE,
        width = 0,
        alpha = 0.6
      ) +
      
      geom_point(
        data = binned,
        aes(x = 10^log10_bin,
            y = mean_haruldane,
            color = combo),
        inherit.aes = FALSE,
        size = 1.5,
        alpha = 0.9
      ) +
      
      labs(
        color = "   geenis + eQTL",
        fill  = "   geenis + eQTL",
        linetype = "   geenis + eQTL"
      )
    
  } else { ## facet plot
    
    p <- ggplot(newdat,
                aes(x = 10^log10_pikkus,
                    y = fit)) +
      
      geom_ribbon(aes(ymin = lower, ymax = upper),
                  alpha = 0.2,
                  fill = "grey70") +
      
      geom_line(size = 1.1) +
      
      geom_errorbar(
        data = binned,
        aes(x = 10^log10_bin,
            ymin = lower,
            ymax = upper),
        inherit.aes = FALSE,
        width = 0,
        alpha = 0.6
      ) +
      
      geom_point(
        data = binned,
        aes(x = 10^log10_bin,
            y = mean_haruldane),
        inherit.aes = FALSE,
        size = 1.5,
        alpha = 0.9
      ) +
      
      facet_wrap(~ geenis + eQTL) +
      
      guides(color = "none",
             fill = "none",
             linetype = "none")
  }
  
  
  
  p <- p +
    geom_vline(xintercept = 5000, linetype = "dashed") +
    scale_x_log10(
      breaks = scales::trans_breaks("log10", function(x) 10^x),
      labels = scales::trans_format("log10", scales::math_format(10^.x))
    ) +
    coord_cartesian(ylim = ylim) +   
    labs(
      x = "SV pikkus",
      y = "Haruldaste SV-de tõenäosus"
    ) +
    theme_classic()
  
  return(p)
}


plot_mudel2_ME = function(data, mudel, facet = FALSE,ylim = c(0,1)){#### identne ülemise funtsiooniga kuid on lisatu mobiilsete elementide jaoks normaaljaotuse tihedusfunktsioon
  
  
  bp <- log10(5000)
  
  
  newdat <- expand.grid(
    log10_pikkus = seq(min(data$log10_pikkus),
                       max(data$log10_pikkus),
                       length.out = 200),
    geenis = levels(data$geenis),
    eQTL  = levels(data$eQTL)
  ) %>%
    mutate(
      after_bp = as.numeric(log10_pikkus > bp),
      combo = interaction(geenis, eQTL, sep = " × "),
      ME_N = dnorm(log10_pikkus, 2.420, 0.1010715)
    )
  
  
  pred_link <- predict(
    mudel,
    newdata = newdat,
    type = "link",
    se.fit = TRUE
  )
  
  newdat <- newdat %>%
    mutate(
      fit_link = pred_link$fit,
      se_link  = pred_link$se.fit,
      fit      = plogis(fit_link),
      lower    = plogis(fit_link - 1.96 * se_link),
      upper    = plogis(fit_link + 1.96 * se_link),
      eQTL  = factor(ifelse(as.character(eQTL) == "TRUE", "JAH", "EI    ")), 
      geenis  = factor(ifelse(as.character(geenis) == "TRUE", "JAH", "EI    ")), 
      combo = interaction(geenis, eQTL, sep = " + ")
    )
  
  
  binned <- data %>%
    mutate(
      log10_bin = round(log10_pikkus, 1),
      eQTL  = factor(ifelse(as.character(eQTL) == "TRUE", "JAH", "EI    ")), 
      geenis  = factor(ifelse(as.character(geenis) == "TRUE", "JAH", "EI    ")), 
      combo = interaction(geenis, eQTL, sep = " + ")
    ) %>%
    group_by(log10_bin, geenis, eQTL, combo) %>%
    summarise(
      mean_haruldane = mean(haruldane, na.rm = TRUE),
      n = n(),
      se = sqrt(mean_haruldane * (1 - mean_haruldane) / n),
      lower = pmax(mean_haruldane - 1.96 * se, 0),
      upper = pmin(mean_haruldane + 1.96 * se, 1),
      .groups = "drop"
    )
  
  
  p <- ggplot(newdat,
              aes(x = 10^log10_pikkus,
                  y = fit,
                  group = combo))
  
  if (!facet) {
    
    p <- p +
      geom_ribbon(aes(ymin = lower,
                      ymax = upper,
                      fill = combo),
                  alpha = 0.2,
                  colour = NA) +
      geom_line(aes(color = combo,
                    linetype = combo),
                size = 1.1) +
      geom_errorbar(
        data = binned,
        aes(x = 10^log10_bin,
            ymin = lower,
            ymax = upper,
            color = combo),
        inherit.aes = FALSE,
        width = 0,
        alpha = 0.6
      ) +
      geom_point(
        data = binned,
        aes(x = 10^log10_bin,
            y = mean_haruldane,
            color = combo),
        inherit.aes = FALSE,
        size = 1.5,
        alpha = 0.9
      )+
      
      labs(
        color = "   geenis + eQTL",
        fill  = "   geenis + eQTL",
        linetype = "   geenis + eQTL"
      )
    
  } else {
    
    p <- p +
      geom_ribbon(aes(ymin = lower, ymax = upper),
                  alpha = 0.2,
                  fill = "grey70",
                  colour = NA) +
      geom_line(size = 1.1) +
      geom_errorbar(
        data = binned,
        aes(x = 10^log10_bin,
            ymin = lower,
            ymax = upper),
        inherit.aes = FALSE,
        width = 0,
        alpha = 0.6
      ) +
      geom_point(
        data = binned,
        aes(x = 10^log10_bin,
            y = mean_haruldane),
        inherit.aes = FALSE,
        size = 1.5,
        alpha = 0.9
      ) +
      facet_wrap(~ geenis + eQTL)
  }
  
  p +
    geom_vline(xintercept = 5000, linetype = "dashed") +
    scale_x_log10(
      breaks = scales::trans_breaks("log10", function(x) 10^x),
      labels = scales::trans_format("log10", scales::math_format(10^.x))
    ) +
    coord_cartesian(ylim = ylim) +  
    labs(x = "SV pikkus",
         y = "Haruldaste SV-de tõenäosus") +
    theme_classic()
}


############
#Järgmised kaks funktsiooni on riskide suhte kujutamiseks
##############


plot_mudel2_RR_del = function(data, mudel, facet = FALSE, max_rr = 2, min_rr = 0.5){
  ### leiab deletsioonidel riskidesuhte
  
  bp <- log10(5000)
  
  
  newdat <- expand.grid( ## teeb valmis võimalike väärtuste hulga
    log10_pikkus = seq(min(data$log10_pikkus),
                       max(data$log10_pikkus),
                       length.out = 200),
    geenis = levels(data$geenis),
    eQTL  = levels(data$eQTL)
  ) %>%
    mutate(
      after_bp = as.numeric(log10_pikkus > bp),
      combo    = interaction(geenis, eQTL, sep = " + ")
    )
  
  
  pred_link <- predict( ## ennsutab
    mudel,
    newdata = newdat,
    type = "link",
    se.fit = TRUE
  )
  
  newdat <- newdat %>% ## lisab ennustus andmed
    mutate(
      fit_link = pred_link$fit,
      se_link  = pred_link$se.fit
    )
  
  
  
  #### võtab baastaseme ennsutused
  ref_df <- newdat %>%
    filter(
      geenis == levels(data$geenis)[1],
      eQTL   == levels(data$eQTL)[1]
    ) %>%
    select(log10_pikkus, ref_link = fit_link, ref_se = se_link)
  
  
  newdat <- newdat %>%
    left_join(ref_df, by = "log10_pikkus") %>% ## lisab kõigi tasemete andmed juurde
    mutate(
      log_or  = fit_link - ref_link, ## arvutab sanside suhte
      OR      = exp(log_or),
      se_diff = sqrt(se_link^2 + ref_se^2),
      se_diff = if_else(
        geenis == levels(data$geenis)[1] & eQTL == levels(data$eQTL)[1],
        0, se_diff
      ),
      p_ref   = plogis(ref_link),
      OR_lo   = exp(log_or - 1.96 * se_diff), 
      OR_hi   = exp(log_or + 1.96 * se_diff),
      RR      = OR      / ((1 - p_ref) + (p_ref * OR)), ### leiab riskidesuhte
      lower   = OR_lo   / ((1 - p_ref) + (p_ref * OR_lo)), 
      upper   = OR_hi   / ((1 - p_ref) + (p_ref * OR_hi)),
      eQTL  = factor(ifelse(as.character(eQTL) == "TRUE", "JAH", "EI    ")), 
      geenis  = factor(ifelse(as.character(geenis) == "TRUE", "JAH", "EI    ")),
      combo    = interaction(geenis, eQTL, sep = " + ")
    )
  
  
  newdat <- newdat %>% ## joonisel ei lub aminna üle piiride
    mutate(
      RR    = pmax(pmin(RR,    max_rr), min_rr),
      lower = pmax(pmin(lower, max_rr), min_rr),
      upper = pmax(pmin(upper, max_rr), min_rr)
    )
   
  if (!facet) { ## plotib
    
    p <- ggplot(newdat,
                aes(x = 10^log10_pikkus, y = RR, color = combo, fill = combo, linetype = combo)) +
      
      geom_ribbon(aes(ymin = lower, ymax = upper),
                  alpha = 0.2,
                  color = NA) +
      
      geom_line(size = 1.1)
    
  } else {
    
    p <- ggplot(newdat,
                aes(x = 10^log10_pikkus, y = RR)) +
      
      geom_ribbon(aes(ymin = lower, ymax = upper),
                  alpha = 0.2,
                  fill = "grey70") +
      
      geom_line(size = 1.1) +
      
      facet_wrap(~ geenis + eQTL)
  }
  
  
  p +
    geom_hline(yintercept = 1, linetype = "dotted") +
    geom_vline(xintercept = 5000, linetype = "dashed") +
    scale_x_log10(
      breaks = scales::trans_breaks("log10", function(x) 10^x),
      labels = scales::trans_format("log10", scales::math_format(10^.x))
    ) +
    scale_y_continuous(
      limits = c(min_rr, max_rr),
      oob = scales::squish
    ) +
    labs(
      x = "SV pikkus",
      y = "Riskide suhe võrreldes grupiga Ei + Ei",
      color = "   geenis + eQTL",
      fill  = "   geenis + eQTL",
      linetype = "   geenis + eQTL"
    ) +
    theme_classic()
}


plot_mudel2_ME_RR = function(data, mudel, facet = FALSE, max_rr = 2, min_rr=0.5){
  #Insertsioonidel identne ülemise funktsiooniga, ainult ME_N on juures
  
  
  
  bp <- log10(5000)
  
  newdat <- expand.grid(
    log10_pikkus = seq(min(data$log10_pikkus),
                       max(data$log10_pikkus),
                       length.out = 200),
    geenis = levels(data$geenis),
    eQTL   = levels(data$eQTL)
  ) %>%
    mutate(
      after_bp = as.numeric(log10_pikkus > bp),
      combo = interaction(geenis, eQTL, sep = " × "),
      ME_N = dnorm(log10_pikkus, 2.420, 0.1010715)
    )
  
  pred_link <- predict(
    mudel,
    newdata = newdat,
    type = "link",
    se.fit = TRUE
  )
  
  newdat <- newdat %>%
    mutate(
      fit_link = pred_link$fit,
      se_link  = pred_link$se.fit
    )
  
  
  
  ref_df <- newdat %>%
    filter(
      geenis == levels(data$geenis)[1],
      eQTL   == levels(data$eQTL)[1]
    ) %>%
    select(log10_pikkus, ref_link = fit_link, ref_se = se_link)
  
  newdat <- newdat %>%
    left_join(ref_df, by = "log10_pikkus") %>%
    mutate(
      log_or  = fit_link - ref_link,
      OR      = exp(log_or),
      se_diff = sqrt(se_link^2 + ref_se^2),
      se_diff = if_else(
        geenis == levels(data$geenis)[1] & eQTL == levels(data$eQTL)[1],
        0, se_diff
      ),
      p_ref   = plogis(ref_link),
      OR_lo   = exp(log_or - 1.96 * se_diff),
      OR_hi   = exp(log_or + 1.96 * se_diff),
      RR      = OR      / ((1 - p_ref) + (p_ref * OR)),
      lower   = OR_lo   / ((1 - p_ref) + (p_ref * OR_lo)),
      upper   = OR_hi   / ((1 - p_ref) + (p_ref * OR_hi)),
      eQTL  = factor(ifelse(as.character(eQTL) == "TRUE", "JAH", "EI    ")), 
      geenis  = factor(ifelse(as.character(geenis) == "TRUE", "JAH", "EI    ")),
      combo    = interaction(geenis, eQTL, sep = " + ")
    )
  

  newdat <- newdat %>%
    mutate(
      RR    = pmax(pmin(RR,    max_rr), min_rr),
      lower = pmax(pmin(lower, max_rr), min_rr),
      upper = pmax(pmin(upper, max_rr), min_rr)
    )
  
  p <- ggplot(newdat,
              aes(x = 10^log10_pikkus,
                  y = RR,
                  group = combo))
  
  if (!facet) {
    
    p <- ggplot(
      newdat,
      aes(
        x = 10^log10_pikkus,
        y = RR,
        color = combo,
        fill = combo,
        linetype = combo
      )
    ) +
      
      geom_ribbon(
        aes(ymin = lower, ymax = upper),
        alpha = 0.2,
        colour = NA
      ) +
      
      geom_line(size = 1.1)
    
  } else {
    
    p <- ggplot(
      newdat,
      aes(
        x = 10^log10_pikkus,
        y = RR
      )
    ) +
      
      geom_ribbon(
        aes(ymin = lower, ymax = upper),
        alpha = 0.2,
        fill = "grey70",
        colour = NA
      ) +
      
      geom_line(size = 1.1) +
      
      facet_wrap(~ geenis + eQTL)
  }
  
  
  p +
    geom_hline(yintercept = 1, linetype = "dotted") +
    geom_vline(xintercept = 5000, linetype = "dashed") +
    scale_x_log10(
      breaks = scales::trans_breaks("log10", function(x) 10^x),
      labels = scales::trans_format("log10", scales::math_format(10^.x))
    ) +
    coord_cartesian(ylim = c(min_rr, max_rr)) +
    labs(
      x = "SV pikkus",
      y = "Riskide suhe võrreldes grupiga Ei + Ei",
      color = "   geenis + eQTL",
      fill = "   geenis + eQTL",
      linetype = "   geenis + eQTL"
    ) +
    theme_classic()
}

