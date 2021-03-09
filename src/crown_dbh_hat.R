
crown_dbh.hat.position <-  brm(DBH ~ s(CROWN_DIAMETER_WIDE_hat) + s(HEIGHT_HAT)
                      # population level features
                      + s(H_VAR) 
                      + s(STATURE)
                      + s(PLOT_MAX_WIDTH)
                      # scale nesting
                      + ((1 + HEIGHT_HAT) | domainID / siteID / 
                           FOREST_TYPE_DESCR / canopyPosition) 
                      #taxon level features
                      + ((1 + CROWN_DIAMETER_WIDE_hat 
                          + HEIGHT_HAT  
                          + STATURE
                          + PLOT_MAX_WIDTH
                          + H_VAR 
                          + HEIGHT_HAT:CROWN_DIAMETER_WIDE_hat  
                          + PLOT_MAX_WIDTH:HEIGHT_HAT  
                          + STATURE:HEIGHT_HAT
                          + H_VAR:HEIGHT_HAT
                      ) |gr(taxonID, cov = PHENOLOGY)) 
                      # landuse level feature
                      ,
                      # define data and grouping correlation structures
                      data = train, family = exgaussian(),
                      data2 = list(PHENOLOGY = PHENOLOGY), 
                      #other parameters
                      chains=2, cores =6, backend = "cmdstanr", threads = threading(3))

crown_dbh.hat.quantile <-  brm(DBH ~ s(CROWN_DIAMETER_WIDE_hat) + s(HEIGHT_HAT)
                          # population level features
                          + s(H_VAR) 
                          + s(STATURE)
                          + s(PLOT_MAX_WIDTH)
                          # scale nesting
                          + ((1 + HEIGHT_HAT) | domainID / siteID) 
                          + (1 | FOREST_TYPE_DESCR)
                          + (1 | canopyPosition) 
                          #taxon level features
                          + ((1 + CROWN_DIAMETER_WIDE_hat 
                              + HEIGHT_HAT  
                              + STATURE
                              + PLOT_MAX_WIDTH
                              + H_VAR 
                              + HEIGHT_HAT:CROWN_DIAMETER_WIDE_hat  
                              + PLOT_MAX_WIDTH:HEIGHT_HAT  
                              + STATURE:HEIGHT_HAT
                              + H_VAR:HEIGHT_HAT
                          ) |gr(taxonID, cov = PHENOLOGY)) 
                          # landuse level feature
                          ,
                          # define data and grouping correlation structures
                          data = train, family = student(),
                          data2 = list(PHENOLOGY = PHENOLOGY), 
                          #other parameters
                          chains=2, cores =6, backend = "cmdstanr", threads = threading(3))
