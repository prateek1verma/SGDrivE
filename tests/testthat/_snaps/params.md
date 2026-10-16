# print() shows groups and marks changed values

    Code
      print(sg_params(adult_substages = 8))
    Output
      <sg_params> preset "default": 52 parameters, 1 differ from the defaults (*)
      
      temperature
          dev_temp_min                  14.5                 degC                  GP2014
          dev_temp_decline              40.5                 degC                  GP2014
          dev_rate_slope                0.00065              h^-1 degC^-1          GP2014
          dev_decline_steepness         1.75                 -                     GP2014
          dev_rate_floor                0                    -                     GP2014
          optimal_survival_temp         27.2                 degC                  GP2014
          oviposition_temp_max          43.5                 degC                  GP2014
          egg_temp_source               "ambient"            -                     adapted
          host_body_temp                38.5                 degC                  GP2014
      
      development
          dd_scale_preoviposition       61                   dd                    GP2014
          egg_duration_dd               10                   dd                    GP2014
          larva_duration_days           6                    d                     GP2014
          pupa_duration_dd              95                   dd                    GP2014
          adult_lifespan_dd             202.8                dd                    GP2014
      
      stages
          egg_substages                 3                    -                     adapted
          larva_substages               6                    -                     adapted
          pupa_substages                10                   -                     adapted
        * adult_substages               8                    -                     adapted
      
      thresholds
          oviposition_age_dd            61                   dd                    GP2014
          female_receptive_age_dd       0                    dd                    adapted
          male_mature_age_dd            0                    dd                    adapted
      
      mortality
          adult_mortality_curvature     0.00036              degC^-2               GP2014
          adult_mortality_min           0.0035               -                     GP2014
          egg_failure_curvature         0.00375              degC^-2               GP2014
          egg_failure_min               0.3                  -                     GP2014
          larva_survival_daily          0.95                 -                     GP2014
          female_background_mortality   0                    d^-1                  adapted
          male_background_mortality     0                    d^-1                  adapted
          larval_crowding_mode          "none"               -                     adapted
          larval_crowding_coef          0                    d^-1 wound larva^-1   adapted
          pupa_drought_mortality        0                    d^-1                  adapted
          drought_free_rainfall         450                  mm yr^-1              GP2014
          survival_floor                1e-05                -                     GP2014
      
      reproduction
          eggs_per_day_max              90                   eggs female^-1 d^-1   GP2014
          wound_search_rate             0.001                dd^-1                 GP2014
          wound_density                 5                    wounds km^-2          GP2014
          wound_finding_max             1                    -                     adapted
          wound_competitors             "searching_females"  -                     adapted
          oviposition_mode              "batches"            -                     adapted
          clutch_size                   NA                   eggs                  adapted
      
      mating
          mating_rate                   0.693147             d^-1                  G2019
          male_density_half_mating      0                    males km^-2           adapted
          mating_mode                   "rate"               -                     MGDrivE2
      
      hosts
          wound_model                   "static"             -                     adapted
          host_density                  NA                   animals km^-2         adapted
          wounded_fraction              NA                   -                     adapted
          wound_healing_days            7                    d                     adapted
          infestation_days              6                    d                     adapted
          infestation_substages         6                    -                     adapted
          infested_wound_attractiveness 1                    -                     adapted
          treatment_rate                0                    d^-1                  adapted
      
      numerical
          extinction_threshold          1                    individuals           adapted

