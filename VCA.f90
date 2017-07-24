MODULE VCA

  USE VCA_INPUT_VARS
  USE VCA_VARS_GLOBAL


  USE VCA_AUX_FUNX, only:                       &
       lso2nnn_reshape,                         &
       nnn2lso_reshape,                         &
       set_Hcluster,                            &
       search_chemical_potential


  USE VCA_IO, only:            &
       vca_get_gimp_matsubara, &
       vca_get_gimp_realaxis,  &
       vca_get_dens,           &
       vca_get_mag,            &
       vca_get_docc,           &
       vca_print_impG,         &
       vca_read_impG



  USE VCA_MAIN, only:     &
       vca_init_solver,   &
       vca_diag


END MODULE VCA
