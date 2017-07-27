MODULE VCA

  USE VCA_INPUT_VARS
  USE VCA_VARS_GLOBAL


  USE VCA_AUX_FUNX, only: &
       los2nnn_reshape,   &
       nnn2los_reshape,   &
       set_Hcluster,      &
       search_chemical_potential


  USE VCA_IO, only:                &
       vca_get_Gcluster_matsubara, &
       vca_get_Gcluster_realaxis,  &
       vca_get_Gsystem_matsubara,  &
       vca_get_Gsystem_realaxis,   &
       vca_get_dens,               &
       vca_get_mag,                &
       vca_get_docc



  USE VCA_MAIN, only:      &
       vca_init_solver,    &
       vca_diag_cluster! ,   &
  ! vca_diag_system


END MODULE VCA
