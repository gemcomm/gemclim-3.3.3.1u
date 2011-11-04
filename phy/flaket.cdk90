! File %M% from Library %Q%
! Version %I% from %G% extracted: %H%
!------------------------------------------------------------------------------

MODULE flaket

!------------------------------------------------------------------------------
!
! Description:
!
!  The main program unit of the lake model FLake,  
!  containing most of the FLake procedures.
!  Most FLake variables and local parameters are declared.
!
!  FLake (Fresh-water Lake) is a lake model capable of predicting the surface temperature 
!  in lakes of various depth on the time scales from a few hours to a year.
!  The model is based on a two-layer parametric representation of
!  the evolving temperature profile, where the structure of the stratified layer between the
!  upper mixed layer and the basin bottom, the lake thermocline,
!  is described using the concept of self-similarity of the temperature-depth curve.
!  The concept was put forward by Kitaigorodskii and Miropolsky (1970) 
!  to describe the vertical temperature structure of the oceanic seasonal thermocline.
!  It has since been successfully used in geophysical applications.
!  The concept of self-similarity of the evolving temperature profile
!  is also used to describe the vertical structure of the thermally active upper layer 
!  of bottom sediments and of the ice and snow cover.
!
!  The lake model incorporates the heat budget equations
!  for the four layers in question, viz., snow, ice, water and bottom sediments,
!  developed with due regard for the vertically distributed character
!  of solar radiation heating.
!  The entrainment equation that incorporates the Zilitinkevich (1975) spin-up term
!  is used to compute the depth of a convectively-mixed layer. 
!  A relaxation-type equation is used
!  to compute the wind-mixed layer depth in stable and neutral stratification,
!  where a multi-limit formulation for the equilibrium mixed-layer depth
!  proposed by Zilitinkevich and Mironov (1996)
!  accounts for the effects of the earth's rotation, of the surface buoyancy flux
!  and of the static stability in the thermocline.
!  The equations for the mixed-layer depth are developed with due regard for  
!  the volumetric character of the radiation heating.
!  Simple thermodynamic arguments are invoked to develop
!  the evolution equations for the ice thickness and for the snow thickness.
!  The heat flux through the water-bottom sediment interface is computed,
!  using a parameterization proposed by Golosov et al. (1998).
!  The heat flux trough the air-water interface 
!  (or through the air-ice or air-snow interface)
!  is provided by the driving atmospheric model.
!
!  Empirical constants and parameters of the lake model
!  are estimated, using independent empirical and numerical data.
!  They should not be re-evaluated when the model is applied to a particular lake.
!  The only lake-specific parameters are the lake depth,
!  the optical characteristics of lake water,
!  the temperature at the bottom of the thermally active layer
!  of bottom sediments and the depth of that layer.
!
!  A detailed description of the lake model is given in
!  Mironov, D. V., 2005:
!  Parameterization of Lakes in Numerical Weather Prediction.
!  Part 1: Description of a Lake Model.
!  Manuscript is available from the author.
!  Dmitrii Mironov 
!  German Weather Service, Kaiserleistr. 29/35, D-63067 Offenbach am Main, Germany.
!  dmitrii.mironov@dwd.de 
!
!  Lines embraced with "!_tmp" contain temporary parts of the code.
!  Lines embraced/marked with "!_dev" may be replaced
!  as improved parameterizations are developed and tested.
!  Lines embraced/marked with "!_dm" are DM's comments
!  that may be helpful to a user.
!  Lines embraced/marked with "!_dbg" are used
!  for debugging purposes only.
!
!
! Current Code Owner: DWD, Dmitrii Mironov
!  Phone:  +49-69-8062 2705
!  Fax:    +49-69-8062 3721
!  E-mail: dmitrii.mironov@dwd.de
!
! History:
! Version    Date       Name
! ---------- ---------- ----
! 1.00       2005/11/17 Dmitrii Mironov 
!  Initial release
! !VERSION!  !DATE!     <Your name>
!  <Modification comments>
!
! Code Description:
! Language: Fortran 90.
! Software Standards: "European Standards for Writing and
! Documenting Exchangeable Fortran 90 Code".
!==============================================================================
!
! Declarations:
!
! Modules used:

USE flake_data_parameters , ONLY : &
  ireals                   , & ! KIND-type parameter for real variables
  iintegers                    ! KIND-type parameter for "normal" integer variables

!==============================================================================

IMPLICIT NONE
PRIVATE

PUBLIC :: flake_radflux, flake_driver, flake_buoypar
PUBLIC :: flake_snowdensity, flake_snowheatconduct

PUBLIC :: &
  T_mnw_p_flk, T_mnw_n_flk,       &
  T_snow_p_flk, T_snow_n_flk,     &
  T_ice_p_flk, T_ice_n_flk,       &
  T_wML_p_flk, T_wML_n_flk,       &
  T_bot_p_flk, T_bot_n_flk,       &
  T_B1_p_flk, T_B1_n_flk,         &
  h_snow_p_flk, h_snow_n_flk,     &
  h_ice_p_flk, h_ice_n_flk,       &
  h_ML_p_flk, h_ML_n_flk,         &
  H_B1_p_flk, H_B1_n_flk,         &
  C_T_p_flk, C_T_n_flk,           &
  Q_snow_flk, Q_ice_flk,          &
  Q_w_flk, Q_bot_flk,             &
  I_atm_flk, I_snow_flk,          &
  I_ice_flk, I_w_flk,             &
  I_h_flk, I_bot_flk,             &
  I_intm_0_h_flk,                 &
  I_intm_h_D_flk,                 &
  Q_star_flk,                     & 
  u_star_w_flk,                   &
  w_star_sfc_flk,                 &
  dMsnowdt_flk

!==============================================================================
!
! Declarations
!
!  The variables declared below
!  are accessible to all program units of the MODULE flake.
!  Some of them should be USEd by the driving routines that call flake routines.
!  These are basically the quantities computed by FLake.
!  All variables declared below have a suffix "flk".

!  FLake variables of type REAL

!  Temperatures at the previous time step ("p") and the updated temperatures ("n") 
REAL (KIND = ireals) ::           &
  T_mnw_p_flk, T_mnw_n_flk      , & ! Mean temperature of the water column [K] 
  T_snow_p_flk, T_snow_n_flk    , & ! Temperature at the air-snow interface [K] 
  T_ice_p_flk, T_ice_n_flk      , & ! Temperature at the snow-ice or air-ice interface [K] 
  T_wML_p_flk, T_wML_n_flk      , & ! Mixed-layer temperature [K] 
  T_bot_p_flk, T_bot_n_flk      , & ! Temperature at the water-bottom sediment interface [K] 
  T_B1_p_flk, T_B1_n_flk            ! Temperature at the bottom of the upper layer of the sediments [K] 

!  Thickness of various layers at the previous time step ("p") and the updated values ("n") 
REAL (KIND = ireals) ::           &
  h_snow_p_flk, h_snow_n_flk    , & ! Snow thickness [m]
  h_ice_p_flk, h_ice_n_flk      , & ! Ice thickness [m]
  h_ML_p_flk, h_ML_n_flk        , & ! Thickness of the mixed-layer [m] 
  H_B1_p_flk, H_B1_n_flk            ! Thickness of the upper layer of bottom sediments [m] 

!  The shape factor(s) at the previous time step ("p") and the updated value(s) ("n") 
REAL (KIND = ireals) ::           &
  C_T_p_flk, C_T_n_flk          , & ! Shape factor (thermocline)
  C_TT_flk                      , & ! Dimensionless parameter (thermocline)
  C_Q_flk                       , & ! Shape factor with respect to the heat flux (thermocline)
  C_I_flk                       , & ! Shape factor (ice)
  C_S_flk                           ! Shape factor (snow)

!  Derivatives of the shape functions
REAL (KIND = ireals) ::           &
  Phi_T_pr0_flk                 , & ! d\Phi_T(0)/d\zeta   (thermocline)
  Phi_I_pr0_flk                 , & ! d\Phi_I(0)/d\zeta_I (ice)
  Phi_I_pr1_flk                 , & ! d\Phi_I(1)/d\zeta_I (ice)
  Phi_S_pr0_flk                     ! d\Phi_S(0)/d\zeta_S (snow)

!  Heat and radiation fluxes
REAL (KIND = ireals) ::           &
  Q_snow_flk                    , & ! Heat flux through the air-snow interface [W m^{-2}]
  Q_ice_flk                     , & ! Heat flux through the snow-ice or air-ice interface [W m^{-2}]
  Q_w_flk                       , & ! Heat flux through the ice-water or air-water interface [W m^{-2}]
  Q_bot_flk                     , & ! Heat flux through the water-bottom sediment interface [W m^{-2}]
  I_atm_flk                     , & ! Radiation flux at the lower boundary of the atmosphere [W m^{-2}],
                                    ! i.e. the incident radiation flux with no regard for the surface albedo.
  I_snow_flk                    , & ! Radiation flux through the air-snow interface [W m^{-2}]
  I_ice_flk                     , & ! Radiation flux through the snow-ice or air-ice interface [W m^{-2}]
  I_w_flk                       , & ! Radiation flux through the ice-water or air-water interface [W m^{-2}]
  I_h_flk                       , & ! Radiation flux through the mixed-layer-thermocline interface [W m^{-2}]
  I_bot_flk                     , & ! Radiation flux through the water-bottom sediment interface [W m^{-2}]
  I_intm_0_h_flk                , & ! Mean radiation flux over the mixed layer [W m^{-1}]
  I_intm_h_D_flk                , & ! Mean radiation flux over the thermocline [W m^{-1}]
  Q_star_flk                        ! A generalized heat flux scale [W m^{-2}]

!  Velocity scales
REAL (KIND = ireals) ::           &
  u_star_w_flk                  , & ! Friction velocity in the surface layer of lake water [m s^{-1}]
  w_star_sfc_flk                    ! Convective velocity scale, 
                                    ! using a generalized heat flux scale [m s^{-1}]

!  The rate of snow accumulation
REAL (KIND = ireals) ::           &
  dMsnowdt_flk                      ! The rate of snow accumulation [kg m^{-2} s^{-1}]

!$omp threadprivate ( &
!$omp T_mnw_p_flk, T_mnw_n_flk,   &
!$omp T_snow_p_flk, T_snow_n_flk, &
!$omp T_ice_p_flk, T_ice_n_flk,   &
!$omp T_wML_p_flk, T_wML_n_flk,   &
!$omp T_bot_p_flk, T_bot_n_flk,   &
!$omp T_B1_p_flk, T_B1_n_flk,     &
!$omp h_snow_p_flk, h_snow_n_flk, &
!$omp h_ice_p_flk, h_ice_n_flk,   &
!$omp h_ML_p_flk, h_ML_n_flk,     &
!$omp H_B1_p_flk, H_B1_n_flk,     &
!$omp C_T_p_flk, C_T_n_flk,       &
!$omp Q_snow_flk, Q_ice_flk,      &
!$omp Q_w_flk, Q_bot_flk,         &
!$omp I_atm_flk, I_snow_flk,      &
!$omp I_ice_flk, I_w_flk,         &
!$omp I_h_flk, I_bot_flk,         &
!$omp I_intm_0_h_flk,             &
!$omp I_intm_h_D_flk,             &
!$omp Q_star_flk,                 & 
!$omp u_star_w_flk,               &
!$omp w_star_sfc_flk,             &
!$omp dMsnowdt_flk )

!==============================================================================
! Procedures 
!==============================================================================

CONTAINS

!==============================================================================
!  The codes of the FLake procedures are stored in separate "*.cdk" files
!  and are included below.
!------------------------------------------------------------------------------

!==============================================================================
include 'flake_radflux.h'

!==============================================================================

!==============================================================================
include 'flake_driver.h'

!==============================================================================

!==============================================================================
include 'flake_buoypar.h'

!==============================================================================

!==============================================================================
include 'flake_snowdensity.h'

!==============================================================================

!==============================================================================
include 'flake_snowheatconduct.h'

!==============================================================================

END MODULE flaket

