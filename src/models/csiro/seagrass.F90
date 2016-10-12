#include "fabm_driver.h"

module csiro_seagrass

   ! This is a FABM implementation of the CSIRO seagrass model described in
   ! Baird et al. (2016) A biophysical representation of seagrass growth for application in a complex shallow-water biogeochemical model. Ecol. Modell. 325, 13–27.
   ! http://dx.doi.org/10.1016/j.ecolmodel.2015.12.011
   !
   ! Original authors: Matthew Adams, Jorn Bruggeman

   use fabm_types

   implicit none

   private

   type,extends(type_base_model),public :: type_csiro_seagrass
      type (type_bottom_state_variable_id) :: id_SG_A
      type (type_bottom_state_variable_id) :: id_SG_B

      type (type_state_variable_id)        :: id_DIC
      type (type_state_variable_id)        :: id_O2
      type (type_state_variable_id)        :: id_D_C_A, id_D_N_A, id_D_P_A
      type (type_bottom_state_variable_id) :: id_D_C_B, id_D_N_B, id_D_P_B
      type (type_bottom_state_variable_id) :: id_N
      type (type_bottom_state_variable_id) :: id_P
      type (type_dependency_id)            :: id_E_par
      type (type_dependency_id)            :: id_T
      
      real(rk) :: mu_max_SG
      real(rk) :: Omega_SG
      real(rk) :: A_par
      real(rk) :: f_below
      real(rk) :: tau_tran
      real(rk) :: K_SG_P
      real(rk) :: K_SG_N
      real(rk) :: E_comp
      real(rk) :: zeta_SG_A
      real(rk) :: zeta_SG_B
      real(rk) :: f_seed
      real(rk) :: z_root
      real(rk) :: sin_beta_blade
      real(rk) :: Q_10
   contains
      procedure :: initialize
      procedure :: do_bottom
   end type

contains

   subroutine initialize(self,configunit)
      class (type_csiro_seagrass),intent(inout),target :: self
      integer,                    intent(in)           :: configunit

      call self%get_parameter(self%mu_max_SG,      'mu_max_SG',     'd-1',                'maximum growth rate of above-ground seagrass',default=0.4_rk)
      call self%get_parameter(self%Omega_SG,       'Omega_SG',      '(g N m-2)-1',        'nitrogen-specific area of seagrass',          default=1.5_rk)
      call self%get_parameter(self%A_par,          'A_par',         '-',                  'leaf absorbance',                             default=0.7_rk)
      call self%get_parameter(self%f_below,        'f_below',       '-',                  'fraction biomass below ground',               default=0.5_rk)
      call self%get_parameter(self%tau_tran,       'tau_tran',      'd-1',                'translocation rate',                          default=0.033_rk)
      call self%get_parameter(self%K_SG_P,         'K_SG_P',        'g P m-3',            'half-saturation P uptake',                    default=0.096_rk)
      call self%get_parameter(self%K_SG_N,         'K_SG_N',        'g N m-3',            'half-saturation N uptake',                    default=0.420_rk)
      call self%get_parameter(self%E_comp,         'E_comp',        'mol photons m-2 d-1','compensation scalar PAR irradiance',          default=4.5_rk)
      call self%get_parameter(self%zeta_SG_A,      'zeta_SG_A',     'd-1',                'leaf loss rate',                              default=0.03_rk)
      call self%get_parameter(self%zeta_SG_B,      'zeta_SG_B',     'd-1',                'root loss rate',                              default=0.004_rk)
      call self%get_parameter(self%f_seed,         'f_seed',        '-',                  'seed biomass as a fraction of 63 % cover',    default=0.01_rk)
      call self%get_parameter(self%z_root,         'z_root',        'm',                  'seagrass root depth',                         default=0.15_rk)
      call self%get_parameter(self%sin_beta_blade, 'sin_beta_blade','-',                  'sine of the nadir blade angle',               default=0.5_rk)
      call self%get_parameter(self%Q_10,           'Q_10',          '-',                  'Q_10 factor',                                 default=2.0_rk)

      call self%register_state_variable(self%id_SG_A,'SG_A','g N m-2','above-ground seagrass biomass',initial_value=self%f_seed/self%Omega_SG*(1.0_rk-self%f_below))
      call self%register_state_variable(self%id_SG_B,'SG_B','g N m-2','below-ground seagrass biomass',initial_value=self%f_seed/self%Omega_SG*self%f_below)

      call self%add_to_aggregate_variable(standard_variables%total_nitrogen,  self%id_SG_A, scale_factor=1000._rk/14._rk)
      call self%add_to_aggregate_variable(standard_variables%total_carbon,    self%id_SG_A, scale_factor=1000._rk/14._rk/30._rk*550._rk)
      call self%add_to_aggregate_variable(standard_variables%total_phosphorus,self%id_SG_A, scale_factor=1000._rk/14._rk/30._rk*  1._rk)
      call self%add_to_aggregate_variable(standard_variables%total_nitrogen,  self%id_SG_B, scale_factor=1000._rk/14._rk)
      call self%add_to_aggregate_variable(standard_variables%total_carbon,    self%id_SG_B, scale_factor=1000._rk/14._rk/30._rk*550._rk)
      call self%add_to_aggregate_variable(standard_variables%total_phosphorus,self%id_SG_B, scale_factor=1000._rk/14._rk/30._rk*  1._rk)

      call self%register_state_dependency(self%id_DIC,'DIC','mmol C m-3','dissolved inorganic carbon')
      call self%register_state_dependency(self%id_O2, 'O2', 'mmol O2 m-3','oxygen')
      call self%register_state_dependency(self%id_N,  'N',  'mmol N m-3','porewater DIN concentration')
      call self%register_state_dependency(self%id_P,  'P',  'mmol P m-3','porewater DIP concentration')

      call self%register_state_dependency(self%id_D_C_A,'D_C_A','mg C m-3','sink for carbon in leaf detritus')
      call self%register_state_dependency(self%id_D_N_A,'D_N_A','mmol N m-3','sink for nitrogen in leaf detritus')
      call self%register_state_dependency(self%id_D_P_A,'D_P_A','mmol P m-3','sink for phosphorus in leaf detritus')
      call self%register_state_dependency(self%id_D_C_B,'D_C_B','mg C m-3','sink for carbon in root detritus')
      call self%register_state_dependency(self%id_D_N_B,'D_N_B','mmol N m-3','sink for nitrogen in root detritus')
      call self%register_state_dependency(self%id_D_P_B,'D_P_B','mmol P m-3','sink for phosphorus in root detritus')

      call self%register_dependency(self%id_E_par,standard_variables%downwelling_photosynthetic_radiative_flux)
      call self%register_dependency(self%id_T,    standard_variables%temperature)
   end subroutine initialize

   subroutine do_bottom(self,_ARGUMENTS_DO_BOTTOM_)
      class (type_csiro_seagrass),intent(in) :: self
      _DECLARE_ARGUMENTS_DO_BOTTOM_

      real(rk) :: SG_A,SG_B,N,P,T,E_par
      real(rk) :: mu_SG_A,k_I,k_resp,Tau,dD_N_A,dD_N_B,T_fac

      real(rk), parameter :: T_ref = 20._rk
      
      _HORIZONTAL_LOOP_BEGIN_

      _GET_HORIZONTAL_(self%id_SG_A,SG_A)
      _GET_HORIZONTAL_(self%id_SG_B,SG_B)
      _GET_HORIZONTAL_(self%id_N,N)
      _GET_HORIZONTAL_(self%id_P,P)
      _GET_(self%id_E_par,E_par)
      _GET_(self%id_T,T)

      ! E_par: from W m-2 to mol m-2 d-1
      E_par = E_par/(563900._rk)*86400._rk

      T_fac = self%Q_10**((T-T_ref)/10._rk) !TODO

      k_resp = 2*(self%E_comp*self%A_par*self%Omega_SG*self%sin_beta_blade - 5500._rk/30._rk/14._rk*self%zeta_SG_A)*SG_A
      k_I = E_par*(1.0_rk-exp(-self%A_par*self%Omega_SG*self%sin_beta_blade*SG_A))
      mu_SG_A = min(self%mu_max_SG*N/(self%K_SG_N+N), self%mu_max_SG*P/(self%K_SG_P+P), 30._rk/5500._rk*14._rk*max(0.0_rk,k_I-k_resp)/SG_A)
      Tau = (self%f_below - SG_B/(SG_A+SG_B))*(SG_A+SG_B)*self%tau_tran

      _SET_BOTTOM_ODE_(self%id_N,                -1._rk/14._rk*mu_SG_A*SG_A/86400._rk*1000._rk)  ! in mmol m-3 s-1
      _SET_BOTTOM_ODE_(self%id_P,         -1._rk/30._rk/14._rk*mu_SG_A*SG_A/86400._rk*1000._rk)  ! in mmol m-3 s-1
      _SET_BOTTOM_EXCHANGE_(self%id_DIC,-550._rk/30._rk/14._rk*mu_SG_A*SG_A/86400._rk*1000._rk)  ! in mmol m-3 s-1
      _SET_BOTTOM_EXCHANGE_(self%id_O2,  716._rk/30._rk/14._rk*mu_SG_A*SG_A/86400._rk*1000._rk)  ! in mmol O2 m-3 s-1
      _SET_BOTTOM_ODE_(self%id_SG_A,(mu_SG_A*SG_A - self%zeta_SG_A*(SG_A-self%f_seed/self%Omega_SG*(1.0_rk-self%f_below)) - Tau)/86400._rk)
      _SET_BOTTOM_ODE_(self%id_SG_B,(             - self%zeta_SG_B*(SG_B-self%f_seed/self%Omega_SG*self%f_below         ) + Tau)/86400._rk)

      ! Detritus production (separate leaf and root)
      ! Also conevrt from d-1 to s-1, and from g N to mmol N
      dD_N_A = (self%zeta_SG_A*(SG_A-self%f_seed/self%Omega_SG*(1.0_rk-self%f_below)))/86400._rk*1000._rk/14._rk
      dD_N_B = (self%zeta_SG_B*(SG_B-self%f_seed/self%Omega_SG*self%f_below         ))/86400._rk*1000._rk/14._rk
      _SET_BOTTOM_EXCHANGE_(self%id_D_N_A,                      dD_N_A)
      _SET_BOTTOM_EXCHANGE_(self%id_D_P_A,  1._rk/30._rk       *dD_N_A)
      _SET_BOTTOM_EXCHANGE_(self%id_D_C_A,550._rk/30._rk*12._rk*dD_N_A)
      _SET_BOTTOM_ODE_(self%id_D_N_B,                      dD_N_B)
      _SET_BOTTOM_ODE_(self%id_D_P_B,  1._rk/30._rk       *dD_N_B)
      _SET_BOTTOM_ODE_(self%id_D_C_B,550._rk/30._rk*12._rk*dD_N_B)

      _HORIZONTAL_LOOP_END_
   end subroutine do_bottom

end module
