#include "fabm_driver.h"

module jzwart_doc_decay

   use fabm_types

   implicit none

   private

   type, extends(type_base_model), public :: type_jzwart_doc_decay
      type (type_dependency_id) :: id_wtr
      type (type_state_variable_id) :: id_doc
      type (type_diagnostic_variable_id) :: id_decayrate
      real(rk) :: d
      
   contains
      procedure :: initialize
      procedure :: do
      
   end type

contains

   subroutine initialize(self,configunit)
      class (type_jzwart_doc_decay), intent(inout), target :: self
      integer,                          intent(in)            :: configunit
 
      call self%reigster_dependency(self%id_wtr,'wtr','degrees_c','water_temperature')
      call self%register_diagnostic_variable(self%id_decayrate,'decayrate','day-1','doc_decay_rate')
      call self%register_state_variable(self%id_doc,'doc','mol_c','dissolved_organic_carbon')
      call self%add_to_aggregate_variable(standard_variable%total_carbon,self%id_doc)
      call self%get_parameter(self%d,'d','day-1','decay_rate_doc',0.002)
      ! Register model parameters and variables here.
   end subroutine initialize

   subroutine do(self,_ARGUMENTS_DO_)
      class (type_jzwart_doc_decay),intent(in) :: self
      _DECLARE_ARGUMENTS_DO_
      
      real(rk) :: doc
      real(rk) :: doc_decay
      real(rk) :: ddoc_dt

      _LOOP_BEGIN_
      
      ! obtain concentration of DOC 
      _GET_(id_doc,doc)
      
      ! compute doc decay based on decay at 20 and water temperature 
      doc_decay = d*1.047^(id_wtr-20)
      ddoc_dt = self%doc_decay*doc
      
      ! send rates of change to FABM
      _SET_ODE_(self%id_doc,-ddoc_dt)
      
      ! send the value of diagnostic variables to FABM
      _SET_DIAGNOSTIC_(self%id_decayrate,doc_decay)
      
      _LOOP_END_
   end subroutine do

end module
