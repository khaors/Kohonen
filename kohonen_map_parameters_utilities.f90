module kohonen_map_parameters_utilities

implicit none

private

type kohonen_map_parameters
  integer :: train_option 
  integer :: number_nodes_nx,number_nodes_ny,number_nodes_nz,number_patterns
  integer :: number_variables1,number_variables2  
  integer :: number_epochs,debug_level,random_seed_ !number_clusters,
  real(kind=8) :: learning_rate
  character(len=40) :: node_type !rectangular, (to be implemented circular,hexagonal)
  character(len=40) :: debug_file,pattern_file,output_file
  character(len=40) :: distance_type !euclidean, manhattan, correlation, correlation2
  character(len=40) :: neighborhood_type !gaussian,bubble
!
  contains
    procedure,public :: print => print_parameters
end type kohonen_map_parameters

public :: kohonen_map_parameters
!
 contains
! 
 subroutine print_parameters(parameters,unit_)
   class(kohonen_map_parameters) :: parameters
   integer,intent(inout),optional :: unit_
!
   integer :: unit1
!
   if(.not. present(unit_)) then 
      unit1=6;
   else
      unit1=unit_:
   endif
   write(unit1,*) 'Kohonen Map Parameters'
   write(unit1,'(I40,A)') parameters%train_option,'!Train option' 
   write(unit1,'(A40,A)') trim(parameters%pattern_file),'!Pattern file';
   write(unit1,'(I40,A)') parameters%number_patterns,'!Number Patterns';
   write(unit1,'(I40,A)') parameters%number_variables1,'!Number Variables1';
   write(unit1,'(I40,A)') parameters%number_variables2,'!Number Variables2';
   write(unit1,'(A40,A)') trim(parameters%output_file),'!Output file';
   write(unit1,'(A40,A)') trim(parameters%debug_file),'!Debug file';
   write(unit1,'(I40,A)') parameters%debug_level,'!Debug level'; 
   write(unit1,'(I40,A)') parameters%number_nodes_nx,'!Number nodes x';
   write(unit1,'(I40,A)') parameters%number_nodes_ny,'!Number nodes y';
   write(unit1,'(I40,A)') parameters%number_nodes_nz,'!Number nodes z';
   write(unit1,'(I40,A)') parameters%number_epochs,'!Number epochs';
   write(unit1,'(f40.5,A)') parameters%learning_rate,'!Learning rate';
   write(unit1,'(I40,A)') parameters%random_seed_,'!Random seed';
!
   return
! 5  format('(A,A)')
! 6  format('(I40,1X,A20)')
! 7  format('(f40.5,A20)')
!
 end subroutine print_parameters
end module kohonen_map_parameters_utilities