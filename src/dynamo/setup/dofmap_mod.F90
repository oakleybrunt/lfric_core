!-------------------------------------------------------------------------------
! (c) The copyright relating to this work is owned jointly by the Crown, 
! Met Office and NERC 2014. 
! However, it has been created with the help of the GungHo Consortium, 
! whose members are identified at https://puma.nerc.ac.uk/trac/GungHo/wiki
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Computes the dofmaps for the 4 element spaces given grid connectivity information
! requires: list of cell next to current cell
!           list of vertices on this cell
!-------------------------------------------------------------------------------
module dofmap_mod

use num_dof_mod
use reference_element_mod
use mesh_generator_mod, only: nedge_h_g, nvert_h_g, face_on_cell, edge_on_cell, vert_on_cell

implicit none

!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------
contains 

subroutine get_dofmap(nlayers, function_space, ndf_entity)
!-----------------------------------------------------------------------------
! Subroutine to read dofmaps (in fact compute them using dofmap_populate)
!-----------------------------------------------------------------------------
  use function_space_mod, only: function_space_type
  
  implicit none
  
  integer, intent(in)       :: nlayers
  integer, intent(in)       :: ndf_entity(0:3)
  type(function_space_type) :: function_space  
  integer, allocatable      :: dofmap(:,:)
  integer :: cell
  
  integer :: ncell, ndf
  
  ncell = function_space%get_ncell()
  ndf   = function_space%get_ndf()
 
  allocate( dofmap(0:ncell,ndf) )
  
  call dofmap_populate(ncell, nlayers,ndf, ndf_entity, dofmap)

  do cell = 0, ncell
     call function_space%populate_cell_dofmap(cell,dofmap(cell,:))
  end do

end subroutine get_dofmap

subroutine dofmap_populate(ncells,nlayers,ndof_sum,ndof_entity,dofmap)
!-----------------------------------------------------------------------------
! Subroutine to populate dofmaps
!-----------------------------------------------------------------------------

  integer, intent(in) :: ncells, nlayers

! number of dofs per entity for this space
  integer, intent(in) :: ndof_entity(0:3)
! total number of dofs associated with each cell  
  integer, intent(in) :: ndof_sum

! output dofmap for this space
  integer, intent(out) :: dofmap(0:ncells,ndof_sum)

! loop counters
  integer :: i, j, k

  integer :: id, id0, jd, jdp, dof_idx
! Number of entities for a single layer  
  integer :: nvert_layer, nedge_layer, nface_layer

! entity dofmaps
  integer, allocatable :: dofmap_d0(:,:), dofmap_d1(:,:), dofmap_d2(:,:), dofmap_d3(:,:)

! dofmaps for a 3D horizontal layer
  nvert_layer = 2*nvert_h_g 
  nedge_layer = 2*nedge_h_g + nvert_h_g
  nface_layer = nedge_h_g + 2*ncells
  
  if ( ndof_entity(0) > 0 ) then
    allocate( dofmap_d0(nvert_layer,ndof_entity(0)) )
  else
    allocate( dofmap_d0(nvert_layer,1) )
  end if
  if ( ndof_entity(1) > 0 ) then  
    allocate( dofmap_d1(nedge_layer,ndof_entity(1)) )
  else
    allocate( dofmap_d1(nedge_layer,1) )
  end if  
  if ( ndof_entity(2) > 0 ) then  
    allocate( dofmap_d2(nface_layer,ndof_entity(2)) )
  else
    allocate( dofmap_d2(nface_layer,1) )
  end if
  if ( ndof_entity(3) > 0 ) then    
    allocate( dofmap_d3(ncells,     ndof_entity(3)) )
  else
    allocate( dofmap_d3(ncells,1) )
  end if 

! initialise entity dofmaps
  dofmap_d0(:,:) = 0
  dofmap_d1(:,:) = 0
  dofmap_d2(:,:) = 0
  dofmap_d3(:,:) = 0

! assume we have all possible global connectivity information
! in practice this requires connectivity
! (3,2) -> faces on cells
! (3,1) -> edges on cells
! (3,0) -> vertices on cells

  id = 1
! loop over 3 entities (cells)
  do i=1,ncells
    dof_idx = 1
! assign dofs for connectivity (3,3) (dofs in cell)
    do j=1,ndof_entity(3)
      dofmap_d3(i,j) = id
      dofmap(i,dof_idx) = dofmap_d3(i,j)
      id = id + nlayers
      dof_idx = dof_idx + 1
    end do
  
! assign dofs for connectivity (3,2) (dofs on faces)
    do j=1,nfaces_h
      jd = face_on_cell(i,j) 
      if ( dofmap_d2(jd,1) == 0 ) then
        do k=1,ndof_entity(2)
          dofmap_d2(jd,k) = id        
          id = id + nlayers
        end do
      end if
      do k=1,ndof_entity(2)
        dofmap(i,dof_idx) = dofmap_d2(jd,k)
        dof_idx = dof_idx + 1
      end do
    end do
    id0 = id
    do j=nfaces_h+1,nfaces
      jd = face_on_cell(i,j) 
      if ( dofmap_d2(jd,1) == 0 ) then
        do k=1,ndof_entity(2)
          dofmap_d2(jd,k) = id        
          id = id + nlayers + 1
        end do
      end if
      do k=1,ndof_entity(2)
        dofmap(i,dof_idx) = dofmap_d2(jd,k)
        dof_idx = dof_idx + 1
      end do
      if (j==nfaces_h+1) then
        id = id0 + 1
      else
        id = id - 1
      end if
    end do
! assign dofs for connectivity (3,1) (dofs on edges)  
    do j=1,nedges_h
      jd  = edge_on_cell(i,j)   
      jdp = edge_on_cell(i,j+nedges-nedges_h)  
      if ( dofmap_d1(jd,1) == 0 ) then
        do k=1,ndof_entity(1)
          dofmap_d1(jd,k)  = id
          dofmap_d1(jdp,k) = id+1
          id = id + nlayers + 1
        end do
      end if
    end do
    do j=5,8
      jd  = edge_on_cell(i,j) 
      if ( dofmap_d1(jd,1) == 0 ) then
        do k=1,ndof_entity(1)
          dofmap_d1(jd,k)  = id
          id = id + nlayers 
        end do
      end if
    end do
    do j=1,nedges
      jd  = edge_on_cell(i,j) 
      do k=1,ndof_entity(1)
        dofmap(i,dof_idx) = dofmap_d1(jd,k)
        dof_idx = dof_idx + 1  
      end do
    end do 
! assign dofs for connectivity (3,0) (dofs on verts)    
    do j=1,nverts_h
      jd  = vert_on_cell(i,j)
      jdp = vert_on_cell(i,j+nverts_h)
      if ( dofmap_d0(jd,1) == 0 ) then
        do k=1,ndof_entity(0)
          dofmap_d0(jd, k)  = id
          dofmap_d0(jdp,k)  = id + 1
          id = id + nlayers + 1  
        end do
      end if
    end do
    do j=1,nverts
      jd  = vert_on_cell(i,j)
      do k=1,ndof_entity(0)
        dofmap(i,dof_idx) = dofmap_d0(jd,k) 
        dof_idx = dof_idx + 1  
      end do
    end do
  end do
  
  dofmap(0,:) = 0

  if (allocated(dofmap_d0) ) deallocate( dofmap_d0 )
  if (allocated(dofmap_d1) ) deallocate( dofmap_d1 )
  if (allocated(dofmap_d2) ) deallocate( dofmap_d2 )
  if (allocated(dofmap_d3) ) deallocate( dofmap_d3 )

end subroutine dofmap_populate

end module dofmap_mod
