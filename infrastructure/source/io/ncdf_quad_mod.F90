!-----------------------------------------------------------------------------
! Copyright (c) 2017,  Met Office, on behalf of HMSO and Queen's Printer
! For further details please refer to the file LICENCE.original which you
! should have received as part of this distribution.
!-----------------------------------------------------------------------------
!>  @brief   File handler for NetCDF ugrid files.
!>  @details Implementation of the ugrid file class for quads in NetCDF format.
!-------------------------------------------------------------------------------
module ncdf_quad_mod

use constants_mod,  only : r_def, i_def, l_def, str_def, str_long
use ugrid_file_mod, only : ugrid_file_type
use netcdf,         only : nf90_max_name, nf90_open, nf90_write, nf90_noerr,   &
                           nf90_strerror, nf90_put_var, nf90_get_var,          &
                           nf90_get_att, nf90_inquire, nf90_inquire_variable,  &
                           nf90_def_var, nf90_inq_varid, nf90_int, nf90_double,&
                           nf90_clobber, nf90_enddef, nf90_inquire_dimension,  &
                           nf90_inq_dimid, nf90_def_dim, nf90_create,          &
                           nf90_inq_attname, nf90_redef,                       &
                           nf90_close, nf90_put_att, nf90_64bit_offset
use log_mod,        only : log_event, log_scratch_space, LOG_LEVEL_ERROR,      &
                           LOG_LEVEL_INFO

implicit none

private

!-------------------------------------------------------------------------------
! Module parameters
!-------------------------------------------------------------------------------

integer(i_def), parameter :: TWO  = 2                  !< Two
integer(i_def), parameter :: FOUR = 4                  !< Four

! Ranks for each variable.
integer(i_def), parameter :: MESH_RANK            = 0
integer(i_def), parameter :: MESH_FACE_NODES_RANK = 2  !< Rank of face-node connectivity arrays
integer(i_def), parameter :: MESH_EDGE_NODES_RANK = 2  !< Rank of edge-node connectivity arrays
integer(i_def), parameter :: MESH_FACE_EDGES_RANK = 2  !< Rank of face-edge connectivity arrays
integer(i_def), parameter :: MESH_FACE_LINKS_RANK = 2  !< Rank of face-face connectivity arrays
integer(i_def), parameter :: MESH_NODE_X_RANK     = 1  !< Rank of node longitude coordinate array
integer(i_def), parameter :: MESH_NODE_Y_RANK     = 1  !< Rank of node latitude  coordinate array


!-------------------------------------------------------------------------------
!> @brief   NetCDF quad file type
!> @details Implements the ugrid file type for NetCDF files storing 2D quads.
!-------------------------------------------------------------------------------

type, public, extends(ugrid_file_type) :: ncdf_quad_type

  private

  integer(i_def)           :: ncid        !< NetCDF file ID
  character(nf90_max_name) :: file_name   !< Filename
  character(nf90_max_name) :: mesh_name
  character(str_def)       :: mesh_class  !< Primitive class of mesh,
                                          !< i.e. sphere, plane

  character(str_long) :: generator_inputs !< Inputs to generator for this mesh

  ! Dimension values
  integer(i_def) :: nmesh_nodes          !< Number of nodes
  integer(i_def) :: nmesh_edges          !< Number of edges
  integer(i_def) :: nmesh_faces          !< Number of faces

  ! Dimension ids
  integer(i_def) :: nmesh_nodes_dim_id   !< NetCDF-assigned ID for number of nodes
  integer(i_def) :: nmesh_edges_dim_id   !< NetCDF-assigned ID for number of edges
  integer(i_def) :: nmesh_faces_dim_id   !< NetCDF-assigned ID for number of faces
  integer(i_def) :: two_dim_id           !< NetCDF-assigned ID for constant two
  integer(i_def) :: four_dim_id          !< NetCDF-assigned ID for constant four

  ! Variable ids
  integer(i_def) :: mesh_id              !< NetCDF-assigned ID this mesh
  integer(i_def) :: mesh_edge_nodes_id   !< NetCDF-assigned ID for the edge-node connectivity
  integer(i_def) :: mesh_face_nodes_id   !< NetCDF-assigned ID for the face-node connectivity
  integer(i_def) :: mesh_face_edges_id   !< NetCDF-assigned ID for the face-edge connectivity
  integer(i_def) :: mesh_face_links_id   !< NetCDF-assigned ID for the face-face connectivity
  integer(i_def) :: mesh_node_x_id       !< NetCDF-assigned ID for node longitude coordinates
  integer(i_def) :: mesh_node_y_id       !< NetCDF-assigned ID for node latitude  coordinates

contains

  procedure :: read_mesh
  procedure :: write_mesh
  procedure :: append_mesh
  procedure :: get_dimensions
  procedure :: get_mesh_names
  procedure :: get_nmeshes
  procedure :: is_mesh_present
  procedure :: file_open
  procedure :: file_close
  procedure :: file_new

end type ncdf_quad_type

!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------
contains

!-------------------------------------------------------------------------------
!>  @brief Open an existing NetCDF file.
!>
!>  @param[in,out]  self      The NetCDF file object.
!>  @param[in]      file_name Name of the file to open.
!-------------------------------------------------------------------------------

subroutine file_open(self, file_name)
  implicit none

  ! Arguments
  class(ncdf_quad_type), intent(inout) :: self
  character(len=*),      intent(in)    :: file_name

  ! Internal variables
  integer(i_def) :: ierr

  self%file_name = file_name

  ierr = nf90_open( trim(self%file_name), nf90_write, self%ncid )
  if (ierr /= nf90_noerr) then
    write (log_scratch_space,*) 'Error in ncdf_open: '                         &
                    // trim(nf90_strerror(ierr)) //' : '// trim(self%file_name)
    call log_event( trim(log_scratch_space), LOG_LEVEL_ERROR )
  end if

  return
end subroutine file_open


!-------------------------------------------------------------------------------
!>  @brief Closes a NetCDF file.
!>
!>  @param[in]  self   The NetCDF file object.
!-------------------------------------------------------------------------------

subroutine file_close(self)
  implicit none

  ! Arguments
  class(ncdf_quad_type), intent(inout) :: self

  ! Internal variables
  integer(i_def) :: ierr

  ierr = nf90_close( self%ncid )
  if (ierr /= nf90_noerr) then
    write (log_scratch_space,*) 'Error in ncdf_close: '                        &
                    // trim(nf90_strerror(ierr)) //' : '// trim(self%file_name)
    call log_event( trim(log_scratch_space), LOG_LEVEL_ERROR )
  end if

  return
end subroutine file_close

!-------------------------------------------------------------------------------
!>  @brief       Create a new NetCDF file.
!>  @description Creates an opens a new, clean NetCDF file. If a file of the
!>               same name already exists, this routine will clobber it.
!>
!>  @param[in,out]  self      The NetCDF file object.
!>  @param[in]      file_name The name of the file to create/open.
!-------------------------------------------------------------------------------

subroutine file_new(self, file_name)
  implicit none

  ! Arguments
  class(ncdf_quad_type), intent(inout) :: self
  character(len=*),      intent(in)    :: file_name

  ! Internal variables
  integer(i_def) :: ierr

  self%file_name = file_name

  ! Create the NetCDF file with 64-bit offsets to support large file sizes
  ierr = nf90_create( path=trim(self%file_name), &
                      cmode=ior(nf90_clobber,nf90_64bit_offset), &
                      ncid=self%ncid )

  if (ierr /= NF90_NOERR) then
    write (log_scratch_space,*) 'Error in ncdf_create: '                       &
                    // trim(nf90_strerror(ierr)) //' : '// trim(self%file_name)
    call log_event( trim(log_scratch_space), LOG_LEVEL_ERROR )
  end if

  return
end subroutine file_new

!-------------------------------------------------------------------------------
!>  @brief   Defines NetCDF dimensions in the NetCDF file.
!>  @details Sets dimension lengths in the NetCDF file, and sets the associated
!>           dimension ids in the NetCDF file object. The dimension lengths are
!>           used for sizes of other arrays within the NetCDF file.
!>
!>  @param[in,out]  self   The NetCDF file object.
!-------------------------------------------------------------------------------

subroutine define_dimensions(self)
  implicit none

  ! Arguments
  class(ncdf_quad_type), intent(inout) :: self

  ! Internal variables
  integer(i_def) :: ierr

  character(str_long) :: routine
  character(str_long) :: cmess
  character(str_long) :: dim_name

  routine = 'define_dimensions'
  cmess = ''

  ! Define dimensions connected to the mesh
  dim_name = 'n'//trim(self%mesh_name)//'_node'
  cmess = 'Defining '//trim(dim_name)
  ierr = nf90_def_dim( self%ncid, trim(dim_name), &
                       self%nmesh_nodes, self%nmesh_nodes_dim_id )
  call check_err(ierr, routine, cmess)

  dim_name = 'n'//trim(self%mesh_name)//'_edge'
  cmess = 'Defining '//trim(dim_name)
  ierr = nf90_def_dim( self%ncid, trim(dim_name), &
                       self%nmesh_edges,self%nmesh_edges_dim_id )
  call check_err(ierr, routine, cmess)

  dim_name = 'n'//trim(self%mesh_name)//'_face'
  cmess = 'Defining '//trim(dim_name)
  ierr = nf90_def_dim( self%ncid, trim(dim_name), &
                       self%nmesh_faces, self%nmesh_faces_dim_id )
  call check_err(ierr, routine, cmess)

  ! If the file is being appended to, constants may already exist
  ! in the NetCDF file. Trying to redefine the same variable name
  ! will throw a error, so check to see if constants are present
  ! already.
  ierr = nf90_inq_dimid (self%ncid, 'Two', self%two_dim_id)
  if (ierr /= nf90_noerr) then
    ierr = nf90_def_dim(self%ncid, 'Two', TWO, self%two_dim_id)
    cmess = 'Defining Two'
    call check_err(ierr, routine, cmess)
  end if

  ierr = nf90_inq_dimid(self%ncid, 'Four', self%four_dim_id)
  if (ierr /= nf90_noerr) then
    cmess = 'Defining Four'
    ierr = nf90_def_dim(self%ncid, 'Four', FOUR, self%four_dim_id)
    call check_err(ierr,routine, cmess)
  end if

  return
end subroutine define_dimensions


!-------------------------------------------------------------------------------
!>  @brief   Defines NetCDF variables in the netCDF file.
!>  @details Tells NetCDF what variables are going to be in the file.
!>           Array lengths are specified via the pre-existing NetCDF dimension
!>           IDs, which were obtained elsewhere in this module.
!>
!>  @param[in,out]  self   The NetCDF file object.
!-------------------------------------------------------------------------------

subroutine define_variables(self)

  implicit none

  ! Arguments
  class(ncdf_quad_type), intent(inout) :: self

  ! Internal variables
  integer(i_def) :: ierr
  integer(i_def) :: zero_sized(0)

  ! Variable shapes
  integer(i_def) :: mesh_face_nodes_dims(MESH_FACE_NODES_RANK)
  integer(i_def) :: mesh_edge_nodes_dims(MESH_EDGE_NODES_RANK)
  integer(i_def) :: mesh_face_edges_dims(MESH_FACE_EDGES_RANK)
  integer(i_def) :: mesh_face_links_dims(MESH_FACE_LINKS_RANK)
  integer(i_def) :: mesh_node_x_dims(MESH_NODE_X_RANK)
  integer(i_def) :: mesh_node_y_dims(MESH_NODE_Y_RANK)

  character(str_long) :: routine
  character(str_long) :: cmess
  character(str_long) :: var_name

  routine = 'define_variables'
  cmess = ''


  cmess = 'Defining '//trim(self%mesh_name)
  ierr = nf90_def_var( self%ncid, trim(self%mesh_name), &
                      nf90_int, zero_sized, self%mesh_id )
  call check_err(ierr, routine, cmess)


  mesh_face_nodes_dims(1) = self%four_dim_id
  mesh_face_nodes_dims(2) = self%nmesh_faces_dim_id
  var_name = trim(self%mesh_name)//'_face_nodes'
  cmess = 'Defining '//trim(var_name)
  ierr = nf90_def_var( self%ncid, trim(var_name),       &
                       nf90_int, mesh_face_nodes_dims,  &
                       self%mesh_face_nodes_id )
  call check_err(ierr, routine, cmess)

  mesh_edge_nodes_dims(1) = self%two_dim_id
  mesh_edge_nodes_dims(2) =self%nmesh_edges_dim_id
  var_name = trim(self%mesh_name)//'_edge_nodes'
  cmess = 'Defining '//trim(var_name)
  ierr = nf90_def_var( self%ncid, trim(var_name),       &
                       nf90_int, mesh_edge_nodes_dims,  &
                       self%mesh_edge_nodes_id )
  call check_err(ierr, routine, cmess)

  mesh_face_edges_dims(1) = self%four_dim_id
  mesh_face_edges_dims(2) = self%nmesh_faces_dim_id
  var_name = trim(self%mesh_name)//'_face_edges'
  cmess = 'Defining '//trim(var_name)
  ierr = nf90_def_var( self%ncid, trim(var_name),       &
                       nf90_int, mesh_face_edges_dims,  &
                       self%mesh_face_edges_id )
  call check_err(ierr, routine, cmess)

  mesh_face_links_dims(1) = self%four_dim_id
  mesh_face_links_dims(2) = self%nmesh_faces_dim_id
  var_name = trim(self%mesh_name)//'_face_links'
  cmess = 'Defining '//trim(var_name)
  ierr = nf90_def_var( self%ncid, trim(var_name),       &
                       nf90_int, mesh_face_links_dims,  &
                       self%mesh_face_links_id )
  call check_err(ierr, routine, cmess)

  mesh_node_x_dims(1) = self%nmesh_nodes_dim_id
  var_name = trim(self%mesh_name)//'_node_x'
  cmess = 'Defining '//trim(var_name)
  ierr = nf90_def_var( self%ncid, trim(var_name),       &
                       nf90_double, mesh_node_x_dims,   &
                       self%mesh_node_x_id )
  call check_err(ierr, routine, cmess)

  mesh_node_y_dims(1) = self%nmesh_nodes_dim_id
  var_name = trim(self%mesh_name)//'_node_y'
  cmess = 'Defining '//trim(var_name)
  ierr = nf90_def_var( self%ncid, trim(var_name),       &
                       nf90_double, mesh_node_y_dims,   &
                       self%mesh_node_y_id )
  call check_err(ierr, routine, cmess)

  return
end subroutine define_variables

!-------------------------------------------------------------------------------
!>  @brief   Assigns attributes to the NetCDF variables.
!>  @details Adds additional information to NetCDF variables that should have
!>           already been defined elsewhere in this module.  Attributes include
!>           variable names and descriptions.
!>
!>  @param[in]   self   The NetCDF file object.
!-------------------------------------------------------------------------------

subroutine assign_attributes(self)
  implicit none

  ! Arguments
  class(ncdf_quad_type), intent(in) :: self

  ! Internal variables
  integer(i_def) :: ierr, id

  character(str_def)  :: std_x_name
  character(str_def)  :: std_y_name
  character(str_def)  :: long_x_name
  character(str_def)  :: long_y_name
  character(str_def)  :: coord_units

  character(str_long) :: routine
  character(str_long) :: cmess

  character(nf90_max_name) :: var_name
  character(str_long) :: attname

  routine = 'assign_attributes'
  cmess = ''


  id = self%mesh_id
  ierr = nf90_inquire_variable( ncid=self%ncid, varid=id, name=var_name )
  call check_err(ierr, routine, cmess)

  !===================================================================
  attname = 'cf_role'
  cmess   = 'Adding attribute "'//trim(attname)// &
            '" to variable "'//trim(var_name)//'"'
  ierr = nf90_put_att( self%ncid, id, trim(attname), 'mesh_topology')
  call check_err(ierr, routine, cmess)


  attname = 'mesh_class'
  cmess   = 'Adding attribute "'//trim(attname)// &
            '" to variable "'//trim(var_name)//'"'
  ierr = nf90_put_att( self%ncid, id, trim(attname), &
                       trim(self%mesh_class))
  call check_err(ierr, routine, cmess)

  attname = 'generator_inputs'
  cmess   = 'Adding attribute "'//trim(attname)// &
            '" to variable "'//trim(var_name)//'"'
  ierr = nf90_put_att( self%ncid, id, trim(attname), &
                       trim(self%generator_inputs) )

  call check_err(ierr, routine, cmess)

  attname = 'long_name'
  cmess   = 'Adding attribute "'//trim(attname)// &
            '" to variable "'//trim(var_name)//'"'
  ierr = nf90_put_att( self%ncid, id, trim(attname), &
                       'Topology data of 2D unstructured mesh')
  call check_err(ierr, routine, cmess)

  attname = 'topology_dimension'
  cmess   = 'Adding attribute "'//trim(attname)// &
            '" to variable "'//trim(var_name)//'"'
  ierr = nf90_put_att( self%ncid, id, trim(attname), [2])
  call check_err(ierr, routine, cmess)

  attname = 'node_coordinates'
  cmess   = 'Adding attribute "'//trim(attname)// &
            '" to variable "'//trim(var_name)//'"'
  ierr = nf90_put_att( self%ncid, id, trim(attname),      &
                       trim(self%mesh_name)//'_node_x '// &
                       trim(self%mesh_name)//'_node_y' )
  call check_err(ierr, routine, cmess)

  attname = 'face_node_connectivity'
  cmess   = 'Adding attribute "'//trim(attname)// &
            '" to variable "'//trim(var_name)//'"'
  ierr = nf90_put_att( self%ncid, id, trim(attname), &
                       trim(self%mesh_name)//'_face_nodes')
  call check_err(ierr, routine, cmess)

  attname = 'edge_node_connectivity'
  cmess   = 'Adding attribute "'//trim(attname)// &
            '" to variable "'//trim(var_name)//'"'
  ierr = nf90_put_att( self%ncid, id, trim(attname), &
                       trim(self%mesh_name)//'_edge_nodes' )
  call check_err(ierr, routine, cmess)

  attname = 'face_edge_connectivity'
  cmess   = 'Adding attribute "'//trim(attname)// &
            '" to variable "'//trim(var_name)//'"'
  ierr = nf90_put_att( self%ncid, id, trim(attname), &
                       trim(self%mesh_name)//'_face_edges' )
  call check_err(ierr, routine, cmess)

  attname = 'face_face_connectivity'
  cmess   = 'Adding attribute "'//trim(attname)// &
            '" to variable "'//trim(var_name)//'"'
  ierr = nf90_put_att( self%ncid, id, trim(attname), &
                       trim(self%mesh_name)//'_face_links' )
  call check_err(ierr, routine, cmess)



  id = self%mesh_face_nodes_id
  ierr = nf90_inquire_variable( ncid=self%ncid, varid=id, name=var_name )
  call check_err(ierr, routine, cmess)
  !===================================================================

  attname = 'cf_role'
  cmess   = 'Adding attribute "'//trim(attname)// &
            '" to variable "'//trim(var_name)//'"'
  ierr    = nf90_put_att( self%ncid, id, trim(attname), &
                          'face_node_connectivity' )
  call check_err(ierr, routine, cmess)

  attname = 'long_name'
  cmess   = 'Adding attribute "'//trim(attname)// &
            '" to variable "'//trim(var_name)//'"'
  ierr = nf90_put_att( self%ncid, id, trim(attname), &
                       'Maps every quadrilateral face to its four corner nodes.')
  call check_err(ierr, routine, cmess)

  attname = 'start_index'
  cmess   = 'Adding attribute "'//trim(attname)// &
            '" to variable "'//trim(var_name)//'"'
  ierr = nf90_put_att( self%ncid, id, trim(attname), [1])
  call check_err(ierr, routine, cmess)



  id = self%mesh_edge_nodes_id
  ierr = nf90_inquire_variable( ncid=self%ncid, varid=id, name=var_name )
  call check_err(ierr, routine, cmess)
  !===================================================================
  attname = 'cf_role'
  cmess   = 'Adding attribute "'//trim(attname)// &
            '" to variable "'//trim(var_name)//'"'
  ierr = nf90_put_att( self%ncid, id, trim(attname), &
                       'edge_node_connectivity' )
  call check_err(ierr, routine, cmess)

  attname = 'long_name'
  cmess   = 'Adding attribute "'//trim(attname)// &
            '" to variable "'//trim(var_name)//'"'
  ierr = nf90_put_att( self%ncid, id, trim(attname), &
                       'Maps every edge to the two nodes that it connects.' )
  call check_err(ierr, routine, cmess)

  attname = 'start_index'
  cmess   = 'Adding attribute "'//trim(attname)// &
            '" to variable "'//trim(var_name)//'"'
  ierr = nf90_put_att( self%ncid, id, trim(attname), [1] )
  call check_err(ierr, routine, cmess)



  id = self%mesh_face_edges_id
  ierr = nf90_inquire_variable( ncid=self%ncid, varid=id, name=var_name )
  call check_err(ierr, routine, cmess)
  !===================================================================
  attname = 'cf_role'
  cmess   = 'Adding attribute "'//trim(attname)// &
            '" to variable "'//trim(var_name)//'"'
  ierr = nf90_put_att( self%ncid, id, trim(attname), &
                       'face_edge_connectivity' )
  call check_err(ierr, routine, cmess)

  attname = 'long_name'
  cmess   = 'Adding attribute "'//trim(attname)// &
            '" to variable "'//trim(var_name)//'"'
  ierr = nf90_put_att( self%ncid, id, trim(attname), &
                       'Maps every quadrilateral face to its four edges.' )
  call check_err(ierr, routine, cmess)

  attname = 'start_index'
  cmess   = 'Adding attribute "'//trim(attname)// &
            '" to variable "'//trim(var_name)//'"'
  ierr = nf90_put_att( self%ncid, id, trim(attname), [1] )
  call check_err(ierr, routine, cmess)



  id   = self%mesh_face_links_id
  ierr = nf90_inquire_variable( ncid=self%ncid, varid=id, name=var_name )
  call check_err(ierr, routine, cmess)
  !===================================================================
  attname = 'cf_role'
  cmess   = 'Adding attribute "'//trim(attname)// &
            '" to variable "'//trim(var_name)//'"'
  ierr = nf90_put_att( self%ncid, id, trim(attname), &
                       'face_face_connectivity' )
  call check_err(ierr, routine, cmess)

  attname = 'long_name'
  cmess   = 'Adding attribute "'//trim(attname)// &
            '" to variable "'//trim(var_name)//'"'
  ierr = nf90_put_att( self%ncid, id, trim(attname), &
                       'Indicates which other faces neighbour each face.' )
  call check_err(ierr, routine, cmess)

  attname = 'start_index'
  cmess   = 'Adding attribute "'//trim(attname)// &
            '" to variable "'//trim(var_name)//'"'
  ierr = nf90_put_att( self%ncid, id, trim(attname), [1] )
  call check_err(ierr, routine, cmess)

  attname = 'flag_values'
  cmess   = 'Adding attribute "'//trim(attname)// &
            '" to variable "'//trim(var_name)//'"'
  ierr = nf90_put_att( self%ncid, id, trim(attname), [-1] )
  call check_err(ierr, routine, cmess)

  attname = 'flag_meanings'
  cmess   = 'Adding attribute "'//trim(attname)// &
            '" to variable "'//trim(var_name)//'"'
  ierr = nf90_put_att( self%ncid, id, trim(attname), 'out_of_mesh' )
  call check_err(ierr, routine, cmess)



  select case (trim(self%mesh_class))
  case ('sphere')
    std_x_name  = 'longitude'
    std_y_name  = 'latitude'
    long_x_name = 'longitude of 2D mesh nodes.'
    long_y_name = 'latitude of 2D mesh nodes.'
    coord_units = 'radians'
  case ('plane')
    std_x_name  = 'projection_x_coordinate'
    std_y_name  = 'projection_y_coordinate'
    long_x_name = 'x coordinate of 2D mesh nodes.'
    long_y_name = 'y coordinate of 2D mesh nodes.'
    coord_units = 'm'
  end select



  id   = self%mesh_node_x_id
  ierr = nf90_inquire_variable( ncid=self%ncid, varid=id, name=var_name )
  call check_err(ierr, routine, cmess)
  !===================================================================
  attname = 'standard_name'
  cmess   = 'Adding attribute "'//trim(attname)// &
            '" to variable "'//trim(var_name)//'"'
  ierr = nf90_put_att( self%ncid, id, trim(attname), std_x_name )
  call check_err(ierr, routine, cmess)

  attname = 'long_name'
  cmess   = 'Adding attribute "'//trim(attname)// &
            '" to variable "'//trim(var_name)//'"'
  ierr = nf90_put_att( self%ncid, id, trim(attname), long_x_name )
  call check_err(ierr, routine, cmess)

  attname = 'units'
  cmess   = 'Adding attribute "'//trim(attname)// &
            '" to variable "'//trim(var_name)//'"'
  ierr = nf90_put_att( self%ncid, id, trim(attname), coord_units )
  call check_err(ierr, routine, cmess)



  id = self%mesh_node_y_id
  ierr = nf90_inquire_variable( ncid=self%ncid, varid=id, name=var_name )
  call check_err(ierr, routine, cmess)
  !===================================================================
  attname = 'standard_name'
  cmess   = 'Adding attribute "'//trim(attname)// &
            '" to variable "'//trim(var_name)//'"'
  ierr = nf90_put_att( self%ncid, id, trim(attname), std_y_name )
  call check_err(ierr, routine, cmess)

  attname = 'long_name'
  cmess   = 'Adding attribute "'//trim(attname)// &
            '" to variable "'//trim(var_name)//'"'
  ierr = nf90_put_att( self%ncid, id, trim(attname), long_y_name )
  call check_err(ierr, routine, cmess)

  attname = 'units'
  cmess   = 'Adding attribute "'//trim(attname)// &
            '" to variable "'//trim(var_name)//'"'
  ierr = nf90_put_att( self%ncid, id, trim(attname), coord_units )
  call check_err(ierr, routine, cmess)


  return
end subroutine assign_attributes

!-------------------------------------------------------------------------------
!>  @brief   Gets dimension ids and variable ids from the open NetCDF file.
!>  @details NetCDF files refer to dimensions and variables by an id, the value
!>           of which is determined by the NetCDF library. This routine finds
!>           dimension and variable ids for all variables of interest in the
!>           open NetCDF file.
!>
!>  @param[in,out] self      The NetCDF file object.
!>  @param[in]     mesh_name Name of mesh topology to get ids for
!-------------------------------------------------------------------------------

subroutine inquire_ids(self, mesh_name)

  implicit none

  ! Arguments
  type(ncdf_quad_type), intent(inout) :: self

  character(str_def), intent(in) :: mesh_name

  ! Internal variables
  integer(i_def) :: ierr

  character(nf90_max_name) :: dim_name
  character(nf90_max_name) :: var_name

  logical(l_def) :: mesh_present

  character(str_long) :: routine
  character(str_long) :: cmess

  routine = 'inquire_ids'
  cmess = ''

  mesh_present = self%is_mesh_present(mesh_name)

  if (.not. mesh_present) then
    write(log_scratch_space,'(A)') &
         'Mesh '//trim(mesh_name)//' not present in file'
    call log_event(trim(log_scratch_space), LOG_LEVEL_ERROR)
  end if

  ierr = nf90_inq_varid( self%ncid, trim(mesh_name), self%mesh_id )

  ! Numbers of entities
  dim_name = 'n'//trim(mesh_name)//'_node'
  cmess = 'Getting id for '//trim(dim_name)
  ierr = nf90_inq_dimid( self%ncid, trim(dim_name), &
                         self%nmesh_nodes_dim_id )
  call check_err(ierr, routine, cmess)

  dim_name = 'n'//trim(mesh_name)//'_edge'
  cmess = 'Getting id for '//trim(dim_name)
  ierr = nf90_inq_dimid( self%ncid, trim(dim_name), &
                        self%nmesh_edges_dim_id )
  call check_err(ierr, routine, cmess)

  dim_name = 'n'//trim(mesh_name)//'_face'
  cmess = 'Getting id for '//trim(dim_name)
  ierr = nf90_inq_dimid( self%ncid, trim(dim_name), &
                         self%nmesh_faces_dim_id )
  call check_err(ierr, routine, cmess)


  ! Node coordinates
  var_name = trim(mesh_name)//'_node_x'
  cmess = 'Getting id for '//trim(var_name)
  ierr = nf90_inq_varid( self%ncid, trim(var_name), &
                         self%mesh_node_x_id )
  call check_err(ierr, routine, cmess)

  var_name = trim(mesh_name)//'_node_y'
  cmess = 'Getting id for '//trim(var_name)
  ierr = nf90_inq_varid( self%ncid, trim(var_name), &
                         self%mesh_node_y_id )
  call check_err(ierr, routine, cmess)

  ! Face node connectivity
  var_name = trim(mesh_name)//'_face_nodes'
  cmess = 'Getting id for '//trim(var_name)
  ierr = nf90_inq_varid( self%ncid, trim(var_name), &
                         self%mesh_face_nodes_id )
  call check_err(ierr, routine, cmess)

  ! Edge node connectivity
  var_name = trim(mesh_name)//'_edge_nodes'
  cmess = 'Getting id for '//trim(var_name)
  ierr = nf90_inq_varid( self%ncid, trim(var_name), &
                         self%mesh_edge_nodes_id )
  call check_err(ierr, routine, cmess)

  ! Face edge connectivity
  var_name = trim(mesh_name)//'_face_edges'
  cmess = 'Getting id for '//trim(var_name)
  ierr = nf90_inq_varid( self%ncid, trim(var_name), &
                         self%mesh_face_edges_id )
  call check_err(ierr, routine, cmess)

  ! Face face connectivity
  var_name = trim(mesh_name)//'_face_links'
  cmess = 'Getting id for '//trim(var_name)
  ierr = nf90_inq_varid( self%ncid, trim(var_name), &
                         self%mesh_face_links_id )
  call check_err(ierr, routine, cmess)


  dim_name = 'Two'
  cmess = 'Getting id for '//trim(dim_name)
  ierr = nf90_inq_dimid(self%ncid, trim(dim_name),  self%two_dim_id)
  call check_err(ierr,routine, cmess)

  dim_name = 'Four'
  cmess = 'Getting id for '//trim(dim_name)
  ierr = nf90_inq_dimid(self%ncid, trim(dim_name), self%four_dim_id)
  call check_err(ierr, routine, cmess)

  return
end subroutine inquire_ids

!-------------------------------------------------------------------------------
!>  @brief   Calls logger on error.
!>  @details Checks the error code returned by the NetCDF file. If an error is
!>           detected, the relevant error message is passed to the logger.
!>
!>  @param[in] ierr    The error code to check.
!>  @param[in] routine The routine name that call the error check
!>  @param[in] cmess   Comment message for the error report
!-------------------------------------------------------------------------------

subroutine check_err(ierr, routine, cmess)
  implicit none

  ! Arguments
  integer(i_def),      intent(in) :: ierr
  character(str_long), intent(in) :: routine
  character(str_long), intent(in) :: cmess

  if (ierr /= NF90_NOERR) then
    write(log_scratch_space,*) 'Error in ncdf_quad ['//trim(routine)//']: '//  &
        trim(cmess) // ' ' // trim(nf90_strerror(ierr))
    call log_event( trim(log_scratch_space), LOG_LEVEL_ERROR )
  end if

  return
end subroutine check_err

!-------------------------------------------------------------------------------
!>  @brief   Returns the number of mesh topologies described in this file.
!>  @details Scans the variable attributes for cf_role = mesh_topology as
!>           a flag that the variable name relates to a mesh.
!>
!>  @return nmeshes    Integer, The number of mesh topologies in the NetCDf file
!-------------------------------------------------------------------------------
function get_nmeshes(self) result (nmeshes)

  implicit none

  class(ncdf_quad_type), intent(in) :: self

  integer(i_def) :: nmeshes

  character(str_def), allocatable :: mesh_names(:)
  integer(i_def), parameter :: max_n_topologies = 20


  allocate( mesh_names (max_n_topologies) )
  call scan_for_topologies(self, mesh_names, nmeshes)

  deallocate(mesh_names)

  return
end function get_nmeshes

!-------------------------------------------------------------------------------
!>  @brief Returns the names of mesh topologies described in this file.
!>
!>  @param[out] mesh_names  Character[:], Names of the mesh_topologies
!>                          in the NetCDF file
!-------------------------------------------------------------------------------
subroutine get_mesh_names(self, mesh_names)

  implicit none

  class(ncdf_quad_type), intent(in)  :: self
  character(len=*),      intent(out) :: mesh_names(:)

  integer(i_def) :: nmeshes

  call scan_for_topologies(self, mesh_names, nmeshes)

  return
end subroutine get_mesh_names

!-------------------------------------------------------------------------------
!>  @brief   Gets dimension information from the NetCDF file, as integers.
!>  @details Calls NetCDF inquiry functions to determine array lengths, such as
!>           the number of nodes.
!>
!>  @param[in,out]   self                   The NetCDF file object.
!>  @param[in]       mesh_name              Name of the mesh topology
!>  @param[out]      num_nodes              The number of nodes on the mesh.
!>  @param[out]      num_edges              The number of edges on the mesh.
!>  @param[out]      num_faces              The number of faces on the mesh.
!>  @param[out]      num_nodes_per_face     The number of nodes per face.
!>  @param[out]      num_edges_per_face     The number of edges per face.
!>  @param[out]      num_nodes_per_edge     The number of nodes per edge.
!>  @param[out]      max_num_faces_per_node The maximum number of faces surrounding a node.
!-------------------------------------------------------------------------------

subroutine get_dimensions( self,               &
                           mesh_name,          &
                           num_nodes,          &
                           num_edges,          &
                           num_faces,          &
                           num_nodes_per_face, &
                           num_edges_per_face, &
                           num_nodes_per_edge, &
                           max_num_faces_per_node )
  implicit none

  ! Arguments
  class(ncdf_quad_type),  intent(inout) :: self

  character(str_def), intent(in)  :: mesh_name
  integer(i_def),     intent(out) :: num_nodes
  integer(i_def),     intent(out) :: num_edges
  integer(i_def),     intent(out) :: num_faces
  integer(i_def),     intent(out) :: num_nodes_per_face
  integer(i_def),     intent(out) :: num_edges_per_face
  integer(i_def),     intent(out) :: num_nodes_per_edge
  integer(i_def),     intent(out) :: max_num_faces_per_node


  integer(i_def) :: ierr

  character(str_long) :: routine
  character(str_long) :: cmess

  call inquire_ids(self, mesh_name)

  routine='get_dimensions'
  cmess=''

  ! Get dimension lengths
  ierr = nf90_inquire_dimension( self%ncid,               &
                                 self%nmesh_nodes_dim_id, &
                                 len=self%nmesh_nodes )
  call check_err(ierr, routine, cmess)

  ierr = nf90_inquire_dimension( self%ncid,               &
                                 self%nmesh_edges_dim_id, &
                                 len=self%nmesh_edges )
  call check_err(ierr, routine, cmess)

  ierr = nf90_inquire_dimension( self%ncid,               &
                                 self%nmesh_faces_dim_id, &
                                 len=self%nmesh_faces )
  call check_err(ierr, routine, cmess)

  num_nodes = self%nmesh_nodes
  num_edges = self%nmesh_edges
  num_faces = self%nmesh_faces

  num_nodes_per_face = 4
  num_edges_per_face = 4
  num_nodes_per_edge = 2
  max_num_faces_per_node = 4

  return
end subroutine get_dimensions

!-------------------------------------------------------------------------------
!>  @brief   Read data from the NetCDF file.
!>  @details Reads coordinate and connectivity information from the NetCDF file.
!>
!>  @param[in,out]  self                     The NetCDF file object.
!>  @param[in]      mesh_name                Name of the mesh topology
!>  @param[in]      mesh_class               Primitive class of mesh
!>  @param[in]      generator_inputs         Inputs used to generate mesh
!>  @param[out]     node_coordinates         long/lat coordinates of each node.
!>  @param[out]     face_node_connectivity   Nodes adjoining each face.
!>  @param[out]     edge_node_connectivity   Nodes adjoining each edge.
!>  @param[out]     face_edge_connectivity   Edges adjoining each face.
!>  @param[out]     face_face_connectivity   Faces adjoining each face (links).
!-------------------------------------------------------------------------------

subroutine read_mesh( self, mesh_name, mesh_class, generator_inputs,  &
                      node_coordinates,                               &
                      face_node_connectivity, edge_node_connectivity, &
                      face_edge_connectivity, face_face_connectivity )
  implicit none

  ! Arguments
  class(ncdf_quad_type),  intent(inout) :: self

  character(str_def),  intent(in)  :: mesh_name
  character(str_def),  intent(out) :: mesh_class
  character(str_long), intent(out) :: generator_inputs
  real(r_def),         intent(out) :: node_coordinates(:,:)
  integer(i_def),      intent(out) :: face_node_connectivity(:,:)
  integer(i_def),      intent(out) :: edge_node_connectivity(:,:)
  integer(i_def),      intent(out) :: face_edge_connectivity(:,:)
  integer(i_def),      intent(out) :: face_face_connectivity(:,:)

  ! Internal variables
  integer(i_def) :: ierr

  character(str_long) :: routine
  character(str_long) :: cmess

  call inquire_ids(self, mesh_name)

  routine = 'read_mesh'
  cmess   = ''


  ! Mesh class
  ierr = nf90_get_att( self%ncid, self%mesh_id, &
                       'mesh_class', mesh_class )

  ! Generator inputs
  ierr = nf90_get_att( self%ncid, self%mesh_id, &
                       'generator_inputs', generator_inputs )

  ! Node coordinates
  ierr = nf90_get_var( self%ncid, self%mesh_node_x_id, &
                       node_coordinates(1,:))
  call check_err(ierr, routine, cmess)

  ierr = nf90_get_var( self%ncid, self%mesh_node_y_id, &
                       node_coordinates(2,:))
  call check_err(ierr, routine, cmess)



  ! Face node connectivity
  ierr = nf90_get_var( self%ncid, self%mesh_face_nodes_id, &
                       face_node_connectivity(:,:) )
  call check_err(ierr, routine, cmess)



  ! Edge node connectivity
  ierr = nf90_get_var( self%ncid, self%mesh_edge_nodes_id, &
                       edge_node_connectivity(:,:) )
  call check_err(ierr, routine, cmess)



  ! Face edge connectivity
  ierr = nf90_get_var( self%ncid, self%mesh_face_edges_id, &
                       face_edge_connectivity(:,:) )
  call check_err(ierr, routine, cmess)



  ! Face face connectivity
  ierr = nf90_get_var( self%ncid, self%mesh_face_links_id, &
                       face_face_connectivity(:,:))
  call check_err(ierr, routine, cmess)

  return
end subroutine read_mesh

!-------------------------------------------------------------------------------
!>  @brief   Writes data to the NetCDF file.
!>  @details Writes dimension, coordinate and connectivity information
!>           to the NetCDF file.
!>
!>  @param[in,out]  self                     The NetCDF file object.
!>  @param[in]      mesh_name                Name of the mesh topology.
!>  @param[in]      mesh_class               Primitive class of mesh.
!>  @param[in]      generator_inputs         Inputs used to create this mesh
!>                                           from the mesh_generator
!>  @param[in]      num_nodes                The number of nodes on the mesh.
!>  @param[in]      num_edges                The number of edges on the mesh.
!>  @param[in]      num_faces                The number of faces on the mesh.
!>  @param[in]      node_coordinates         long/lat coordinates of each node.
!>  @param[in]      face_node_connectivity   Nodes adjoining each face.
!>  @param[in]      edge_node_connectivity   Nodes adjoining each edge.
!>  @param[in]      face_edge_connectivity   Edges adjoining each face.
!>  @param[in]      face_face_connectivity   Faces adjoining each face (links).
!-------------------------------------------------------------------------------

subroutine write_mesh( self, mesh_name, mesh_class, generator_inputs,     &
                       num_nodes, num_edges, num_faces, node_coordinates, &
                       face_node_connectivity, edge_node_connectivity,    &
                       face_edge_connectivity, face_face_connectivity )
  implicit none

  ! Arguments
  class(ncdf_quad_type),  intent(inout) :: self

  character(str_def),  intent(in) :: mesh_name
  character(str_def),  intent(in) :: mesh_class
  character(str_long), intent(in) :: generator_inputs
  integer(i_def),      intent(in) :: num_nodes
  integer(i_def),      intent(in) :: num_edges
  integer(i_def),      intent(in) :: num_faces
  real(r_def),         intent(in) :: node_coordinates(:,:)
  integer(i_def),      intent(in) :: face_node_connectivity(:,:)
  integer(i_def),      intent(in) :: edge_node_connectivity(:,:)
  integer(i_def),      intent(in) :: face_edge_connectivity(:,:)
  integer(i_def),      intent(in) :: face_face_connectivity(:,:)

  ! Internal variables
  integer(i_def) :: ierr
  character(str_long) :: routine
  character(str_long) :: cmess


  routine = 'write_mesh'
  cmess   = ''

  self%mesh_name        = mesh_name
  self%mesh_class       = mesh_class
  self%generator_inputs = generator_inputs

  self%nmesh_nodes = num_nodes
  self%nmesh_edges = num_edges
  self%nmesh_faces = num_faces


  ! Set up NetCDF header
  call define_dimensions (self)
  call define_variables  (self)
  call assign_attributes (self)


  ! End definitions before putting data in.
  ierr = nf90_enddef(self%ncid)
  call check_err(ierr, routine, cmess)


  ! Node coordinates
  ierr = nf90_put_var( self%ncid, self%mesh_node_x_id, node_coordinates(1,:) )
  call check_err(ierr, routine, cmess)

  ierr = nf90_put_var( self%ncid, self%mesh_node_y_id, node_coordinates(2,:) )
  call check_err(ierr, routine, cmess)

  ! Face node connectivity
  ierr = nf90_put_var( self%ncid, self%mesh_face_nodes_id, &
                       face_node_connectivity(:,:) )
  call check_err(ierr, routine, cmess)

  ! Edge node connectivity
  ierr = nf90_put_var( self%ncid, self%mesh_edge_nodes_id, &
                       edge_node_connectivity(:,:) )
  call check_err(ierr, routine, cmess)

  ! Face edge connectivity
  ierr = nf90_put_var( self%ncid, self%mesh_face_edges_id, &
                       face_edge_connectivity(:,:) )
  call check_err(ierr, routine, cmess)

  ! Face face connectivity
  ierr = nf90_put_var( self%ncid, self%mesh_face_links_id, &
                       face_face_connectivity(:,:) )
  call check_err(ierr, routine, cmess)

  return
end subroutine write_mesh

!-------------------------------------------------------------------------------
!>  @brief Function to determine if mesh is present in NetCDF ugrid file
!>
!>  @param[in]       mesh_name   Name of the mesh topology
!>  @return          answer      .True. if mesh_name is present in file
!-------------------------------------------------------------------------------
function is_mesh_present(self, mesh_name) result(answer)

  implicit none

  class(ncdf_quad_type), intent(in) :: self
  character(str_def),    intent(in) :: mesh_name

  logical(l_def) :: answer

  character(nf90_max_name), allocatable :: mesh_names(:)

  integer(i_def) :: nmeshes
  integer(i_def) :: i

  nmeshes = self%get_nmeshes()
  allocate(mesh_names(nmeshes))
  mesh_names = ''

  call get_mesh_names(self, mesh_names)

  answer = .false.
  do i=1, nmeshes
    if ( trim(mesh_names(i)) == trim(mesh_name) ) then
      answer = .true.
      exit
    end if
  end do

  deallocate(mesh_names)

  return

end function is_mesh_present

!-------------------------------------------------------------------------------
!>  @brief Adds a mesh to an existing NetCDF ugrid file
!>
!>  @param[in,out]  self                     The NetCDF file object.
!>  @param[in]      mesh_name                Name of the mesh topology
!>  @param[in]      mesh_class               Primitive class of mesh.
!>  @param[in]      generator_inputs         Inputs used to create this mesh
!>                                           from the mesh_generator
!>  @param[in]      num_nodes                The number of nodes on the mesh.
!>  @param[in]      num_edges                The number of edges on the mesh.
!>  @param[in]      num_faces                The number of faces on the mesh.
!>  @param[in]      node_coordinates         long/lat coordinates of each node.
!>  @param[in]      face_node_connectivity   Nodes adjoining each face.
!>  @param[in]      edge_node_connectivity   Nodes adjoining each edge.
!>  @param[in]      face_edge_connectivity   Edges adjoining each face.
!>  @param[in]      face_face_connectivity   Faces adjoining each face (links).
!-------------------------------------------------------------------------------
subroutine append_mesh( self, mesh_name, mesh_class, generator_inputs,     &
                        num_nodes, num_edges, num_faces, node_coordinates, &
                        face_node_connectivity, edge_node_connectivity,    &
                        face_edge_connectivity, face_face_connectivity)
  implicit none

  ! Arguments
  class(ncdf_quad_type), intent(inout) :: self
  character(str_def),    intent(in)    :: mesh_class
  character(str_def),    intent(in)    :: mesh_name
  character(str_long),   intent(in)    :: generator_inputs
  integer(i_def),        intent(in)    :: num_nodes
  integer(i_def),        intent(in)    :: num_edges
  integer(i_def),        intent(in)    :: num_faces
  real(r_def),           intent(in)    :: node_coordinates(:,:)
  integer(i_def),        intent(in)    :: face_node_connectivity(:,:)
  integer(i_def),        intent(in)    :: edge_node_connectivity(:,:)
  integer(i_def),        intent(in)    :: face_edge_connectivity(:,:)
  integer(i_def),        intent(in)    :: face_face_connectivity(:,:)

  ! Internal variables
  integer(i_def) :: ierr

  character(str_long) :: routine
  character(str_long) :: cmess

  logical(l_def) :: mesh_present

  routine='append_mesh'
  cmess=''

  mesh_present = self%is_mesh_present(mesh_name)

  if (mesh_present) then
    write(log_scratch_space,'(A)') &
        'Mesh '//trim(mesh_name)//' already used or is not unique.'
    call log_event(log_scratch_space, LOG_LEVEL_ERROR)
    return
  end if

  ierr = nf90_redef(self%ncid)
  call check_err(ierr, routine, cmess)

  call self%write_mesh(                                &
      mesh_name  = mesh_name,                          &
      mesh_class = mesh_class,                         &
      generator_inputs = generator_inputs,             &
      num_nodes  = num_nodes,                          &
      num_edges  = num_edges,                          &
      num_faces  = num_faces,                          &
      node_coordinates       = node_coordinates,       &
      face_node_connectivity = face_node_connectivity, &
      edge_node_connectivity = edge_node_connectivity, &
      face_edge_connectivity = face_edge_connectivity, &
      face_face_connectivity = face_face_connectivity )

  return
end subroutine append_mesh

!-------------------------------------------------------------------------------
!>  @brief   Returns the NetCDF variable names in the NetCDF file which are
!>           ugrid mesh topologies.
!>  @details Scans the variable attributes for cf_role = mesh_topology as
!>           a flag that the variable name relates to a mesh.
!>
!>  @param[in]  self        ncdf_quad_type object associated with NetCDF file
!>
!>  @param[out] mesh_names  Character[:], Names of the mesh_topologies
!>                          in the NetCDF file
!>  @param[out] nmeshes     Integer, The number of mesh topologies
!>                          in the NetCDf file <<optional>>
!-------------------------------------------------------------------------------
subroutine scan_for_topologies(self, mesh_names, nmeshes)

  implicit none

  class(ncdf_quad_type), intent(in)  :: self
  character(len=*),      intent(out) :: mesh_names(:)
  integer(i_def),        intent(out) :: nmeshes

  character(nf90_max_name), allocatable :: var_names(:)
  integer(i_def),           allocatable :: var_n_attributes(:)
  logical,                  allocatable :: is_mesh_topology(:)

  integer(i_def) :: n_variables
  integer(i_def) :: i, j, counter, ierr

  character(nf90_max_name) :: attribute_name
  character(nf90_max_name) :: attribute_value

  integer(i_def) :: n_mesh_topologies

  character(str_long) :: routine
  character(str_long) :: cmess

  routine = 'get_mesh_names'
  cmess = ''

  ierr = nf90_inquire(ncid=self%ncid, nvariables=n_variables )

  allocate(var_names(n_variables))
  allocate(var_n_attributes(n_variables))
  allocate(is_mesh_topology(n_variables))

  is_mesh_topology = .false.

  do i=1, n_variables
    ierr = nf90_inquire_variable( ncid=self%ncid, varid=i, &
                                  name=var_names(i),       &
                                  natts=var_n_attributes(i) )
    write(cmess,'(2(A,I0))')       &
        'Invalid variable id:', i, &
        'or NetCDF file id:', self%ncid
    call check_err (ierr, routine, cmess)
    do j=1, var_n_attributes(i)
      ierr = nf90_inq_attname( ncid=self%ncid,    &
                               varid=i, attnum=j, &
                               name=attribute_name )
      write(cmess,'(A,I0,A)')              &
          'Invalid attribute number: ', j, &
          ' for variable '//trim(var_names(i))
      call check_err (ierr, routine, cmess)

      if (trim(attribute_name) == 'cf_role') then
        ierr = nf90_get_att( ncid=self%ncid,            &
                             varid=i,                   &
                             name=trim(attribute_name), &
                             values=attribute_value )
        write(cmess,'(A)')                                             &
            'Unable to get value of attribute cf_role for variable '// &
            trim(var_names(i))
        call check_err (ierr, routine, cmess)
        if (trim(attribute_value) == 'mesh_topology') then
          is_mesh_topology(i) = .true.
          exit
        end if
      end if

    end do

  end do

  n_mesh_topologies = count(is_mesh_topology)

  if ( size( mesh_names ) < n_mesh_topologies ) then
    write(log_scratch_space,'(I0,A,I0)') n_mesh_topologies,   &
        ' mesh topologies found but output array provided'//  &
        ' is only length ', size(mesh_names)
    call log_event(trim(log_scratch_space), LOG_LEVEL_ERROR)
  end if

  counter=0
  do i=1, n_variables
    if (is_mesh_topology(i)) then
      counter=counter+1
      mesh_names(counter) = trim(var_names(i))
    end if
  end do

  nmeshes = n_mesh_topologies

  deallocate(var_names)
  deallocate(var_n_attributes)
  deallocate(is_mesh_topology)

  return
end subroutine scan_for_topologies

end module ncdf_quad_mod

