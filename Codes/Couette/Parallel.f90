!==================================================================================================
MODULE Parallel		! Defines Parallel (MPI) Variables
							! Contains Parallel (MPI) Subroutines (MPI_Sub_Info, MPI_Initialize, FillSendArrays, MPI_Transfer, SendData, RecvData)
							! Written by Yanxing Wang (2008)
			 				! Modified by Gino Banco (2008-2009)
!================================================================================================== 
USE SetPrecision
USE Setup
USE LBM
USE PassiveScalar
USE MPI					! Intrinsic MPI definitions module

IMPLICIT NONE 

CONTAINS

!--------------------------------------------------------------------------------------------------
SUBROUTINE MPI_Setup	! setup the MPI (parallel) component
!--------------------------------------------------------------------------------------------------
IMPLICIT NONE

INTEGER(lng) :: mpierr! MPI standard error object
! Initialize variables and arrays
f_Comps				= 0_lng		! specifies the components of the distribution functions to transfer in each MPI communication direction
Corner_SendIndex	= 0_lng		! i, j, and k indices for each corner
Corner_RecvIndex	= 0_lng		! i, j, and k indices for each corner (phantom node for recieving data)
Z_SendIndex			= 0_lng		! i and j indices for each Z side 
Z_RecvIndex			= 0_lng		! i and j indices for each Z side (phantom node for recieving data)
X_SendIndex			= 0_lng		! j and k indices for each X side 
X_RecvIndex			= 0_lng		! j and k indices for each X side (phantom node for recieving data)
Y_SendIndex			= 0_lng		! i and k indices for each Y side 
Y_RecvIndex			= 0_lng		! i and k indices for each Y side (phantom node for recieving data)
YZ_SendIndex		= 0_lng		! i index for each YZ face 
YZ_RecvIndex		= 0_lng		! i index for each YZ face (phantom node for recieving data)
ZX_SendIndex		= 0_lng		! j index for each ZX face 
ZX_RecvIndex		= 0_lng		! j index for each ZX face (phantom node for recieving data)
XY_SendIndex		= 0_lng		! k index for each XY face 
XY_RecvIndex		= 0_lng		! k index for each XY face (phantom node for recieving data)
OppCommDir 			= 0_lng		! opposite MPI communication directions (like bounceback) 
CommDataStart_f	= 0_lng		! array of starting indices in the send arrays for the distribution functions from each communication direction 
CommDataStart_rho	= 0_lng		! array of starting indices in the send arrays for the density from each communication direction
CommDataStart_phi	= 0_lng		! array of starting indices in the send arrays for the scalar from each communication direction
CommDataStart_u	= 0_lng		! array of starting indices in the send arrays for the scalar from each communication direction
CommDataStart_v	= 0_lng		! array of starting indices in the send arrays for the scalar from each communication direction
CommDataStart_w	= 0_lng		! array of starting indices in the send arrays for the scalar from each communication direction
CommDataStart_node	= 0_lng		! array of starting indices in the send arrays for the scalar from each communication direction
fSize				= 0_lng		! array of the number of elements sent for each communication direction (distribution functions)
dsSize				= 0_lng		! array of the number of elements sent for each communication direction (density and scalar)
uvwSize				= 0_lng		! array of the number of elements sent for each communication direction (density and scalar)
nodeSize			= 0_lng		! array of the number of elements sent for each communication direction (density and scalar)

! Initialize/setup derived data types
!der_block_len = (/1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1/)
!der_block_types=	(/ MPI_INTEGER, &
!			MPI_DOUBLE_PRECISION, &
!			MPI_DOUBLE_PRECISION, &
!			MPI_DOUBLE_PRECISION, &
!			MPI_DOUBLE_PRECISION, &
!			MPI_DOUBLE_PRECISION, &
!			MPI_DOUBLE_PRECISION, &
!			MPI_DOUBLE_PRECISION, &
!			MPI_DOUBLE_PRECISION, &
!			MPI_DOUBLE_PRECISION, &
!			MPI_DOUBLE_PRECISION, &
!			MPI_DOUBLE_PRECISION, &
!			MPI_DOUBLE_PRECISION, &
!			MPI_DOUBLE_PRECISION, &
!			MPI_DOUBLE_PRECISION, &
!			MPI_DOUBLE_PRECISION, &
!			MPI_DOUBLE_PRECISION, &
!			MPI_DOUBLE_PRECISION, &
!			MPI_DOUBLE_PRECISION, &
!			MPI_DOUBLE_PRECISION, &
!			MPI_INTEGER, &
!			MPI_INTEGER /)

!der_block_len = (/1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1/)
!der_block_types=	(/ MPI_INTEGER, &
!			MPI_INTEGER, &
!			MPI_INTEGER, &
!			MPI_DOUBLE_PRECISION, &
!			MPI_DOUBLE_PRECISION, &
!			MPI_DOUBLE_PRECISION, &
!			MPI_DOUBLE_PRECISION, &
!			MPI_DOUBLE_PRECISION, &
!			MPI_DOUBLE_PRECISION, &
!			MPI_DOUBLE_PRECISION, &
!			MPI_DOUBLE_PRECISION, &
!			MPI_DOUBLE_PRECISION, &
!			MPI_DOUBLE_PRECISION, &
!			MPI_DOUBLE_PRECISION, &
!			MPI_DOUBLE_PRECISION, &
!			MPI_DOUBLE_PRECISION, &
!			MPI_DOUBLE_PRECISION, &
!			MPI_DOUBLE_PRECISION, &
!			MPI_DOUBLE_PRECISION, &
!			MPI_DOUBLE_PRECISION, &
!			MPI_DOUBLE_PRECISION, &
!			MPI_DOUBLE_PRECISION/)
!CALL MPI_TYPE_EXTENT(MPI_DOUBLE_PRECISION,mpidblextent,mpierr)
!CALL MPI_TYPE_EXTENT(MPI_INTEGER,mpiintextent,mpierr)
!der_block_offsets=	(/ 0, &
!			mpiintextent, &
!			2*mpiintextent, &
!			3*mpiintextent+0*mpidblextent, &
!			3*mpiintextent+1*mpidblextent, &
!			3*mpiintextent+2*mpidblextent, &
!			3*mpiintextent+3*mpidblextent, &
!			3*mpiintextent+4*mpidblextent, &
!			3*mpiintextent+5*mpidblextent, &
!			3*mpiintextent+6*mpidblextent, &
!			3*mpiintextent+7*mpidblextent, &
!			3*mpiintextent+8*mpidblextent, &
!			3*mpiintextent+9*mpidblextent, &
!			3*mpiintextent+10*mpidblextent, &
!			3*mpiintextent+11*mpidblextent, &
!			3*mpiintextent+12*mpidblextent, &
!			3*mpiintextent+13*mpidblextent, &
!			3*mpiintextent+14*mpidblextent, &
!			3*mpiintextent+15*mpidblextent, &
!			3*mpiintextent+16*mpidblextent, &
!			3*mpiintextent+17*mpidblextent, &
!			3*mpiintextent+18*mpidblextent/)


der_block_len = (/1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1/)
der_block_types=	(/ MPI_INTEGER, &
			MPI_INTEGER, &
			MPI_INTEGER, &
			MPI_DOUBLE_PRECISION, &
			MPI_DOUBLE_PRECISION, &
			MPI_DOUBLE_PRECISION, &
			MPI_DOUBLE_PRECISION, &
			MPI_DOUBLE_PRECISION, &
			MPI_DOUBLE_PRECISION, &
			MPI_DOUBLE_PRECISION, &
			MPI_DOUBLE_PRECISION, &
			MPI_DOUBLE_PRECISION, &
			MPI_DOUBLE_PRECISION, &
			MPI_DOUBLE_PRECISION, &
			MPI_DOUBLE_PRECISION, &
			MPI_DOUBLE_PRECISION, &
			MPI_DOUBLE_PRECISION, &
			MPI_DOUBLE_PRECISION, &
			MPI_DOUBLE_PRECISION, &
			MPI_DOUBLE_PRECISION, &
			MPI_DOUBLE_PRECISION, &
			MPI_DOUBLE_PRECISION, &
			MPI_DOUBLE_PRECISION, &
			MPI_DOUBLE_PRECISION, &
			MPI_DOUBLE_PRECISION, &
			MPI_DOUBLE_PRECISION/)
CALL MPI_TYPE_EXTENT(MPI_DOUBLE_PRECISION,mpidblextent,mpierr)
CALL MPI_TYPE_EXTENT(MPI_INTEGER,mpiintextent,mpierr)
der_block_offsets=	(/ 0, &
			mpiintextent, &
			2*mpiintextent, &
			3*mpiintextent+0*mpidblextent, &
			3*mpiintextent+1*mpidblextent, &
			3*mpiintextent+2*mpidblextent, &
			3*mpiintextent+3*mpidblextent, &
			3*mpiintextent+4*mpidblextent, &
			3*mpiintextent+5*mpidblextent, &
			3*mpiintextent+6*mpidblextent, &
			3*mpiintextent+7*mpidblextent, &
			3*mpiintextent+8*mpidblextent, &
			3*mpiintextent+9*mpidblextent, &
			3*mpiintextent+10*mpidblextent, &
			3*mpiintextent+11*mpidblextent, &
			3*mpiintextent+12*mpidblextent, &
			3*mpiintextent+13*mpidblextent, &
			3*mpiintextent+14*mpidblextent, &
			3*mpiintextent+15*mpidblextent, &
			3*mpiintextent+16*mpidblextent, &
			3*mpiintextent+17*mpidblextent, &
			3*mpiintextent+18*mpidblextent, &
			3*mpiintextent+19*mpidblextent, &
			3*mpiintextent+20*mpidblextent, &
			3*mpiintextent+21*mpidblextent, &
			3*mpiintextent+22*mpidblextent/)

CALL MPI_TYPE_STRUCT(der_type_count, der_block_len, der_block_offsets, der_block_types,mpipartransfertype, mpierr)
!CALL MPI_TYPE_CREATE_STRUCT(der_type_count, der_block_len, der_block_offsets, der_block_types,mpipartransfertype,mpierr)
CALL MPI_TYPE_COMMIT(mpipartransfertype, mpierr)


! Fill out the MPI arrays
CALL MPI_Initialize

!------------------------------------------------
END SUBROUTINE MPI_Setup
!------------------------------------------------

!--------------------------------------------------------------------------------------------------
SUBROUTINE MPI_Initialize	! initialize the MPI arrays and variables
!--------------------------------------------------------------------------------------------------
IMPLICIT NONE

! Define local variables
INTEGER(lng) :: mpierr! MPI standard error object
INTEGER(lng) :: iComm												! index variable
INTEGER(lng) :: YZ_FaceSize, ZX_FaceSize, XY_FaceSize		! number of nodes on a subdomain face oriented on the respective faces
INTEGER(lng) :: f_SendSize, ds_SendSize, node_Sendsize, uvw_SendSize, total_SendSize	! sizes of the distribution function, density, velocity, and scalar data transfer arrays respectively
!INTEGER(lng) :: msgSze												! size of the message send/recv arrays (maximum size)

! Opposite directions for communication - same concept as bounceback directions, but there are 26 communication directions whereas there are only 14 non-stationary distribution function directions
OppCommDir(1)	=	2_lng
OppCommDir(2)	=	1_lng

OppCommDir(3)	=	4_lng
OppCommDir(4)	=	3_lng

OppCommDir(5)	=	6_lng
OppCommDir(6)	=	5_lng

OppCommDir(7)	=	8_lng
OppCommDir(8)	=	7_lng

OppCommDir(9)	=	10_lng
OppCommDir(10)	=	9_lng

OppCommDir(11)	=	12_lng
OppCommDir(12)	=	11_lng

OppCommDir(13)	=	14_lng
OppCommDir(14)	=	13_lng

OppCommDir(15)	=	16_lng
OppCommDir(16)	=	15_lng

OppCommDir(17)	=	18_lng
OppCommDir(18)	=	17_lng

OppCommDir(19)	=	20_lng
OppCommDir(20)	=	19_lng

OppCommDir(21)	=	22_lng
OppCommDir(22)	=	21_lng

OppCommDir(23)	=	24_lng
OppCommDir(24)	=	23_lng

OppCommDir(25)	=	26_lng
OppCommDir(26)	=	25_lng


! Distribution function components transferred
f_Comps = 0_lng		! initialize the array

! Faces
f_Comps(1,1) = 1
f_Comps(1,2) = 7
f_Comps(1,3) = 9
f_Comps(1,4) = 12
f_Comps(1,5) = 13

f_Comps(2,1) = 2
f_Comps(2,2) = 8
f_Comps(2,3) = 10
f_Comps(2,4) = 11
f_Comps(2,5) = 14

f_Comps(3,1) = 3
f_Comps(3,2) = 7
f_Comps(3,3) = 9
f_Comps(3,4) = 11
f_Comps(3,5) = 14

f_Comps(4,1) = 4
f_Comps(4,2) = 8
f_Comps(4,3) = 10
f_Comps(4,4) = 12
f_Comps(4,5) = 13

f_Comps(5,1) = 5
f_Comps(5,2) = 7
f_Comps(5,3) = 10
f_Comps(5,4) = 11
f_Comps(5,5) = 13

f_Comps(6,1) = 6
f_Comps(6,2) = 8
f_Comps(6,3) = 9
f_Comps(6,4) = 12
f_Comps(6,5) = 14

! Sides
f_Comps(7,1) = 7
f_Comps(7,2) = 9

f_Comps(8,1) = 8
f_Comps(8,2) = 10

f_Comps(9,1) = 12
f_Comps(9,2) = 13 

f_Comps(10,1) = 11
f_Comps(10,2) = 14

f_Comps(11,1) = 7
f_Comps(11,2) = 11

f_Comps(12,1) = 8
f_Comps(12,2) = 12

f_Comps(13,1) = 9
f_Comps(13,2) = 14

f_Comps(14,1) = 10
f_Comps(14,2) = 13

f_Comps(15,1) = 7
f_Comps(15,2) = 13

f_Comps(16,1) = 8
f_Comps(16,2) = 14

f_Comps(17,1) = 10
f_Comps(17,2) = 11

f_Comps(18,1) = 9
f_Comps(18,2) = 12

! Corners
f_Comps(19,1) = 7
f_Comps(20,1) = 8
f_Comps(21,1) = 9
f_Comps(22,1) = 10
f_Comps(23,1) = 11
f_Comps(24,1) = 12
f_Comps(25,1) = 13
f_Comps(26,1) = 14


! Fill out the size arrays
YZ_FaceSize		= nySub*nzSub							! number of nodes on a subdomain face oriented in the ZY plane
ZX_FaceSize		= nzSub*nxSub							! number of nodes on a subdomain face oriented in the ZX plane
XY_FaceSize		= nxSub*nySub							! number of nodes on a subdomain face oriented in the XY plane

fSize 			= 0_lng									! initialize array
dsSize 			= 0_lng									! initialize array
uvwSize			= 0_lng									! initialize array
nodeSize		= 0_lng									! initialize array

fSize(1:2)		= YZ_FaceSize*NumFs_face			! YZ faces
fSize(3:4)		= ZX_FaceSize*NumFs_face			! ZX faces
fSize(5:6)		= XY_FaceSize*NumFs_face			! XY faces
fSize(7:10)		= nzSub*NumFs_side					! Z sides
fSize(11:14)	= nxSub*NumFs_side					! X sides
fSize(15:18)	= nySub*NumFs_side					! Y sides
fSize(19:26)	= 1_lng*NumFs_corner					! corners

dsSize(1:2)		= YZ_FaceSize							! YZ faces
dsSize(3:4)		= ZX_FaceSize							! ZX faces
dsSize(5:6)		= XY_FaceSize							! XY faces
dsSize(7:10)	= nzSub									! Z sides
dsSize(11:14)	= nxSub									! X sides
dsSize(15:18)	= nySub									! Y sides
dsSize(19:26)	= 1_lng									! corners


uvwSize(1:2)		= YZ_FaceSize							! YZ faces
uvwSize(3:4)		= ZX_FaceSize							! ZX faces
uvwSize(5:6)		= XY_FaceSize							! XY faces
uvwSize(7:10)	= nzSub									! Z sides
uvwSize(11:14)	= nxSub									! X sides
uvwSize(15:18)	= nySub									! Y sides
uvwSize(19:26)	= 1_lng									! corners

nodeSize(1:2)		= YZ_FaceSize							! YZ faces
nodeSize(3:4)		= ZX_FaceSize							! ZX faces
nodeSize(5:6)		= XY_FaceSize							! XY faces
nodeSize(7:10)	= nzSub									! Z sides
nodeSize(11:14)	= nxSub									! X sides
nodeSize(15:18)	= nySub									! Y sides
nodeSize(19:26)	= 1_lng									! corners

msgSize(:)		= fSize(:) + 2_lng*(dsSize(:))+3_lng*(uvwSize(:))+1_lng*(nodeSize(:))	! total message sizes

f_SendSize	= SUM(fSize)
ds_SendSize	= SUM(dsSize)
uvw_SendSize	= SUM(uvwSize)
node_SendSize	= SUM(nodeSize)
total_SendSize  = SUM(msgSize)

ALLOCATE(msgSend(total_SendSize))						
ALLOCATE(msgRecv(total_SendSize))

! Fill out the '3D_Index' arrays (for converting back from 1D array)
! Faces
YZ_SendIndex(1)	= nxSub
YZ_RecvIndex(1)	= nxSub+1_lng		

YZ_SendIndex(2)	= 1_lng	
YZ_RecvIndex(2)	= 0_lng	

ZX_SendIndex(3)	= nySub		
ZX_RecvIndex(3)	= nySub+1_lng		

ZX_SendIndex(4)	= 1_lng	
ZX_RecvIndex(4)	= 0_lng	

XY_SendIndex(5)	= nzSub
XY_RecvIndex(5)	= nzSub+1_lng		
		
XY_SendIndex(6)	= 1_lng
XY_RecvIndex(6)	= 0_lng	


! Sides
! Z Sides
Z_SendIndex(7,1)	= nxSub
Z_RecvIndex(7,1)	= nxSub+1_lng

Z_SendIndex(7,2)	= nySub
Z_RecvIndex(7,2)	= nySub+1_lng


Z_SendIndex(8,1)	= 1_lng
Z_RecvIndex(8,1)	= 0_lng

Z_SendIndex(8,2)	= 1_lng
Z_RecvIndex(8,2)	= 0_lng


Z_SendIndex(9,1)	= nxSub
Z_RecvIndex(9,1)	= nxSub+1_lng

Z_SendIndex(9,2)	= 1_lng
Z_RecvIndex(9,2)	= 0_lng


Z_SendIndex(10,1)	= 1_lng
Z_RecvIndex(10,1)	= 0_lng

Z_SendIndex(10,2)	= nySub
Z_RecvIndex(10,2)	= nySub+1_lng


! X Sides
X_SendIndex(11,1)	= nySub
X_RecvIndex(11,1)	= nySub+1_lng

X_SendIndex(11,2)	= nzSub
X_RecvIndex(11,2)	= nzSub+1_lng


X_SendIndex(12,1)	= 1_lng
X_RecvIndex(12,1)	= 0_lng

X_SendIndex(12,2)	= 1_lng
X_RecvIndex(12,2)	= 0_lng


X_SendIndex(13,1)	= nySub
X_RecvIndex(13,1)	= nySub+1_lng

X_SendIndex(13,2)	= 1_lng
X_RecvIndex(13,2)	= 0_lng


X_SendIndex(14,1)	= 1_lng
X_RecvIndex(14,1)	= 0_lng

X_SendIndex(14,2)	= nzSub
X_RecvIndex(14,2)	= nzSub+1_lng


! Y Sides
Y_SendIndex(15,1)	= nxSub
Y_RecvIndex(15,1)	= nxSub+1_lng

Y_SendIndex(15,2)	= nzSub
Y_RecvIndex(15,2)	= nzSub+1_lng


Y_SendIndex(16,1)	= 1_lng
Y_RecvIndex(16,1)	= 0_lng

Y_SendIndex(16,2)	= 1_lng
Y_RecvIndex(16,2)	= 0_lng


Y_SendIndex(17,1)	= 1
Y_RecvIndex(17,1)	= 0_lng

Y_SendIndex(17,2)	= nzSub
Y_RecvIndex(17,2)	= nzSub+1_lng


Y_SendIndex(18,1)	= nxSub
Y_RecvIndex(18,1)	= nxSub+1_lng

Y_SendIndex(18,2)	= 1_lng
Y_RecvIndex(18,2)	= 0_lng


! Corners
Corner_SendIndex(19,1) = nxSub
Corner_SendIndex(19,2) = nySub
Corner_SendIndex(19,3) = nzSub

Corner_RecvIndex(19,1) = nxSub+1_lng
Corner_RecvIndex(19,2) = nySub+1_lng
Corner_RecvIndex(19,3) = nzSub+1_lng


Corner_SendIndex(20,1) = 1_lng
Corner_SendIndex(20,2) = 1_lng
Corner_SendIndex(20,3) = 1_lng

Corner_RecvIndex(20,1) = 0_lng
Corner_RecvIndex(20,2) = 0_lng
Corner_RecvIndex(20,3) = 0_lng


Corner_SendIndex(21,1) = nxSub
Corner_SendIndex(21,2) = nySub
Corner_SendIndex(21,3) = 1_lng

Corner_RecvIndex(21,1) = nxSub+1_lng
Corner_RecvIndex(21,2) = nySub+1_lng
Corner_RecvIndex(21,3) = 0_lng


Corner_SendIndex(22,1) = 1_lng
Corner_SendIndex(22,2) = 1_lng
Corner_SendIndex(22,3) = nzSub

Corner_RecvIndex(22,1) = 0_lng
Corner_RecvIndex(22,2) = 0_lng
Corner_RecvIndex(22,3) = nzSub+1_lng


Corner_SendIndex(23,1) = 1_lng
Corner_SendIndex(23,2) = nySub
Corner_SendIndex(23,3) = nzSub

Corner_RecvIndex(23,1) = 0_lng
Corner_RecvIndex(23,2) = nySub+1_lng
Corner_RecvIndex(23,3) = nzSub+1_lng


Corner_SendIndex(24,1) = nxSub
Corner_SendIndex(24,2) = 1_lng
Corner_SendIndex(24,3) = 1_lng

Corner_RecvIndex(24,1) = nxSub+1_lng
Corner_RecvIndex(24,2) = 0_lng
Corner_RecvIndex(24,3) = 0_lng


Corner_SendIndex(25,1) = nxSub
Corner_SendIndex(25,2) = 1_lng
Corner_SendIndex(25,3) = nzSub

Corner_RecvIndex(25,1) = nxSub+1_lng
Corner_RecvIndex(25,2) = 0_lng
Corner_RecvIndex(25,3) = nzSub+1_lng


Corner_SendIndex(26,1) = 1_lng
Corner_SendIndex(26,2) = nySub
Corner_SendIndex(26,3) = 1_lng

Corner_RecvIndex(26,1) = 0_lng
Corner_RecvIndex(26,2) = nySub+1_lng
Corner_RecvIndex(26,3) = 0_lng


! Fill out the 'CommDataStart' arrays
! Initialize arrays
CommDataStart_f  	= 0_lng					! distribution functions
CommDataStart_rho  	= 0_lng					! density
CommDataStart_phi  	= 0_lng					! scalar
CommDataStart_u  	= 0_lng					! velocity
CommDataStart_v  	= 0_lng					! velocity
CommDataStart_w  	= 0_lng					! velocity
CommDataStart_node  	= 0_lng					! node

CommDataStart_f(1)	= 1_lng	
CommDataStart_rho(1) = CommDataStart_f(1)	+ fSize(1)
CommDataStart_phi(1) = CommDataStart_rho(1)	+ dsSize(1)
CommDataStart_u(1) = CommDataStart_phi(1)	+ dsSize(1)
CommDataStart_v(1) = CommDataStart_u(1)		+ uvwSize(1)
CommDataStart_w(1) = CommDataStart_v(1)		+ uvwSize(1)
CommDataStart_node(1) = CommDataStart_w(1)	+ uvwSize(1)

DO iComm=2,NumCommDirs							! fill out for communication directions 2-NumCommDirs
  CommDataStart_f(iComm)	= CommDataStart_node(iComm-1) 	+ nodeSize(iComm-1)
  !CommDataStart_f(iComm)	= CommDataStart_w(iComm-1) 	+ uvwSize(iComm-1)
  CommDataStart_rho(iComm)	= CommDataStart_f(iComm)	+ fSize(iComm)
  CommDataStart_phi(iComm)	= CommDataStart_rho(iComm) 	+ dsSize(iComm)
  CommDataStart_u(iComm)	= CommDataStart_phi(iComm) 	+ dsSize(iComm)
  CommDataStart_v(iComm)	= CommDataStart_u(iComm) 	+ uvwSize(iComm)
  CommDataStart_w(iComm)	= CommDataStart_v(iComm) 	+ uvwSize(iComm)
  CommDataStart_node(iComm)	= CommDataStart_w(iComm) 	+ uvwSize(iComm)
END DO

!WRITE(6678,*) 'f_SendSize', total_SendSize
!WRITE(6678,*) 'ds_SendSize', total_SendSize
!WRITE(6678,*) 'total_SendSize', total_SendSize
!DO iComm=1,NumCommDirs
!  WRITE(6678,*) 'CommDataStart_f(iComm)', CommDataStart_f(iComm)
!  WRITE(6678,*) 'CommDataStart_rho(iComm)', CommDataStart_rho(iComm) 
!  WRITE(6678,*) 'CommDataStart_phi(iComm)', CommDataStart_phi(iComm) 
!END DO
!STOP

! Allocate the MPI_WAITALL status array
ALLOCATE(waitStat(MPI_STATUS_SIZE,2*NumCommDirs))

! ALLOCATE Particle MPI transfer variables
!ALLOCATE(partransfersend(NumCommDirsPar,NumParVar))
ALLOCATE(partransfersend(NumParVar,NumCommDirsPar))
!ALLOCATE(partransferrecv(NumCommDirsPar,NumParVar))
ALLOCATE(partransferrecv(NumParVar,NumCommDirsPar))
ALLOCATE(numpartransfer(NumCommDirsPar))
ALLOCATE(parreqid(2*NumCommDirsPar))
ALLOCATE(parwtstat(MPI_STATUS_SIZE,2*NumCommDirsPar))
ALLOCATE(probestat(MPI_STATUS_SIZE))
!ALLOCATE(ParSendArray(1,NumCommDirsPar))
!ALLOCATE(ParRecvArray(1,NumCommDirsPar))

ParInit%parid = 0_lng
ParInit%cur_part = 0_lng
ParInit%new_part = 0_lng
ParInit%xp = 0.0_dbl
ParInit%yp = 0.0_dbl
ParInit%zp = 0.0_dbl
ParInit%up = 0.0_dbl
ParInit%vp = 0.0_dbl
ParInit%wp = 0.0_dbl
ParInit%rp = 0.0_dbl
ParInit%xpold = 0.0_dbl
ParInit%ypold = 0.0_dbl
ParInit%zpold = 0.0_dbl
ParInit%upold = 0.0_dbl
ParInit%vpold = 0.0_dbl
ParInit%wpold = 0.0_dbl
ParInit%rpold = 0.0_dbl
ParInit%par_conc = 0.0_dbl
ParInit%bulk_conc = 0.0_dbl
ParInit%sh = 0.0_dbl
ParInit%S = 0.0_dbl
ParInit%Sst = 0.0_dbl
ParInit%Veff = 0.0_dbl
ParInit%Nbj = 0.0_dbl
ParInit%gamma_cont = 0.0_dbl
ParInit%delNBbyCV = 0.0_dbl

!ALLOCATE(partransfersend(NumParVar))
!ALLOCATE(partransferrecv(NumParVar))
!ALLOCATE(parreqid(2))
!ALLOCATE(parwtstat(MPI_STATUS_SIZE,2))
!ALLOCATE(probestat(MPI_STATUS_SIZE))

!------------------------------------------------
END SUBROUTINE MPI_Initialize
!------------------------------------------------

!--------------------------------------------------------------------------------------------------
SUBROUTINE PackData	! transfer the distribution functions between neighboring initialize the MPI arrays and variables
!--------------------------------------------------------------------------------------------------
IMPLICIT NONE

! Define local variables
INTEGER(lng) :: i,j,k,m,iComm,ii	! index variables

! Initialize array
msgSend		= 0.0_dbl

! Fill Arrays
! FACES
! YZ Faces
DO iComm=1,2

  IF(SubID(iComm) .NE. 0) THEN							! only fill if the neighbor contains fluid (otherwise, no data transfer is necessary)

    i = YZ_SendIndex(iComm)								! i index for assigning the proper component of the distribution function to the 1D send array
  
    DO k=1,nzSub
      DO j=1,nySub

        ! density and scalar
        ii = CommDataStart_rho(iComm)			&		! start location
           + (j - 1_lng)							&		! convert 3D-coordinate into proper 1D-array coordinate
           + (k - 1_lng)*nySub							! convert 3D-coordinate into proper 1D-array coordinate

        msgSend(ii) 			= rho(i,j,k)	! store the proper density in the send array
        msgSend(ii+dsSize(iComm))	= phi(i,j,k)	! store the proper scalar quantity in the send array
        msgSend(ii+2_lng*dsSize(iComm))	= u(i,j,k)	! store the proper velocity quantity in the send array
        msgSend(ii+2_lng*dsSize(iComm)+uvwSize(icomm))= v(i,j,k)	! store the proper velocity quantity in the send array
        msgSend(ii+2_lng*dsSize(iComm)+2_lng*uvwSize(icomm))= w(i,j,k)	! store the proper velocity quantity in the send array
        msgSend(ii+2_lng*dsSize(iComm)+3_lng*uvwSize(icomm))= REAL(node(i,j,k),dbl)	! store the proper node flag quantity in the send array

        ! distribution functions
        DO m=1,NumFs_face
     
          ii = CommDataStart_f(iComm)			&		! start location
             + (m - 1_lng)							&		! current distribution function location
             + (j - 1_lng)*NumFs_face			&		! convert 3D-coordinate into proper 1D-array coordinate
             + (k - 1_lng)*NumFs_face*nySub			! convert 3D-coordinate into proper 1D-array coordinate

          msgSend(ii) = f(f_Comps(iComm,m),i,j,k)	! store the proper component of the distribution function in the send array

        END DO   

      END DO
    END DO

  END IF

END DO


! ZX Faces
DO iComm=3,4

  IF(SubID(iComm) .NE. 0) THEN							! only fill if the neighbor contains fluid (otherwise, no data transfer is necessary)

    j = ZX_SendIndex(iComm)								! j index for assigning the proper component of the distribution function to the 1D send array
  
    DO i=1,nxSub
      DO k=1,nzSub

        ! density and scalar
        ii = CommDataStart_rho(iComm)			&		! start location
           + (k - 1_lng)							&		! convert 3D-coordinate into proper 1D-array coordinate
           + (i - 1_lng)*nzSub							! convert 3D-coordinate into proper 1D-array coordinate
    
        msgSend(ii) 						= rho(i,j,k)	! store the proper density in the send array
        msgSend(ii+dsSize(iComm))	= phi(i,j,k)	! store the proper scalar quantity in the send array
        msgSend(ii+2_lng*dsSize(iComm))	= u(i,j,k)	! store the proper velocity quantity in the send array
        msgSend(ii+2_lng*dsSize(iComm)+uvwSize(icomm))= v(i,j,k)	! store the proper velocity quantity in the send array
        msgSend(ii+2_lng*dsSize(iComm)+2_lng*uvwSize(icomm))= w(i,j,k)	! store the proper velocity quantity in the send array
        msgSend(ii+2_lng*dsSize(iComm)+3_lng*uvwSize(icomm))= REAL(node(i,j,k),dbl)	! store the proper node flag quantity in the send array

        ! distribution functions
        DO m=1,NumFs_face
    
          ii = CommDataStart_f(iComm)			&		! start location
             + (m - 1_lng)							&		! current distribution function location
             + (k - 1_lng)*NumFs_face			&		! convert 3D-coordinate into proper 1D-array coordinate
             + (i - 1_lng)*NumFs_face*nzSub			! convert 3D-coordinate into proper 1D-array coordinate

          msgSend(ii) = f(f_Comps(iComm,m),i,j,k)	! store the proper component of the distribution function in the send array

        END DO      

      END DO
    END DO

  END IF  

END DO

! XY Faces
DO iComm=5,6

  IF(SubID(iComm) .NE. 0) THEN							! only fill if the neighbor contains fluid (otherwise, no data transfer is necessary)  

    k = XY_SendIndex(iComm)								! j index for assigning the proper component of the distribution function to the 1D send array
  
    DO j=1,nySub
      DO i=1,nxSub

        ! density and scalar
        ii = CommDataStart_rho(iComm)			&		! start location
           + (i - 1_lng)							&		! convert 3D-coordinate into proper 1D-array coordinate
           + (j - 1_lng)*nxSub							! convert 3D-coordinate into proper 1D-array coordinate     
    
        msgSend(ii) 						= rho(i,j,k)	! store the proper density in the send array
        msgSend(ii+dsSize(iComm))	= phi(i,j,k)	! store the proper scalar quantity in the send array
        msgSend(ii+2_lng*dsSize(iComm))	= u(i,j,k)	! store the proper velocity quantity in the send array
        msgSend(ii+2_lng*dsSize(iComm)+uvwSize(icomm))= v(i,j,k)	! store the proper velocity quantity in the send array
        msgSend(ii+2_lng*dsSize(iComm)+2_lng*uvwSize(icomm))= w(i,j,k)	! store the proper velocity quantity in the send array
        msgSend(ii+2_lng*dsSize(iComm)+3_lng*uvwSize(icomm))= REAL(node(i,j,k),dbl)	! store the proper node flag quantity in the send array
      
        DO m=1,NumFs_face
      
          ii = CommDataStart_f(iComm)		&			! start location
             + (m - 1_lng)						&			! current distribution function location
             + (i - 1_lng)*NumFs_face		&			! convert 3D-coordinate into proper 1D-array coordinate
             + (j - 1_lng)*NumFs_face*nxSub			! convert 3D-coordinate into proper 1D-array coordinate

          msgSend(ii) = f(f_Comps(iComm,m),i,j,k)	! store the proper component of the distribution function in the send array
  
        END DO
      
      END DO
    END DO
  
  END IF

END DO


! SIDES
! Z Sides
DO iComm=7,10

  IF(SubID(iComm) .NE. 0) THEN							! only fill if the neighbor contains fluid (otherwise, no data transfer is necessary)    

    i = Z_SendIndex(iComm,1)								! i index for assigning the proper component of the distribution function to the 1D send array        
    j = Z_SendIndex(iComm,2)								! j index for assigning the proper component of the distribution function to the 1D send array
  
    DO k=1,nzSub

      	! density and scalar
	ii = CommDataStart_rho(iComm)			&			! start location
           + (k - 1_lng)										! convert 3D-coordinate into proper 1D-array coordinate
    
	msgSend(ii) 					= rho(i,j,k)		! store the proper density in the send array
	msgSend(ii+dsSize(iComm))	= phi(i,j,k)		! store the proper scalar quantity in the send array
        msgSend(ii+2_lng*dsSize(iComm))	= u(i,j,k)	! store the proper velocity quantity in the send array
        msgSend(ii+2_lng*dsSize(iComm)+uvwSize(icomm))= v(i,j,k)	! store the proper velocity quantity in the send array
        msgSend(ii+2_lng*dsSize(iComm)+2_lng*uvwSize(icomm))= w(i,j,k)	! store the proper velocity quantity in the send array
        msgSend(ii+2_lng*dsSize(iComm)+3_lng*uvwSize(icomm))= REAL(node(i,j,k),dbl)	! store the proper node flag quantity in the send array
    
      DO m=1,NumFs_side
  
        ! distribution functions   
        ii = CommDataStart_f(iComm)			&			! start location
           + (m - 1_lng)						&			! current distribution function location
           + (k - 1_lng)*NumFs_side						! convert 3D-coordinate into proper 1D-array coordinate

        msgSend(ii) = f(f_Comps(iComm,m),i,j,k)		! store the proper component of the distribution function in the send array
 
      END DO
  
    END DO

  END IF  
  
END DO

! X Sides
DO iComm=11,14

  IF(SubID(iComm) .NE. 0) THEN							! only fill if the neighbor contains fluid (otherwise, no data transfer is necessary)     

    j = X_SendIndex(iComm,1)								! j index for assigning the proper component of the distribution function to the 1D send array        
    k = X_SendIndex(iComm,2)								! k index for assigning the proper component of the distribution function to the 1D send array
  
    DO i=1,nxSub

      ! density and scalar
      ii = CommDataStart_rho(iComm)			&			! start location
         + (i - 1_lng)										! convert 3D-coordinate into proper 1D-array coordinate
    
      msgSend(ii) 					= rho(i,j,k)		! store the proper density in the send array
      msgSend(ii+dsSize(iComm))	= phi(i,j,k)		! store the proper scalar quantity in the send array 
      msgSend(ii+2_lng*dsSize(iComm))	= u(i,j,k)	! store the proper velocity quantity in the send array
      msgSend(ii+2_lng*dsSize(iComm)+uvwSize(icomm))= v(i,j,k)	! store the proper velocity quantity in the send array
      msgSend(ii+2_lng*dsSize(iComm)+2_lng*uvwSize(icomm))= w(i,j,k)	! store the proper velocity quantity in the send array
      msgSend(ii+2_lng*dsSize(iComm)+3_lng*uvwSize(icomm))= REAL(node(i,j,k),dbl)	! store the proper node flag quantity in the send array
    
      DO m=1,NumFs_side

        ! distribution functions      
        ii = CommDataStart_f(iComm)			&			! start location
           + (m - 1_lng)						&			! current distribution function location
           + (i - 1_lng)*NumFs_side						! convert 3D-coordinate into proper 1D-array coordinate

        msgSend(ii) = f(f_Comps(iComm,m),i,j,k)		! store the proper component of the distribution function in the send array

      END DO
    
    END DO

  END IF
  
END DO

! Y Sides
DO iComm=15,18

  IF(SubID(iComm) .NE. 0) THEN							! only fill if the neighbor contains fluid (otherwise, no data transfer is necessary)   

    i = Y_SendIndex(iComm,1)								! i index for assigning the proper component of the distribution function to the 1D send array        
    k = Y_SendIndex(iComm,2)								! k index for assigning the proper component of the distribution function to the 1D send array
  
    DO j=1,nySub

      ! density and scalar
      ii = CommDataStart_rho(iComm)			&			! start location
         + (j - 1_lng)										! convert 3D-coordinate into proper 1D-array coordinate
    
        msgSend(ii) 						= rho(i,j,k)	! store the proper density in the send array
        msgSend(ii+dsSize(iComm))	= phi(i,j,k)	! store the proper scalar quantity in the send array
        msgSend(ii+2_lng*dsSize(iComm))	= u(i,j,k)	! store the proper velocity quantity in the send array
        msgSend(ii+2_lng*dsSize(iComm)+uvwSize(icomm))= v(i,j,k)	! store the proper velocity quantity in the send array
        msgSend(ii+2_lng*dsSize(iComm)+2_lng*uvwSize(icomm))= w(i,j,k)	! store the proper velocity quantity in the send array
        msgSend(ii+2_lng*dsSize(iComm)+3_lng*uvwSize(icomm))= REAL(node(i,j,k),dbl)	! store the proper node flag quantity in the send array
    
      DO m=1,NumFs_side

        ! distribution functions
        ii = CommDataStart_f(iComm)			&			! start location
           + (m - 1_lng)						&			! current distribution function locationcurrent distribution function location
           + (j - 1_lng)*NumFs_side						! convert 3D-coordinate into proper 1D-array coordinate

        msgSend(ii) = f(f_Comps(iComm,m),i,j,k)		! store the proper component of the distribution function in the send array

      END DO
     
    END DO
  
  END IF

END DO


! CORNERS
DO iComm=19,26
  
  IF(SubID(iComm) .NE. 0) THEN							! only fill if the neighbor contains fluid (otherwise, no data transfer is necessary)   
  
    i = Corner_SendIndex(iComm,1)						! j index for assigning the proper component of the distribution function to the 1D send array 	
    j = Corner_SendIndex(iComm,2)						! j index for assigning the proper component of the distribution function to the 1D send array        
    k = Corner_SendIndex(iComm,3)						! k index for assigning the proper component of the distribution function to the 1D send array

    ii = CommDataStart_rho(iComm)						! start location
    msgSend(ii) 		= rho(i,j,k)			! store the proper density in the send array
    msgSend(ii+dsSize(iComm))	= phi(i,j,k)			! store the proper scalar quantity in the send array
    msgSend(ii+2_lng*dsSize(iComm))	= u(i,j,k)	! store the proper velocity quantity in the send array
    msgSend(ii+2_lng*dsSize(iComm)+uvwSize(icomm))= v(i,j,k)	! store the proper velocity quantity in the send array
    msgSend(ii+2_lng*dsSize(iComm)+2_lng*uvwSize(icomm))= w(i,j,k)	! store the proper velocity quantity in the send array
    msgSend(ii+2_lng*dsSize(iComm)+3_lng*uvwSize(icomm))= REAL(node(i,j,k),dbl)	! store the proper node flag quantity in the send array
  
    ii = CommDataStart_f(iComm)							! start location
    msgSend(ii) = f(f_Comps(iComm,1),i,j,k)			! store the proper component of the distribution function in the send array

  END IF
  	
END DO

!------------------------------------------------
END SUBROUTINE PackData
!------------------------------------------------
!--------------------------------------------------------------------------------------------------
SUBROUTINE  Collect_Distribute_Global_Bulk_Scalar_Conc
! This subroutine does the
!parallel communication needed to collect all the bulk concentration in each of
!the processes and computeѕ an average, which is the distributed to the various
!process for computing the drug concentration. 
!--------------------------------------------------------------------------------------------------
IMPLICIT NONE
REAL(dbl) :: Cb_global_temp
INTEGER(lng):: Cb_numFluids_temp,mpierr
!CALL MPI_BARRIER(MPI_COMM_WORLD,mpierr)					! synchronize all processing units before next loop [Intrinsic]
!CALL MPI_GATHER......
!write(*,*)  'In Collect and Distribute before reduce', iter,mySub,Cb_global,Cb_numFluids
Cb_global_temp = 0.0_dbl
Cb_numFluids_temp = 0_lng
CALL MPI_REDUCE(Cb_global,Cb_global_temp,1,MPI_DOUBLE_PRECISION,MPI_SUM,master,MPI_COMM_WORLD,mpierr)
CALL MPI_REDUCE(Cb_numFluids,Cb_numFluids_temp,1,MPI_INTEGER,MPI_SUM,master,MPI_COMM_WORLD,mpierr)
Cb_global = Cb_global_temp/Cb_numFluids_temp
Cb_numFLuids = Cb_numFluids_temp
!write(*,*) 'In Collect and Distribute before bcas',iter,mySub,Cb_global,Cb_global_temp,Cb_numFluids,Cb_numFluids_temp
CALL MPI_BCast(Cb_global,1,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,mpierr)
CALL MPI_BCast(Cb_numFluids,1,MPI_INTEGER,master,MPI_COMM_WORLD,mpierr)
!CALL MPI_BARRIER(MPI_COMM_WORLD,mpierr)					! synchronize all processing units before next loop [Intrinsic]
!write(*,*) 'In Collect and Distribute after Bcast',iter,mySub,Cb_global,Cb_numFluids,Cb_numFluids*xcf*ycf*zcf,Cb_numFluids_temp

!--------------------------------------------------------------------------------------------------
END SUBROUTINE  Collect_Distribute_Global_Bulk_Scalar_Conc
!--------------------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------------------
SUBROUTINE Particle_MPI_Transfer		! transfer the particles to neighbouring partitions
!--------------------------------------------------------------------------------------------------
IMPLICIT NONE
TYPE(ParRecord), POINTER :: current  => NULL()
TYPE(ParRecord), POINTER :: next => NULL()
INTEGER(lng) :: nreqs,source,sendtag,dest,recvtag		! number of send/recv requests
INTEGER(lng) :: mpierr,commdir,iComm,par_num!,numparticlesDomain			! MPI standard error object
LOGICAL :: probeflag
INTEGER(lng) :: parsendindex(NumCommDirsPar), parrecvindex(NumCommDirsPar)

!write(*,*) 'In Particle_MPI	','FLAG=',ParticleTransfer,'SUBID = ',mySub,'iter = ',iter

! IMPORTANT NOTE: It is important for the partransfersend and partransferrecv arrays to store the particle data in columns, i.e. of 
! dimenions NumParVar x NumCommDirsPar as against NumCommDirsPar x  NumParVar as we are planning to send packets of data 
! corresponding to an entire particle across processors. In fortran and MPI, this requires that the entire particle data 
!(all NumParVar  variables given by partransfersend(:,iComm)) need to be in column format so that they can be passed as a single array.
! For some yet-to-be understood reason, MPI dopes not accept passing an entire row of data, say, partransfersend(iComm,:) 
!as a contiguous array into MPI_SEND or MPI_RECV. Hence the convetion used below. 

nreqs = 0_lng
partransfersend = 0.0_dbl
partransferrecv = 0.0_dbl
sendtag = 0_lng



!Identify how many particles need to be transfererd in each direction
numpartransfer = 0_lng
current => ParListHead%next
DO WHILE (ASSOCIATED(current))
	next => current%next ! copy pointer of next node
	IF (current%pardata%cur_part /= current%pardata%new_part) THEN ! Transfer only if the particle has moved to another partition
		DO iComm = 1,NumCommDirsPar
			IF (SubID(iComm).EQ.current%pardata%new_part) THEN
				numpartransfer(iComm) = numpartransfer(iComm)+1_lng
			END IF
		END DO
	END IF
	! point to next node in the list
	current => next
	!write(*,*) current%parid
ENDDO
! Now 'numpartransfer' contains the number of particles in each direction


!numparticlesDomain = 1_lng!MAXVAL(numpartransfer)
ALLOCATE(ParSendArray(numparticlesDomain,NumCommDirsPar))
ALLOCATE(ParRecvArray(numparticlesDomain,NumCommDirsPar))
!ALLOCATE(ParSendArray(1,NumCommDirsPar))
!ALLOCATE(ParRecvArray(1,NumCommDirsPar))
DO iComm = 1,NumCommDirsPar
	DO par_num = 1,numparticlesDomain
		ParSendArray(par_num,iComm) = ParInit
		ParRecvArray(par_num,iComm) = ParInit
	END DO
END DO

!write(*,*) 'In MPI_Particle_transfer - allocated particle send and recv arrays',mySub,iter
!CALL MPI_BARRIER(MPI_COMM_WORLD,mpierr)					! synchronize all processing units before next loop [Intrinsic]

! Pack data
parsendindex = 0_lng
current => ParListHead%next
DO WHILE (ASSOCIATED(current))
	next => current%next ! copy pointer of next node
	IF (current%pardata%cur_part /= current%pardata%new_part) THEN ! Transfer only if the particle has moved to another partition
		DO iComm = 1,NumCommDirsPar
			IF (SubID(iComm).EQ.current%pardata%new_part) THEN
				commdir = iComm
				parsendindex(commdir) = parsendindex(commdir) + 1_lng
				!write(*,*) 'data being packed in ',mySub,'for direction ',commdir
			END IF
		END DO
		! PACK DATA INTO AN ARRAY

!		ParSendArray(parsendindex(commdir),commdir)%parid 	= current%pardata%parid
!		ParSendArray(parsendindex(commdir),commdir)%cur_part 	= current%pardata%cur_part
!		ParSendArray(parsendindex(commdir),commdir)%new_part 	= current%pardata%new_part
!		ParSendArray(parsendindex(commdir),commdir)%xp 		= current%pardata%xp
!		ParSendArray(parsendindex(commdir),commdir)%yp 		= current%pardata%yp
!		ParSendArray(parsendindex(commdir),commdir)%zp 		= current%pardata%zp
!		ParSendArray(parsendindex(commdir),commdir)%up 		= current%pardata%up
!		ParSendArray(parsendindex(commdir),commdir)%vp 		= current%pardata%vp
!		ParSendArray(parsendindex(commdir),commdir)%wp 		= current%pardata%wp
!		ParSendArray(parsendindex(commdir),commdir)%rp 		= current%pardata%rp
!		ParSendArray(parsendindex(commdir),commdir)%delNBbyCV 	= current%pardata%delNBbyCV
!		ParSendArray(parsendindex(commdir),commdir)%par_conc 	= current%pardata%par_conc
!		ParSendArray(parsendindex(commdir),commdir)%bulk_conc 	= current%pardata%bulk_conc
!		ParSendArray(parsendindex(commdir),commdir)%xpold 	= current%pardata%xpold
!		ParSendArray(parsendindex(commdir),commdir)%ypold 	= current%pardata%ypold
!		ParSendArray(parsendindex(commdir),commdir)%zpold 	= current%pardata%zpold
!		ParSendArray(parsendindex(commdir),commdir)%upold 	= current%pardata%upold
!		ParSendArray(parsendindex(commdir),commdir)%vpold 	= current%pardata%vpold
!		ParSendArray(parsendindex(commdir),commdir)%wpold 	= current%pardata%wpold
!		ParSendArray(parsendindex(commdir),commdir)%rpold 	= current%pardata%rpold
!		ParSendArray(parsendindex(commdir),commdir)%sh 		= current%pardata%sh
!		ParSendArray(parsendindex(commdir),commdir)%gamma_cont 	= current%pardata%gamma_cont
!		!write(*,*) 'data being packed in ',mySub,'for direction ',commdir

		ParSendArray(parsendindex(commdir),commdir)		= current%pardata




!		!sendtag = current%parid
!		partransfersend(1,commdir) = REAL(current%pardata%parid,dbl)
!		partransfersend(2,commdir) = REAL(current%pardata%xp,dbl)
!		partransfersend(3,commdir) = REAL(current%pardata%yp,dbl)
!		partransfersend(4,commdir) = REAL(current%pardata%zp,dbl)
!		partransfersend(5,commdir) = REAL(current%pardata%up,dbl)
!		partransfersend(6,commdir) = REAL(current%pardata%vp,dbl)
!		partransfersend(7,commdir) = REAL(current%pardata%wp,dbl)
!		partransfersend(8,commdir) = REAL(current%pardata%rp,dbl)
!		partransfersend(9,commdir) = REAL(current%pardata%delNBbyCV,dbl)
!		partransfersend(10,commdir) = REAL(current%pardata%par_conc,dbl)
!		partransfersend(11,commdir) = REAL(current%pardata%bulk_conc,dbl)
!		partransfersend(12,commdir) = REAL(current%pardata%rpold,dbl)
!		partransfersend(13,commdir) = REAL(current%pardata%sh,dbl)
!		partransfersend(14,commdir) = REAL(current%pardata%gamma_cont,dbl)
!		partransfersend(15,commdir) = REAL(current%pardata%cur_part,dbl)
!		partransfersend(16,commdir) = REAL(current%pardata%new_part,dbl)
!		!write(*,*) current%pardata%parid,current%pardata%cur_part,current%pardata%new_part
!		!write(*,*) 'data being packed in ',mySub,'for direction ',commdir
	END IF
	! point to next node in the list
	current => next
	!write(*,*) current%parid
ENDDO

!write(*,*) 'In MPI_Particle_transfer - packed data into send array',mySub,iter
!CALL MPI_BARRIER(MPI_COMM_WORLD,mpierr)					! synchronize all processing units before next loop [Intrinsic]

!! Post receives
!parrecvindex = 0_lng
!DO iComm = 1,NumCommDirsPar
!	!write(*,*) 'In post receive ', mySub,iComm,SubID(OppCommDir(iComm))
!	IF(SubID(OppCommDir(iComm)) .NE. 0) THEN
!		nreqs 	= nreqs + 1_lng
!		source	= SubID(OppCommDir(iComm)) - 1_lng ! rank of processing unit sending message TO this processing unit
!		recvtag = iComm + 20_lng
!10		CALL MPI_IPROBE(source,recvtag,MPI_COMM_WORLD,probeflag,probestat,mpierr)
!		CALL MPI_GET_COUNT(probestat,mpipartransfertype,parrecvindex(iComm),mpierr)	
!		IF (parrecvindex(iComm).LT.1_lng) THEN
!			GOTO 10
!		END IF
!
!		!!CALL MPI_IPROBE(MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,probeflag,probestat,mpierr)	
!		!!nreqs = nreqs + 1_lng
!		!!recvtag = probestat(MPI_TAG)
!		!!source =  probestat(MPI_SOURCE)
!!10		CALL MPI_IPROBE(source,recvtag,MPI_COMM_WORLD,probeflag,probestat,mpierr)
!!		CALL MPI_GET_COUNT(probestat,mpipartransfertype,parrecvindex(iComm),mpierr)	
!!		IF (parrecvindex(iComm).LT.1_lng) THEN
!!			GOTO 10
!!		END IF
!
!		!write(*,*) 'what we get from direction iComm ',parrecvindex(iComm),iComm,mySub,source,iter
!
!		!CALL MPI_IRECV(recvarray,NumParVar,MPI_DOUBLE_PRECISION,source,recvtag,MPI_COMM_WORLD,parreqid(nreqs),mpierr)
!		!CALL MPI_IRECV(partransferrecv(1:NumParVar,iComm),NumParVar,MPI_DOUBLE_PRECISION,source,recvtag,MPI_COMM_WORLD,parreqid(nreqs),mpierr)
!		!CALL MPI_IRECV(partransferrecv(:,iComm),NumParVar,MPI_DOUBLE_PRECISION,source,recvtag,MPI_COMM_WORLD,parreqid(nreqs),mpierr)
!		!CALL MPI_IRECV(ParRecvArray(1,iComm),1,mpipartransfertype,source,recvtag,MPI_COMM_WORLD,parreqid(nreqs),mpierr)
!		CALL MPI_IRECV(ParRecvArray(1:parrecvindex(iComm),iComm),parrecvindex(iComm),mpipartransfertype,source,recvtag,MPI_COMM_WORLD,parreqid(nreqs),mpierr)
!		!partransferrecv(iComm,:) = recvarray
!		write(*,*) '********recd  ','with tag ',recvtag,'from ',source+1_lng,'to ',mySub,'at iter ',iter,parrecvindex(iComm),' items'
!		!write(*,*) iComm,partransferrecv(iComm,1),recvarray(1)
!	END IF
!END DO

! Post send
DO iComm = 1,NumCommDirsPar
	!write(*,*) 'In post send ',mySub,iComm,SubID(iComm)
	IF(SubID(iComm) .NE. 0) THEN
		nreqs 	= nreqs + 1_lng
		dest	= SubID(iComm) - 1_lng ! rank of processing unit sending message TO this processing unit
		sendtag = iComm + 20_lng
		!sendtag = mySub*1000_lng + iComm*100_lng + 20_lng
		!sendarray = partransfersend(iComm,:)
		!CALL MPI_ISEND(sendarray,NumParVar,MPI_DOUBLE_PRECISION,dest,sendtag,MPI_COMM_WORLD,parreqid(nreqs),mpierr)	! send data
		!CALL MPI_ISEND(partransfersend(1:NumParVar,iComm),NumParVar,MPI_DOUBLE_PRECISION,dest,sendtag,MPI_COMM_WORLD,parreqid(nreqs),mpierr)	! send data
		!CALL MPI_ISEND(partransfersend(:,iComm),NumParVar,MPI_DOUBLE_PRECISION,dest,sendtag,MPI_COMM_WORLD,parreqid(nreqs),mpierr)	! send data
		!CALL MPI_ISEND(ParSendArray(1,iComm),1,mpipartransfertype,dest,sendtag,MPI_COMM_WORLD,parreqid(nreqs),mpierr)	! send data
		CALL MPI_ISEND(ParSendArray(1:max(parsendindex(iComm),1),iComm),max(parsendindex(iComm),1),mpipartransfertype,dest,sendtag,MPI_COMM_WORLD,parreqid(nreqs),mpierr)	! send data
		!write(*,*) '**********sent   ','with tag ',sendtag,'from ',mySub,'to ',dest+1_lng,'at iter ',iter,parsendindex(iComm),' items'
		!write(*,*) iComm,partransfersend(iComm,1),sendarray(1)
	END IF
END DO

!write(*,*) 'In MPI_Particle_transfer - post send',mySub,iter
!CALL MPI_BARRIER(MPI_COMM_WORLD,mpierr)					! synchronize all processing units before next loop [Intrinsic]



! Post receives
parrecvindex = 0_lng
DO iComm = 1,NumCommDirsPar
	!write(*,*) 'In post receive ', mySub,iComm,SubID(OppCommDir(iComm))
	IF(SubID(OppCommDir(iComm)) .NE. 0) THEN
		nreqs 	= nreqs + 1_lng
		source	= SubID(OppCommDir(iComm)) - 1_lng ! rank of processing unit sending message TO this processing unit
		recvtag = iComm + 20_lng
		!recvtag = source*1000_lng + iComm*1000_lng + 20_lng
!10		CALL MPI_IPROBE(source,recvtag,MPI_COMM_WORLD,probeflag,probestat,mpierr)
10		CALL MPI_PROBE(source,recvtag,MPI_COMM_WORLD,probestat,mpierr)
		CALL MPI_GET_COUNT(probestat,mpipartransfertype,parrecvindex(iComm),mpierr)	
		IF (parrecvindex(iComm).LT.1_lng) THEN
			GOTO 10
		END IF

		!!CALL MPI_IPROBE(MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,probeflag,probestat,mpierr)	
		!!nreqs = nreqs + 1_lng
		!!recvtag = probestat(MPI_TAG)
		!!source =  probestat(MPI_SOURCE)
!10		CALL MPI_IPROBE(source,recvtag,MPI_COMM_WORLD,probeflag,probestat,mpierr)
!		CALL MPI_GET_COUNT(probestat,mpipartransfertype,parrecvindex(iComm),mpierr)	
!		IF (parrecvindex(iComm).LT.1_lng) THEN
!			GOTO 10
!		END IF

		!write(*,*) 'what we get from direction iComm ',parrecvindex(iComm),iComm,mySub,source,iter

		!CALL MPI_IRECV(recvarray,NumParVar,MPI_DOUBLE_PRECISION,source,recvtag,MPI_COMM_WORLD,parreqid(nreqs),mpierr)
		!CALL MPI_IRECV(partransferrecv(1:NumParVar,iComm),NumParVar,MPI_DOUBLE_PRECISION,source,recvtag,MPI_COMM_WORLD,parreqid(nreqs),mpierr)
		!CALL MPI_IRECV(partransferrecv(:,iComm),NumParVar,MPI_DOUBLE_PRECISION,source,recvtag,MPI_COMM_WORLD,parreqid(nreqs),mpierr)
		!CALL MPI_IRECV(ParRecvArray(1,iComm),1,mpipartransfertype,source,recvtag,MPI_COMM_WORLD,parreqid(nreqs),mpierr)
		CALL MPI_IRECV(ParRecvArray(1:parrecvindex(iComm),iComm),parrecvindex(iComm),mpipartransfertype,source,recvtag,MPI_COMM_WORLD,parreqid(nreqs),mpierr)
		!partransferrecv(iComm,:) = recvarray
		!write(*,*) '********recd  ','with tag ',recvtag,'from ',source+1_lng,'to ',mySub,'at iter ',iter,parrecvindex(iComm),' items'
		!write(*,*) iComm,partransferrecv(iComm,1),recvarray(1)
	END IF
END DO

!write(*,*) 'In MPI_Particle_transfer - post receives',mySub,iter
!CALL MPI_BARRIER(MPI_COMM_WORLD,mpierr)					! synchronize all processing units before next loop [Intrinsic]

! Wait for all communication to complete
CALL MPI_WAITALL(nreqs,parreqid(1:nreqs),parwtstat,mpierr)


! Unpack data if recd.
DO iComm = 1,NumCommDirsPar
	DO par_num = 1,parrecvindex(iComm)
		!IF (partransferrecv(1,iComm).GT.0_lng) THEN
		IF (ParRecvArray(par_num,iComm)%parid.GT.0_lng) THEN
			CALL list_insert(ParListHead)
			current => ParListHead%next
	
			! UNPACK DATA
 
!			current%pardata%parid 		=ParRecvArray(par_num,iComm)%parid
!			current%pardata%cur_part 	=ParRecvArray(par_num,iComm)%cur_part
!			current%pardata%new_part 	=ParRecvArray(par_num,iComm)%new_part
!			current%pardata%xp 		=ParRecvArray(par_num,iComm)%xp
!			current%pardata%yp 		=ParRecvArray(par_num,iComm)%yp
!			current%pardata%zp 		=ParRecvArray(par_num,iComm)%zp
!			current%pardata%up 		=ParRecvArray(par_num,iComm)%up
!			current%pardata%vp 		=ParRecvArray(par_num,iComm)%vp
!			current%pardata%wp 		=ParRecvArray(par_num,iComm)%wp
!			current%pardata%rp 		=ParRecvArray(par_num,iComm)%rp
!			current%pardata%delNBbyCV 	=ParRecvArray(par_num,iComm)%delNBbyCV
!			current%pardata%par_conc 	=ParRecvArray(par_num,iComm)%par_conc
!			current%pardata%bulk_conc 	=ParRecvArray(par_num,iComm)%bulk_conc
!			current%pardata%xpold 		=ParRecvArray(par_num,iComm)%xpold
!			current%pardata%ypold 		=ParRecvArray(par_num,iComm)%ypold
!			current%pardata%zpold 		=ParRecvArray(par_num,iComm)%zpold
!			current%pardata%upold 		=ParRecvArray(par_num,iComm)%upold
!			current%pardata%vpold 		=ParRecvArray(par_num,iComm)%vpold
!			current%pardata%wpold 		=ParRecvArray(par_num,iComm)%wpold
!			current%pardata%rpold 		=ParRecvArray(par_num,iComm)%rpold
!			current%pardata%sh 		=ParRecvArray(par_num,iComm)%sh
!			current%pardata%gamma_cont 	=ParRecvArray(par_num,iComm)%gamma_cont
!			!write(*,*) 'unpacking data ',mySub,current%pardata%parid,iComm,current%pardata%cur_part,current%pardata%new_part
!			!PAUSE
	
			current%pardata 		=ParRecvArray(par_num,iComm)

	
!			current%parid 		= INT(partransferrecv(1,iComm),lng)
!			!write(*,*) 'unpacking data ',mySub,current%pardata%parid
!			!PAUSE
!			current%pardata%xp 		= partransferrecv(2,iComm)
!			current%pardata%yp 		= partransferrecv(3,iComm)
!			current%pardata%zp		= partransferrecv(4,iComm)
!			current%pardata%up 		= partransferrecv(5,iComm)
!			current%pardata%vp 		= partransferrecv(6,iComm)
!			current%pardata%wp 		= partransferrecv(7,iComm)
!			current%pardata%rp 		= partransferrecv(8,iComm)
!			current%pardata%delNBbyCV 	= partransferrecv(9,iComm)
!			current%pardata%par_conc 	= partransferrecv(10,iComm)
!			current%pardata%bulk_conc 	= partransferrecv(11,iComm)
!			current%pardata%rpold 		= partransferrecv(12,iComm)
!			current%pardata%sh 		= partransferrecv(13,iComm)
!			current%pardata%gamma_cont 	= partransferrecv(14,iComm)
!			current%pardata%cur_part 	= INT(partransferrecv(15,iComm),lng)
!			current%pardata%new_part 	= INT(partransferrecv(16,iComm),lng)
	
			current%pardata%cur_part 	= mySub!current%pardata%new_part
		END IF
	END DO
END DO

! If data was sent then remove sent particle from list
current => ParListHead%next
DO WHILE (ASSOCIATED(current))
	next => current%next ! copy pointer of next node
	!write(*,*) 'traverse node',mySub,current%pardata%parid,current%pardata%cur_part,current%pardata%new_part,iter
	!write(*,*) current%prev%next%pardata%parid,current%pardata%parid,iter
	IF (current%pardata%cur_part .NE. current%pardata%new_part) THEN ! Transfer only if the particle has moved to another partition
		!write(*,*) 'In delete node',mySub,current%pardata%parid,current%pardata%cur_part,current%pardata%new_part,iter
		!PAUSE
		CALL list_delete(current)

		!current%prev%next => next
		!IF (ASSOCIATED(next)) THEN
		!	write(*,*) 'current%next is associated'
		!	next%prev => current%prev
		!ENDIF
		!NULLIFY(current%prev)
		!NULLIFY(current%next)
		!DEALLOCATE(current)

		NULLIFY(current)
	END IF
	! point to next node in the list
	current => next
	!IF (ASSOCIATED(current)) THEN
	!	write(*,*) 'test association prev ',ASSOCIATED(current%prev),iter
	!	write(*,*) 'test association next ',ASSOCIATED(current%next),iter
	!END IF
	!write(*,*) current%pardata%parid
ENDDO
!write(*,*) mySub,ParListHead%next%pardata%parid
numparticlesSub = 0_lng
current => ParListHead%next
DO WHILE (ASSOCIATED(current))
	numparticlesSub = numparticlesSub + 1_lng
	next => current%next ! copy pointer of next node
	!write(*,*) 'leftover node',mySub,current%pardata%parid,current%pardata%cur_part,current%pardata%new_part
	! point to next node in the list
	current => next
	!write(*,*) current%pardata%parid
ENDDO
!write(*,*) 'In MPIParticleTransfer  ',iter,mySub,numparticlesSub

DEALLOCATE(ParSendArray)
DEALLOCATE(ParRecvArray)

!write(*,*) mySub,'____________________________________________________________________'



! Rough outline of the algorithm
! PACK DATA INTO OPTIMAL NUMBER OF PACKETS
! POST RECEIVE
! SEND DATA
! MPI_WAITALL - make sure all the data transfer is completed
! UNPACK DATA AND ADD DATA TO THE PARTICLE LINKED LIST

!------------------------------------------------
END SUBROUTINE Particle_MPI_Transfer
!------------------------------------------------
!! Working for single particle
!!--------------------------------------------------------------------------------------------------
!SUBROUTINE Particle_MPI_Transfer		! transfer the particles to neighbouring partitions
!!--------------------------------------------------------------------------------------------------
!IMPLICIT NONE
!TYPE(ParRecord), POINTER :: current
!TYPE(ParRecord), POINTER :: next
!INTEGER(lng) :: nreqs,source,sendtag,dest,recvtag		! number of send/recv requests
!INTEGER(lng) :: mpierr,commdir,iComm			! MPI standard error object
!LOGICAL :: probeflag
!REAL(dbl) :: sendarray(NumParVar),recvarray(NumParVar)
!
!
!nreqs = 0_lng
!partransfersend = 0.0_dbl
!partransferrecv = 0.0_dbl
!sendarray = 0.0_dbl
!recvarray = 0.0_dbl
!sendtag = 0_lng
!!write(*,*) 'In Particle_MPI	','FLAG=',ParticleTransfer,'SUBID = ',mySub,'iter = ',iter
!
!! IMPORTANT NOTE: It is important for the partransfersend and partransferrecv arrays to store the particle data in columns, i.e. of 
!! dimenions NumParVar x NumCommDirsPar as against NumCommDirsPar x  NumParVar as we are planning to send packets of data 
!! corresponding to an entire particle across processors. In fortran and MPI, this requires that the entire particle data 
!!(all NumParVar  variables given by partransfersend(:,iComm)) need to be in column format so that they can be passed as a single array.
!! For some yet-to-be understood reason, MPI dopes not accept passing an entire row of data, say, partransfersend(iComm,:) 
!!as a contiguous array into MPI_SEND or MPI_RECV. Hence the convetion used below. 
!
!
!! Pack data
!current => ParListHead%next
!DO WHILE (ASSOCIATED(current))
!	next => current%next ! copy pointer of next node
!	IF (current%cur_part /= current%new_part) THEN ! Transfer only if the particle has moved to another partition
!		DO iComm = 1,NumCommDirsPar
!			IF (SubID(iComm).EQ.current%new_part) THEN
!				commdir = iComm	
!			END IF
!		END DO
!		! PACK DATA INTO AN ARRAY
!		!sendtag = current%parid
!		partransfersend(1,commdir) = REAL(current%parid,dbl)
!		partransfersend(2,commdir) = REAL(current%xp,dbl)
!		partransfersend(3,commdir) = REAL(current%yp,dbl)
!		partransfersend(4,commdir) = REAL(current%zp,dbl)
!		partransfersend(5,commdir) = REAL(current%up,dbl)
!		partransfersend(6,commdir) = REAL(current%vp,dbl)
!		partransfersend(7,commdir) = REAL(current%wp,dbl)
!		partransfersend(8,commdir) = REAL(current%rp,dbl)
!		partransfersend(9,commdir) = REAL(current%delNBbyCV,dbl)
!		partransfersend(10,commdir) = REAL(current%par_conc,dbl)
!		partransfersend(11,commdir) = REAL(current%bulk_conc,dbl)
!		partransfersend(12,commdir) = REAL(current%rpold,dbl)
!		partransfersend(13,commdir) = REAL(current%sh,dbl)
!		partransfersend(14,commdir) = REAL(current%gamma_cont,dbl)
!		partransfersend(15,commdir) = REAL(current%cur_part,dbl)
!		partransfersend(16,commdir) = REAL(current%new_part,dbl)
!		!write(*,*) current%parid,current%cur_part,current%new_part
!		!write(*,*) 'data being packed in ',mySub,'for direction ',commdir
!	END IF
!	! point to next node in the list
!	current => next
!	!write(*,*) current%parid
!ENDDO
!
!! Post receives
!DO iComm = 1,NumCommDirsPar
!	!write(*,*) 'In post receive ', mySub,iComm,SubID(OppCommDir(iComm))
!	IF(SubID(OppCommDir(iComm)) .NE. 0) THEN
!		nreqs 	= nreqs + 1_lng
!		source	= SubID(OppCommDir(iComm)) - 1_lng ! rank of processing unit sending message TO this processing unit
!		recvtag = iComm + 20_lng
!		!CALL MPI_IRECV(recvarray,NumParVar,MPI_DOUBLE_PRECISION,source,recvtag,MPI_COMM_WORLD,parreqid(nreqs),mpierr)
!		!CALL MPI_IRECV(partransferrecv(1:NumParVar,iComm),NumParVar,MPI_DOUBLE_PRECISION,source,recvtag,MPI_COMM_WORLD,parreqid(nreqs),mpierr)
!		CALL MPI_IRECV(partransferrecv(:,iComm),NumParVar,MPI_DOUBLE_PRECISION,source,recvtag,MPI_COMM_WORLD,parreqid(nreqs),mpierr)
!		!partransferrecv(iComm,:) = recvarray
!		!write(*,*) '********recd  ','with tag ',recvtag,'from ',source+1_lng,'to ',mySub,'at iter ',iter
!		!write(*,*) iComm,partransferrecv(iComm,1),recvarray(1)
!	END IF
!END DO
!! Post send
!DO iComm = 1,NumCommDirsPar
!	!write(*,*) 'In post send ',mySub,iComm,SubID(iComm)
!	IF(SubID(iComm) .NE. 0) THEN
!		nreqs 	= nreqs + 1_lng
!		dest	= SubID(iComm) - 1_lng ! rank of processing unit sending message TO this processing unit
!		sendtag = iComm + 20_lng
!		!sendarray = partransfersend(iComm,:)
!		!CALL MPI_ISEND(sendarray,NumParVar,MPI_DOUBLE_PRECISION,dest,sendtag,MPI_COMM_WORLD,parreqid(nreqs),mpierr)	! send data
!		!CALL MPI_ISEND(partransfersend(1:NumParVar,iComm),NumParVar,MPI_DOUBLE_PRECISION,dest,sendtag,MPI_COMM_WORLD,parreqid(nreqs),mpierr)	! send data
!		CALL MPI_ISEND(partransfersend(:,iComm),NumParVar,MPI_DOUBLE_PRECISION,dest,sendtag,MPI_COMM_WORLD,parreqid(nreqs),mpierr)	! send data
!		!write(*,*) '**********sent   ','with tag ',sendtag,'from ',mySub,'to ',dest+1_lng,'at iter ',iter
!		!write(*,*) iComm,partransfersend(iComm,1),sendarray(1)
!	END IF
!END DO
!
!! Wait for all communication to complete
!CALL MPI_WAITALL(nreqs,parreqid(1:nreqs),parwtstat,mpierr)
!
!
!! Unpack data if recd.
!DO iComm = 1,NumCommDirsPar
!	IF (partransferrecv(1,iComm).GT.0_lng) THEN
!		CALL list_insert(ParListHead)
!		current => ParListHead%next
!		! PACK DATA INTO AN ARRAY
!		current%parid 		= INT(partransferrecv(1,iComm),lng)
!		!write(*,*) 'unpacking data ',mySub,current%parid
!		!PAUSE
!		current%xp 		= partransferrecv(2,iComm)
!		current%yp 		= partransferrecv(3,iComm)
!		current%zp		= partransferrecv(4,iComm)
!		current%up 		= partransferrecv(5,iComm)
!		current%vp 		= partransferrecv(6,iComm)
!		current%wp 		= partransferrecv(7,iComm)
!		current%rp 		= partransferrecv(8,iComm)
!		current%delNBbyCV 	= partransferrecv(9,iComm)
!		current%par_conc 	= partransferrecv(10,iComm)
!		current%bulk_conc 	= partransferrecv(11,iComm)
!		current%rpold 		= partransferrecv(12,iComm)
!		current%sh 		= partransferrecv(13,iComm)
!		current%gamma_cont 	= partransferrecv(14,iComm)
!		current%cur_part 	= INT(partransferrecv(15,iComm),lng)
!		current%new_part 	= INT(partransferrecv(16,iComm),lng)
!	
!		current%cur_part 	= current%new_part
!	END IF
!END DO
!
!!!!!! OLD PIECE OF CODE
!!nreqs = nreqs + 1_lng
!!IF (myid.EQ.1_lng) THEN
!! source = 0_lng
!! recvtag = 21_lng 
!! iComm = 1_lng
!!ELSE
!! source = 1_lng
!! recvtag = 12_lng 
!! iComm = 2_lng
!!END IF
!!CALL MPI_IRECV(recvarray,16_lng,MPI_DOUBLE_PRECISION,source,recvtag,MPI_COMM_WORLD,parreqid(nreqs),mpierr)
!!write(*,*) '********recd  ','parid ',recvtag,'from ',source+1_lng,'to ',mySub,'at iter ',iter
!!!partransferrecv(iComm,:) = recvarray
!!
!!! Post send
!!nreqs = nreqs + 1_lng
!!IF (myid.EQ.1_lng) THEN
!! dest = 0_lng 
!! sendtag = 12_lng
!! iComm = 2_lng
!!ELSE
!! dest = 1_lng
!! sendtag = 21_lng 
!! iComm = 1_lng
!!END IF
!!sendarray = partransfersend(iComm,:)
!!CALL MPI_ISEND(sendarray,16_lng,MPI_DOUBLE_PRECISION,dest,sendtag,MPI_COMM_WORLD,parreqid(nreqs),mpierr)	! send data
!!write(*,*) '**********sent   ','parid ',sendtag,'from ',mySub,'to ',dest+1_lng,'at iter ',iter
!!
!!
!!! Wait for all communication to complete
!!CALL MPI_WAITALL(nreqs,parreqid(1:nreqs),parwtstat,mpierr)
!!
!!! Unpack data if recd. 
!!IF (recvarray(1).GT.0_lng) THEN
!!	CALL list_insert(ParListHead)
!!	current => ParListHead%next
!!	! PACK DATA INTO AN ARRAY
!!	current%parid 		= INT(recvarray(1),lng)
!!	!write(*,*) mySub,current%parid,recvarray(1),recvarray(5),recvarray(6)
!!	!PAUSE
!!	current%xp 		= recvarray(2)
!!	current%yp 		= recvarray(3)
!!	current%zp		= recvarray(4)
!!	current%up 		= recvarray(5)
!!	current%vp 		= recvarray(6)
!!	current%wp 		= recvarray(7)
!!	current%rp 		= recvarray(8)
!!	current%delNBbyCV 	= recvarray(9)
!!	current%par_conc 	= recvarray(10)
!!	current%bulk_conc 	= recvarray(11)
!!	current%rpold 		= recvarray(12)
!!	current%sh 		= recvarray(13)
!!	current%gamma_cont 	= recvarray(14)
!!	current%cur_part 	= INT(recvarray(15),lng)
!!	current%new_part 	= INT(recvarray(16),lng)
!!
!!	current%cur_part 	= current%new_part
!!END IF
!
!! If data was sent then remove sent particle from list
!current => ParListHead%next
!DO WHILE (ASSOCIATED(current))
!	next => current%next ! copy pointer of next node
!	IF (current%cur_part /= current%new_part) THEN ! Transfer only if the particle has moved to another partition
!	!write(*,*) current%parid,current%cur_part,current%new_part
!	CALL list_delete(current)
!	END IF
!	! point to next node in the list
!	current => next
!	!write(*,*) current%parid
!ENDDO
!
!write(*,*) '____________________________________________________________________'
!
!
!
!! Rough outline of the algorithm
!! PACK DATA INTO OPTIMAL NUMBER OF PACKETS
!! POST RECEIVE
!! SEND DATA
!! MPI_WAITALL - make sure all the data transfer is completed
!! UNPACK DATA AND ADD DATA TO THE PARTICLE LINKED LIST
!
!!------------------------------------------------
!END SUBROUTINE Particle_MPI_Transfer
!!------------------------------------------------
!!--------------------------------------------------------------------------------------------------
!SUBROUTINE Particle_MPI_Transfer	! transfer the particles to neighbouring partitions
!!--------------------------------------------------------------------------------------------------
!IMPLICIT NONE
!TYPE(ParRecord), POINTER :: current
!TYPE(ParRecord), POINTER :: next
!INTEGER(lng) :: nreqs,source,sendtag,dest,recvtag		! number of send/recv requests
!INTEGER(lng) :: mpierr			! MPI standard error object
!LOGICAL :: probeflag
!
!
!nreqs = 0_lng
!partransfersend = 0.0_dbl
!partransferrecv = 0.0_dbl
!sendtag = 0_lng
!write(*,*) 'In Particle_MPI	','FLAG=',ParticleTransfer,'SUBID = ',mySub,'iter = ',iter
!
!! Pack data
!current => ParListHead%next
!DO WHILE (ASSOCIATED(current))
!	next => current%next ! copy pointer of next node
!	IF (current%cur_part /= current%new_part) THEN ! Transfer only if the particle has moved to another partition
!		! PACK DATA INTO AN ARRAY
!		sendtag = current%parid
!		partransfersend(1) = REAL(current%parid,dbl)
!		!write(*,*) mySub,current%parid,partransfersend(1)
!		partransfersend(2) = REAL(current%xp,dbl)
!		partransfersend(3) = REAL(current%yp,dbl)
!		partransfersend(4) = REAL(current%zp,dbl)
!		partransfersend(5) = REAL(current%up,dbl)
!		partransfersend(6) = REAL(current%vp,dbl)
!		partransfersend(7) = REAL(current%wp,dbl)
!		partransfersend(8) = REAL(current%rp,dbl)
!		partransfersend(9) = REAL(current%delNBbyCV,dbl)
!		partransfersend(10) = REAL(current%par_conc,dbl)
!		partransfersend(11) = REAL(current%bulk_conc,dbl)
!		partransfersend(12) = REAL(current%rpold,dbl)
!		partransfersend(13) = REAL(current%sh,dbl)
!		partransfersend(14) = REAL(current%gamma_cont,dbl)
!		partransfersend(15) = REAL(current%cur_part,dbl)
!		partransfersend(16) = REAL(current%new_part,dbl)
!		!write(*,*) current%parid,current%cur_part,current%new_part
!	END IF
!	! point to next node in the list
!	current => next
!	!write(*,*) current%parid
!ENDDO
!
!! Post receives
!!CALL MPI_IPROBE(MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,probeflag,probestat,mpierr)	
!!IF (probeflag) THEN
!!	nreqs = nreqs + 1_lng  
!!	!nreqs = 2_lng  
!!	!CALL MPI_IRECV(partransferrecv,16_lng,MPI_DOUBLE_PRECISION,probestat(MPI_SOURCE),probestat(MPI_TAG),MPI_COMM_WORLD,parreqid(nreqs),mpierr)
!!	CALL MPI_IRECV(partransferrecv,16_lng,MPI_DOUBLE_PRECISION,probestat(MPI_SOURCE),probestat(MPI_TAG),MPI_COMM_WORLD,parreqid(nreqs),mpierr)
!!END IF
!
!!CALL MPI_IPROBE(MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,probeflag,probestat,mpierr)	
!!nreqs = nreqs + 1_lng
!!recvtag = probestat(MPI_TAG)
!!source =  probestat(MPI_SOURCE)
!nreqs = nreqs + 1_lng
!IF (myid.EQ.1_lng) THEN
! source = 0_lng
! recvtag = 21_lng 
!ELSE
! source = 1_lng
! recvtag = 12_lng 
!END IF
!CALL MPI_IRECV(partransferrecv,16_lng,MPI_DOUBLE_PRECISION,source,recvtag,MPI_COMM_WORLD,parreqid(nreqs),mpierr)
!write(*,*) '********recd  ','parid ',recvtag,'from ',source+1_lng,'to ',mySub,'at iter ',iter
!
!! Post send
!nreqs = nreqs + 1_lng
!IF (myid.EQ.1_lng) THEN
! dest = 0_lng 
! sendtag = 12_lng 
!ELSE
! dest = 1_lng
! sendtag = 21_lng 
!END IF
!CALL MPI_ISEND(partransfersend,16_lng,MPI_DOUBLE_PRECISION,dest,sendtag,MPI_COMM_WORLD,parreqid(nreqs),mpierr)	! send data
!write(*,*) '**********sent   ','parid ',sendtag,'from ',mySub,'to ',dest+1_lng,'at iter ',iter
!
!! Wait for all communication to complete
!CALL MPI_WAITALL(nreqs,parreqid(1:nreqs),parwtstat,mpierr)
!
!! Unpack data if recd. 
!IF (partransferrecv(1).GT.0_lng) THEN
!	CALL list_insert(ParListHead)
!	current => ParListHead%next
!	! PACK DATA INTO AN ARRAY
!	current%parid 		= INT(partransferrecv(1),lng)
!	!write(*,*) mySub,current%parid,partransferrecv(1),partransferrecv(5),partransferrecv(6)
!	!PAUSE
!	current%xp 		= partransferrecv(2)
!	current%yp 		= partransferrecv(3)
!	current%zp		= partransferrecv(4)
!	current%up 		= partransferrecv(5)
!	current%vp 		= partransferrecv(6)
!	current%wp 		= partransferrecv(7)
!	current%rp 		= partransferrecv(8)
!	current%delNBbyCV 	= partransferrecv(9)
!	current%par_conc 	= partransferrecv(10)
!	current%bulk_conc 	= partransferrecv(11)
!	current%rpold 		= partransferrecv(12)
!	current%sh 		= partransferrecv(13)
!	current%gamma_cont 	= partransferrecv(14)
!	current%cur_part 	= INT(partransferrecv(15),lng)
!	current%new_part 	= INT(partransferrecv(16),lng)
!
!	current%cur_part 	= current%new_part
!END IF
!
!! If data was sent then remove sent particle from list
!current => ParListHead%next
!DO WHILE (ASSOCIATED(current))
!	next => current%next ! copy pointer of next node
!	IF (current%cur_part /= current%new_part) THEN ! Transfer only if the particle has moved to another partition
!	!write(*,*) current%parid,current%cur_part,current%new_part
!	CALL list_delete(current)
!	END IF
!	! point to next node in the list
!	current => next
!	!write(*,*) current%parid
!ENDDO
!
!
!
!! Rough outline of the algorithm
!! PACK DATA INTO OPTIMAL NUMBER OF PACKETS
!! POST RECEIVE
!! SEND DATA
!! MPI_WAITALL - make sure all the data transfer is completed
!! UNPACK DATA AND ADD DATA TO THE PARTICLE LINKED LIST
!
!!------------------------------------------------
!END SUBROUTINE Particle_MPI_Transfer
!!------------------------------------------------
!!--------------------------------------------------------------------------------------------------
!SUBROUTINE Particle_MPI_Transfer	! transfer the distribution functions between neighboring initialize the MPI arrays and variables
!!--------------------------------------------------------------------------------------------------
!IMPLICIT NONE
!TYPE(ParRecord), POINTER :: current
!TYPE(ParRecord), POINTER :: next
!INTEGER(lng) :: nreqs		! number of send/recv requests
!INTEGER(lng) :: mpierr			! MPI standard error object
!LOGICAL :: probeflag
!
!
!nreqs = 0_lng
!
!write(*,*) 'In Particle_MPI	','FLAG=',ParticleTransfer,'SUBID = ',mySub,'iter = ',iter
!
!IF (ParticleTransfer) THEN
!	! Post all the send requests
!	current => ParListHead%next
!	DO WHILE (ASSOCIATED(current))
!		next => current%next ! copy pointer of next node
!		IF (current%cur_part /= current%new_part) THEN ! Transfer only if the particle has moved to another partition
!		! PACK DATA INTO AN ARRAY
!		partransfersend(1) = REAL(current%parid,dbl)
!		!write(*,*) mySub,current%parid,partransfersend(1)
!		partransfersend(2) = REAL(current%xp,dbl)
!		partransfersend(3) = REAL(current%yp,dbl)
!		partransfersend(4) = REAL(current%zp,dbl)
!		partransfersend(5) = REAL(current%up,dbl)
!		partransfersend(6) = REAL(current%vp,dbl)
!		partransfersend(7) = REAL(current%wp,dbl)
!		partransfersend(8) = REAL(current%rp,dbl)
!		partransfersend(9) = REAL(current%delNBbyCV,dbl)
!		partransfersend(10) = REAL(current%par_conc,dbl)
!		partransfersend(11) = REAL(current%bulk_conc,dbl)
!		partransfersend(12) = REAL(current%rpold,dbl)
!		partransfersend(13) = REAL(current%sh,dbl)
!		partransfersend(14) = REAL(current%gamma_cont,dbl)
!		partransfersend(15) = REAL(current%cur_part,dbl)
!		partransfersend(16) = REAL(current%new_part,dbl)
!    		nreqs = nreqs + 1_lng  
!    		!nreqs = 1_lng  
!		! SEND DATA
!		!CALL MPI_ISEND(partransfersend,16_lng,MPI_DOUBLE_PRECISION &
!		!,current%new_part-1_lng,current%parid,MPI_COMM_WORLD, &
!		!parreqid(nreqs),mpierr)	! send data
!		CALL MPI_SEND(partransfersend,16_lng,MPI_DOUBLE_PRECISION &
!		,current%new_part-1_lng,current%parid,MPI_COMM_WORLD, &
!		parreqid(nreqs),mpierr)	! send data
!		write(*,*) '**********sent','parid = ',current%parid,'from ',current%cur_part,'to ',current%new_part,'at iter ',iter,'in SUBID ',mySub
!		CALL list_delete(current)
!		END IF
!		!write(*,*) current%parid,current%cur_part,current%new_part
!		! point to next node in the list
!		current => next
!		!write(*,*) current%parid
!	ENDDO
!!	! delete all the nodes that are sent out
!!	current => ParListHead%next
!!	DO WHILE (ASSOCIATED(current))
!!		next => current%next ! copy pointer of next node
!!		IF (current%cur_part /= current%new_part) THEN ! Transfer only if the particle has moved to another partition
!!			CALL list_delete(current)
!!		END IF
!!		!write(*,*) current%parid,current%cur_part,current%new_part
!!		! point to next node in the list
!!		current => next
!!		!write(*,*) current%parid
!!	ENDDO
!
!ELSE
!! Now, we are not sure if this process wil rceive any data, so probe if there is any pending data. 
!	CALL MPI_IPROBE(MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,probeflag,probestat,mpierr)	
!	write(*,*) 'In MPI	','TAG=',probestat(MPI_TAG),'SOURCE =',probestat(MPI_SOURCE)+1_lng,'SUBID = ',mySub,'iter = ',iter
!	IF (probeflag) THEN
!    		nreqs = nreqs + 1_lng  
!	    	!nreqs = 2_lng  
!		!CALL MPI_IRECV(partransferrecv,16_lng,MPI_DOUBLE_PRECISION,probestat(MPI_SOURCE),probestat(MPI_TAG),MPI_COMM_WORLD,parreqid(nreqs),mpierr)
!		CALL MPI_RECV(partransferrecv,16_lng,MPI_DOUBLE_PRECISION,probestat(MPI_SOURCE),probestat(MPI_TAG),MPI_COMM_WORLD,parreqid(nreqs),mpierr)
!		!write(*,*) 'recd',probestat(MPI_TAG),probestat(MPI_SOURCE)+1_lng,mySub
!		CALL MPI_WAITALL(nreqs,parreqid(1:nreqs),parwtstat,mpierr)
!		write(*,*) '********recd  ','parid ',probestat(MPI_TAG),'from partition ',probestat(MPI_SOURCE)+1_lng,'in SUBID = ',mySub,'at iter ',iter
!		CALL list_insert(ParListHead)
!		current => ParListHead%next
!		! PACK DATA INTO AN ARRAY
!		current%parid 		= INT(partransferrecv(1),lng)
!		!write(*,*) mySub,current%parid,partransferrecv(1),partransferrecv(5),partransferrecv(6)
!		!PAUSE
!		current%xp 		= partransferrecv(2)
!		current%yp 		= partransferrecv(3)
!		current%zp		= partransferrecv(4)
!		current%up 		= partransferrecv(5)
!		current%vp 		= partransferrecv(6)
!		current%wp 		= partransferrecv(7)
!		current%rp 		= partransferrecv(8)
!		current%delNBbyCV 	= partransferrecv(9)
!		current%par_conc 	= partransferrecv(10)
!		current%bulk_conc 	= partransferrecv(11)
!		current%rpold 		= partransferrecv(12)
!		current%sh 		= partransferrecv(13)
!		current%gamma_cont 	= partransferrecv(14)
!		current%cur_part 	= INT(partransferrecv(15),lng)
!		current%new_part 	= INT(partransferrecv(16),lng)
!	
!		current%cur_part 	= current%new_part
!		write(*,*) partransferrecv(15),partransferrecv(16)
!		!	current%xpold 		= current%xp
!		!	current%ypold 		= current%yp
!		!	current%zpold 		= current%zp
!		!	current%upold 		= current%up
!		!	current%vpold 		= current%vp
!		!	current%wpold 		= current%wp
!	END IF
!
!END IF
!
!
!!CALL MPI_WAITALL(numReqs,req(1:numReqs),waitStat,mpierr)
!! POST RECEIVE
!		!CALL MPI_IRECV(partransferrecv,1_lng,MPI_DOUBLE_PRECISION,src,tag,MPI_COMM_WORLD,req(numReqs),mpierr)
!! MPI_WAITALL - make sure all the data transfer is completed
!! UNPACK DATA AND ADD DATA TO THE PARTICLE LINKED LIST
!	!CALL MPI_IRECV(partransferrecv,16_lng,MPI_DOUBLE_PRECISION,current%cur_part,tag,MPI_COMM_WORLD,reqid,mpierr) 
!
!
!
!! Rough outline of the algorithm
!! PACK DATA INTO OPTIMAL NUMBER OF PACKETS
!! POST RECEIVE
!! SEND DATA
!! MPI_WAITALL - make sure all the data transfer is completed
!! UNPACK DATA AND ADD DATA TO THE PARTICLE LINKED LIST
!
!!------------------------------------------------
!END SUBROUTINE Particle_MPI_Transfer
!!------------------------------------------------
!!--------------------------------------------------------------------------------------------------
!SUBROUTINE Particle_MPI_Transfer	! transfer the distribution functions between neighboring initialize the MPI arrays and variables
!!--------------------------------------------------------------------------------------------------
!IMPLICIT NONE
!REAL(dbl), ALLOCATABLE  :: partransfersend(:),partransferrecv(:)
!TYPE(ParRecord), POINTER :: current
!TYPE(ParRecord), POINTER :: next
!INTEGER(lng) :: nreqs		! number of send/recv requests
!INTEGER(lng),ALLOCATABLE :: reqid(:),wtstat(:,:)		! number of send/recv requests
!INTEGER(lng) :: mpierr			! MPI standard error object
!INTEGER(lng) :: probestat(MPI_STATUS_SIZE),recvstat(MPI_STATUS_SIZE)	! MPI status object
!LOGICAL :: probeflag
!
!ALLOCATE(partransfersend(16))
!ALLOCATE(partransferrecv(16))
!nreqs = 0_lng
!ALLOCATE(reqid(2))
!ALLOCATE(wtstat(MPI_STATUS_SIZE,2))
!
!!CALL MPI_IPROBE(MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,probeflag,probestat,mpierr)	
!!IF (probeflag) THEN
!!    	nreqs = nreqs + 1_lng  
!!	CALL MPI_IRECV(partransferrecv,16_lng,MPI_DOUBLE_PRECISION,probestat(MPI_SOURCE),probestat(MPI_TAG),MPI_COMM_WORLD,reqid(nreqs),mpierr)
!!END IF
!
!
!IF (ParticleTransfer) THEN
!	! Post all the send requests
!	current => ParListHead%next
!	DO WHILE (ASSOCIATED(current))
!		next => current%next ! copy pointer of next node
!		IF (current%cur_part /= current%new_part) THEN ! Transfer only if the particle has moved to another partition
!		! PACK DATA INTO AN ARRAY
!		partransfersend(1) = REAL(current%parid,dbl)
!		!write(*,*) mySub,current%parid,partransfersend(1)
!		partransfersend(2) = REAL(current%xp,dbl)
!		partransfersend(3) = REAL(current%yp,dbl)
!		partransfersend(4) = REAL(current%zp,dbl)
!		partransfersend(5) = REAL(current%up,dbl)
!		partransfersend(6) = REAL(current%vp,dbl)
!		partransfersend(7) = REAL(current%wp,dbl)
!		partransfersend(8) = REAL(current%rp,dbl)
!		partransfersend(9) = REAL(current%delNBbyCV,dbl)
!		partransfersend(10) = REAL(current%par_conc,dbl)
!		partransfersend(11) = REAL(current%bulk_conc,dbl)
!		partransfersend(12) = REAL(current%rpold,dbl)
!		partransfersend(13) = REAL(current%sh,dbl)
!		partransfersend(14) = REAL(current%gamma_cont,dbl)
!		partransfersend(15) = REAL(current%cur_part,dbl)
!		partransfersend(16) = REAL(current%new_part,dbl)
!    		nreqs = nreqs + 1_lng  
!    		!nreqs = 1_lng  
!		! SEND DATA
!	CALL MPI_ISEND(partransfersend,16_lng,MPI_DOUBLE_PRECISION &
!	,current%new_part-1_lng,current%parid,MPI_COMM_WORLD, &
!	reqid(nreqs),mpierr)	! send data
!		write(*,*) 'sent',current%parid,current%cur_part,current%new_part,iter
!		CALL list_delete(current)
!		END IF
!		!write(*,*) current%parid,current%cur_part,current%new_part
!		! point to next node in the list
!		current => next
!		!write(*,*) current%parid
!	ENDDO
!!	! delete all the nodes that are sent out
!!	current => ParListHead%next
!!	DO WHILE (ASSOCIATED(current))
!!		next => current%next ! copy pointer of next node
!!		IF (current%cur_part /= current%new_part) THEN ! Transfer only if the particle has moved to another partition
!!			CALL list_delete(current)
!!		END IF
!!		!write(*,*) current%parid,current%cur_part,current%new_part
!!		! point to next node in the list
!!		current => next
!!		!write(*,*) current%parid
!!	ENDDO
!END IF
!
!! Now, we are not sure if this process wil rceive any data, so probe if there is any pending data. 
!CALL MPI_IPROBE(MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,probeflag,probestat,mpierr)	
!IF (probeflag) THEN
!    	nreqs = nreqs + 1_lng  
!    	!nreqs = 2_lng  
!	CALL MPI_IRECV(partransferrecv,16_lng,MPI_DOUBLE_PRECISION,probestat(MPI_SOURCE),probestat(MPI_TAG),MPI_COMM_WORLD,reqid(nreqs),mpierr)
!	!write(*,*) 'recd',probestat(MPI_TAG),probestat(MPI_SOURCE)+1_lng,mySub
!	CALL MPI_WAITALL(nreqs,reqid(1:nreqs),wtstat,mpierr)
!	write(*,*) 'recd',probestat(MPI_TAG),probestat(MPI_SOURCE)+1_lng,mySub,iter
!	CALL list_insert(ParListHead)
!	current => ParListHead%next
!	! PACK DATA INTO AN ARRAY
!	current%parid 		= INT(partransferrecv(1),lng)
!	!write(*,*) mySub,current%parid,partransferrecv(1),partransferrecv(5),partransferrecv(6)
!	!PAUSE
!	current%xp 		= partransferrecv(2)
!	current%yp 		= partransferrecv(3)
!	current%zp		= partransferrecv(4)
!	current%up 		= partransferrecv(5)
!	current%vp 		= partransferrecv(6)
!	current%wp 		= partransferrecv(7)
!	current%rp 		= partransferrecv(8)
!	current%delNBbyCV 	= partransferrecv(9)
!	current%par_conc 	= partransferrecv(10)
!	current%bulk_conc 	= partransferrecv(11)
!	current%rpold 		= partransferrecv(12)
!	current%sh 		= partransferrecv(13)
!	current%gamma_cont 	= partransferrecv(14)
!	current%cur_part 	= INT(partransferrecv(15),lng)
!	current%new_part 	= INT(partransferrecv(16),lng)
!
!	current%cur_part 	= current%new_part
!	write(*,*) partransferrecv(15),partransferrecv(16)
!!	current%xpold 		= current%xp
!!	current%ypold 		= current%yp
!!	current%zpold 		= current%zp
!!	current%upold 		= current%up
!!	current%vpold 		= current%vp
!!	current%wpold 		= current%wp
!ENDIF
!
!DEALLOCATE(partransfersend)
!DEALLOCATE(partransferrecv)
!DEALLOCATE(reqid)
!DEALLOCATE(wtstat)
!
!!CALL MPI_WAITALL(numReqs,req(1:numReqs),waitStat,mpierr)
!! POST RECEIVE
!		!CALL MPI_IRECV(partransferrecv,1_lng,MPI_DOUBLE_PRECISION,src,tag,MPI_COMM_WORLD,req(numReqs),mpierr)
!! MPI_WAITALL - make sure all the data transfer is completed
!! UNPACK DATA AND ADD DATA TO THE PARTICLE LINKED LIST
!	!CALL MPI_IRECV(partransferrecv,16_lng,MPI_DOUBLE_PRECISION,current%cur_part,tag,MPI_COMM_WORLD,reqid,mpierr) 
!
!
!
!! Rough outline of the algorithm
!! PACK DATA INTO OPTIMAL NUMBER OF PACKETS
!! POST RECEIVE
!! SEND DATA
!! MPI_WAITALL - make sure all the data transfer is completed
!! UNPACK DATA AND ADD DATA TO THE PARTICLE LINKED LIST
!
!!------------------------------------------------
!END SUBROUTINE Particle_MPI_Transfer
!!------------------------------------------------

!--------------------------------------------------------------------------------------------------
SUBROUTINE MPI_Transfer	! transfer the distribution functions between neighboring initialize the MPI arrays and variables
!--------------------------------------------------------------------------------------------------
IMPLICIT NONE

! Define local variables
INTEGER(lng) :: iComm											! index variable
INTEGER(lng) :: numReqs											! number of send/recv requests
INTEGER(lng) :: mpierr											! MPI standard error object

CALL PackData														! fill out the arrays to be transferred 

numReqs = 0_lng

! Post the receives
DO iComm = 1,NumCommDirs

  IF(SubID(OppCommDir(iComm)) .NE. 0) THEN 
    numReqs = numReqs + 1_lng   
    CALL PostRecv(iComm,numReqs)										! receive data
  END IF

END DO

! Send the data
DO iComm = 1,NumCommDirs

  IF(SubID(iComm) .NE. 0) THEN
    numReqs = numReqs + 1_lng 
    CALL SendData(iComm,numReqs)						  				! send data
  END IF

END DO

CALL MPI_WAITALL(numReqs,req(1:numReqs),waitStat,mpierr)

! Store the Data
DO iComm = 1,NumCommDirs

  IF(SubID(OppCommDir(iComm)) .NE. 0) THEN    
    CALL UnPackData(iComm)										! store the sent data in the proper location at the destination
  END IF

END DO

!------------------------------------------------
END SUBROUTINE MPI_Transfer
!------------------------------------------------

!--------------------------------------------------------------------------------------------------
SUBROUTINE PostRecv(iComm,numReqs)	! receives information from a neighboring subdomain
!--------------------------------------------------------------------------------------------------
IMPLICIT NONE

INTEGER(lng), INTENT(IN) :: iComm								! communication direction from which to receive
INTEGER(lng), INTENT(IN) :: numReqs								! number of send/recv requests
INTEGER(lng) :: msgStart,msgEnd									! start and end indices of the message
INTEGER(lng) :: src													! source processing unit
INTEGER(lng) :: tag													! message tag
INTEGER(lng) :: mpierr												! MPI standard error variable

! starting/ending indices of the message in the mgsRecv array
msgStart = CommDataStart_f(OppCommDir(iComm))
msgEnd	= msgStart + msgSize(OppCommDir(iComm))

src	= SubID(OppCommDir(iComm)) - 1_lng						! rank of processing unit sending message TO this processing unit
tag	= iComm + 100_lng												! message tag 

CALL MPI_IRECV(msgRecv(msgStart:msgEnd),msgSize(OppCommDir(iComm)),MPI_DOUBLE_PRECISION,src,tag,MPI_COMM_WORLD,req(numReqs),mpierr)		! receive data

!------------------------------------------------
END SUBROUTINE PostRecv
!------------------------------------------------

!--------------------------------------------------------------------------------------------------
SUBROUTINE SendData(iComm,numReqs)											! sends information to a neighboring subdomain
!--------------------------------------------------------------------------------------------------
IMPLICIT NONE

INTEGER(lng), INTENT(IN) :: iComm								! communication direction
INTEGER(lng), INTENT(IN) :: numReqs								! number of send/recv requests
INTEGER(lng) :: msgStart,msgEnd									! start and end indices of the distribution function data
INTEGER(lng) :: dest													! rank of destination processing unit
INTEGER(lng) :: tag													! message tag
INTEGER(lng) :: mpierr												! MPI standard error variable 

! starting/ending indices of the message in the mgsSend array
msgStart		= CommDataStart_f(iComm)	
msgEnd		= msgStart + msgSize(iComm)

dest 			= SubID(iComm) - 1_lng								! rank of processing unit receiving message from the current processing unit (-1 to correspond to rank (myid))
tag			= iComm + 100_lng										! message tag 

CALL MPI_ISEND(msgSend(msgStart:msgEnd),msgSize(iComm),MPI_DOUBLE_PRECISION,dest,tag,MPI_COMM_WORLD,req(numReqs),mpierr)					! send data

!------------------------------------------------
END SUBROUTINE SendData
!------------------------------------------------

!--------------------------------------------------------------------------------------------------
SUBROUTINE UnPackData(iComm)	! store the recieved data in the proper locations
!--------------------------------------------------------------------------------------------------
IMPLICIT NONE

INTEGER(lng), INTENT(IN) :: iComm																				! communication direction
INTEGER(lng) :: i,j,k,m,ii																							! index variables
INTEGER(lng) :: stat(MPI_STATUS_SIZE)																			! MPI status object
INTEGER(lng) :: mpierr																								! MPI standard error variable

SELECT CASE(OppCommDir(iComm))
 
  CASE(1,2)			! YZ Faces

    i = YZ_RecvIndex(OppCommDir(iComm))																		! i index for obtaining the proper information from the 1D transfer array 	

    DO k=1,nzSub
      DO j=1,nySub

          ii = j + (k-1)*nySub																					! location of density and scalar function to recieve

          rho(i,j,k) = msgRecv((CommDataStart_rho(OppCommDir(iComm))-1) + ii)						! store recieved density in proper place
          phi(i,j,k) = msgRecv((CommDataStart_phi(OppCommDir(iComm))-1) + ii)						! store recieved scalar quantity in proper place          
          u(i,j,k) = msgRecv((CommDataStart_u(OppCommDir(iComm))-1) + ii)						! store recieved velocity quantity in proper place          
          v(i,j,k) = msgRecv((CommDataStart_v(OppCommDir(iComm))-1) + ii)						! store recieved velocity quantity in proper place          
          w(i,j,k) = msgRecv((CommDataStart_w(OppCommDir(iComm))-1) + ii)						! store recieved velocity quantity in proper place          
          node(i,j,k) = INT(msgRecv((CommDataStart_node(OppCommDir(iComm))-1) + ii),lng)					! store recieved node flag quantity in proper place          
        
        DO m=1,NumFs_Face
          
          ii = m 							&																			! location of distribution function to recieve
             + (j-1)*NumFs_Face 		&
             + (k-1)*nySub*NumFs_Face
             
          f(f_Comps(iComm,m),i,j,k) = msgRecv((CommDataStart_f(OppCommDir(iComm))-1) + ii)	! store the recieved distribution function component in the proper place in the f array

        END DO
        
      END DO
    END DO

  CASE(3,4)			! ZX Faces

    j = ZX_RecvIndex(OppCommDir(iComm))																		! j index for obtaining the proper information from the 1D transfer array 

    DO i=1,nxSub
      DO k=1,nzSub

          ii = k + (i-1)*nzSub																					! location of density and scalar function to recieve

          rho(i,j,k) = msgRecv((CommDataStart_rho(OppCommDir(iComm))-1) + ii)						! store recieved density in proper place
          phi(i,j,k) = msgRecv((CommDataStart_phi(OppCommDir(iComm))-1) + ii)						! store recieved scalar quantity in proper place              
          u(i,j,k) = msgRecv((CommDataStart_u(OppCommDir(iComm))-1) + ii)						! store recieved velocity quantity in proper place          
          v(i,j,k) = msgRecv((CommDataStart_v(OppCommDir(iComm))-1) + ii)						! store recieved velocity quantity in proper place          
          w(i,j,k) = msgRecv((CommDataStart_w(OppCommDir(iComm))-1) + ii)						! store recieved velocity quantity in proper place          
          node(i,j,k) = INT(msgRecv((CommDataStart_node(OppCommDir(iComm))-1) + ii),lng)					! store recieved node flag quantity in proper place          
        
        DO m=1,NumFs_Face
          
          ii = m 							&																			! location of distribution function to recieve
             + (k-1)*NumFs_Face 		&																			! location of density and scalar function to recieve
             + (i-1)*nzSub*NumFs_Face
             
          f(f_Comps(iComm,m),i,j,k) = msgRecv((CommDataStart_f(OppCommDir(iComm))-1) + ii)	! store the recieved distribution function component in the proper place in the f array

        END DO
        
      END DO
    END DO

  CASE(5,6)			! XY Faces

    k = XY_RecvIndex(OppCommDir(iComm))																		! k index for obtaining the proper information from the 1D transfer array 

    DO j=1,nySub
      DO i=1,nxSub

        ii = i + (j-1)*nxSub																						! location of density and scalar function to recieve
             	
          rho(i,j,k) = msgRecv((CommDataStart_rho(OppCommDir(iComm))-1) + ii)						! store recieved density in proper place
          phi(i,j,k) = msgRecv((CommDataStart_phi(OppCommDir(iComm))-1) + ii)						! store recieved scalar quantity in proper place         
          u(i,j,k) = msgRecv((CommDataStart_u(OppCommDir(iComm))-1) + ii)						! store recieved velocity quantity in proper place          
          v(i,j,k) = msgRecv((CommDataStart_v(OppCommDir(iComm))-1) + ii)						! store recieved velocity quantity in proper place          
          w(i,j,k) = msgRecv((CommDataStart_w(OppCommDir(iComm))-1) + ii)						! store recieved velocity quantity in proper place          
          node(i,j,k) = INT(msgRecv((CommDataStart_node(OppCommDir(iComm))-1) + ii),lng)					! store recieved node flag quantity in proper place          
        
        DO m=1,NumFs_Face
          
          ii = m 							&																			! location of distribution function to recieve
             + (i-1)*NumFs_Face 		&																			! location of density and scalar function to recieve
             + (j-1)*nxSub*NumFs_Face
             
          f(f_Comps(iComm,m),i,j,k) = msgRecv((CommDataStart_f(OppCommDir(iComm))-1) + ii)	! store the recieved distribution function component in the proper place in the f array
  
        END DO
        
      END DO
    END DO

  CASE(7:10)		! Z Sides

    i = Z_RecvIndex(OppCommDir(iComm),1)																		! i index for obtaining the proper information from the 1D transfer array 	
    j = Z_RecvIndex(OppCommDir(iComm),2)																		! j index for obtaining the proper information from the 1D transfer array 

    DO k=1,nzSub

      rho(i,j,k) = msgRecv((CommDataStart_rho(OppCommDir(iComm))-1) + k)							! store recieved density in proper place
      phi(i,j,k) = msgRecv((CommDataStart_phi(OppCommDir(iComm))-1) + k)							! store recieved scalar quantity in proper place            
      u(i,j,k) = msgRecv((CommDataStart_u(OppCommDir(iComm))-1) + k)						! store recieved velocity quantity in proper place          
      v(i,j,k) = msgRecv((CommDataStart_v(OppCommDir(iComm))-1) + k)						! store recieved velocity quantity in proper place          
      w(i,j,k) = msgRecv((CommDataStart_w(OppCommDir(iComm))-1) + k)						! store recieved velocity quantity in proper place          
      node(i,j,k) = INT(msgRecv((CommDataStart_node(OppCommDir(iComm))-1) + k),lng)					! store recieved node flag quantity in proper place          
        
      DO m=1,NumFs_Side
          
        ii = m + (k-1)*NumFs_Side																				! location of distribution function to recieve
             
        f(f_Comps(iComm,m),i,j,k) = msgRecv((CommDataStart_f(OppCommDir(iComm))-1) + ii)		! store the recieved distribution function component in the proper place in the f array

      END DO

    END DO

  CASE(11:14)		! X Sides

    j = X_RecvIndex(OppCommDir(iComm),1)																		! j index for obtaining the proper information from the 1D transfer array 	
    k = X_RecvIndex(OppCommDir(iComm),2)																		! k index for obtaining the proper information from the 1D transfer array 

    DO i=1,nxSub

      rho(i,j,k) = msgRecv((CommDataStart_rho(OppCommDir(iComm))-1) + i)							! store recieved density in proper place
      phi(i,j,k) = msgRecv((CommDataStart_phi(OppCommDir(iComm))-1) + i)							! store recieved scalar quantity in proper place  
      u(i,j,k) = msgRecv((CommDataStart_u(OppCommDir(iComm))-1) + i)						! store recieved velocity quantity in proper place          
      v(i,j,k) = msgRecv((CommDataStart_v(OppCommDir(iComm))-1) + i)						! store recieved velocity quantity in proper place          
      w(i,j,k) = msgRecv((CommDataStart_w(OppCommDir(iComm))-1) + i)						! store recieved velocity quantity in proper place          
      node(i,j,k) = INT(msgRecv((CommDataStart_node(OppCommDir(iComm))-1) + i),lng)					! store recieved node flag quantity in proper place          
        
      DO m=1,NumFs_Side
          
        ii = m + (i-1)*NumFs_Side																				! location of distribution function to recieve
             
        f(f_Comps(iComm,m),i,j,k) = msgRecv((CommDataStart_f(OppCommDir(iComm))-1) + ii)		! store the recieved distribution function component in the proper place in the f array

      END DO

    END DO

  CASE(15:18)		! Y Sides

    i = Y_RecvIndex(OppCommDir(iComm),1)																		! i index for obtaining the proper information from the 1D transfer array 	
    k = Y_RecvIndex(OppCommDir(iComm),2)																		! k index for obtaining the proper information from the 1D transfer array 

    DO j=1,nySub

      rho(i,j,k) = msgRecv((CommDataStart_rho(OppCommDir(iComm))-1) + j)							! store recieved density in proper place
      phi (i,j,k) = msgRecv((CommDataStart_phi(OppCommDir(iComm))-1) + j)							! store recieved scalar quantity in proper place  
      u(i,j,k) = msgRecv((CommDataStart_u(OppCommDir(iComm))-1) + j)						! store recieved velocity quantity in proper place          
      v(i,j,k) = msgRecv((CommDataStart_v(OppCommDir(iComm))-1) + j)						! store recieved velocity quantity in proper place          
      w(i,j,k) = msgRecv((CommDataStart_w(OppCommDir(iComm))-1) + j)						! store recieved velocity quantity in proper place          
      node(i,j,k) = INT(msgRecv((CommDataStart_node(OppCommDir(iComm))-1) + j),lng)					! store recieved node flag quantity in proper place          
        
      DO m=1,NumFs_Side
          
        ii = m + (j-1)*NumFs_Side																				! location of distribution function to recieve
             
        f(f_Comps(iComm,m),i,j,k) = msgRecv((CommDataStart_f(OppCommDir(iComm))-1) + ii)		! store the recieved distribution function component in the proper place in the f array

      END DO

    END DO

  CASE(19:26)		! Corners

    i = Corner_RecvIndex(OppCommDir(iComm),1) 																! i index for obtaining the proper information from the 1D transfer array 
    j = Corner_RecvIndex(OppCommDir(iComm),2)															 	! j index for obtaining the proper information from the 1D transfer array 	
    k = Corner_RecvIndex(OppCommDir(iComm),3)																! k index for obtaining the proper information from the 1D transfer array    

    rho(i,j,k) = msgRecv(CommDataStart_rho(OppCommDir(iComm)))											! store recieved density in proper place
    phi(i,j,k) = msgRecv(CommDataStart_phi(OppCommDir(iComm)))											! store recieved scalar quantity in proper place  
    u(i,j,k) = msgRecv(CommDataStart_u(OppCommDir(iComm)))						! store recieved velocity quantity in proper place          
    v(i,j,k) = msgRecv(CommDataStart_v(OppCommDir(iComm)))						! store recieved velocity quantity in proper place          
    w(i,j,k) = msgRecv(CommDataStart_w(OppCommDir(iComm)))						! store recieved velocity quantity in proper place          
    node(i,j,k) = INT(msgRecv(CommDataStart_node(OppCommDir(iComm))),lng)					! store recieved node flag quantity in proper place          
    f(f_Comps(iComm,1),i,j,k) = msgRecv(CommDataStart_f(OppCommDir(iComm)))						! store the recieved distribution function component in the proper place in the f array

  CASE DEFAULT

    OPEN(1000,FILE="error.txt")
    WRITE(1000,*) "Error in UnPackData in Parallel.f90: iComm is not 1-26..."
    WRITE(1000,*) "iComm",iComm
    CLOSE(1000)
    STOP

END SELECT

!------------------------------------------------
END SUBROUTINE UnPackData
!------------------------------------------------

!================================================
END MODULE Parallel
!================================================
