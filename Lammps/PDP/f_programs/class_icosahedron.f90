MODULE class_icosahedron

   !Tools for creation of global finite element grids by subdivision of the
   !facets of the icosahedron into quasi-equilateral spherical triangles.
   !
   !by Peter Bird, UCLA; written in BASIC ~1982; translated to Fortran 90 2001.11
   !----------------------------------------------------------------------------
     PRIVATE
   !
   ! USER ROUTINES:
   !
                   PUBLIC Make_Global_Grid
                   PUBLIC Write_Global_Grid
   !
   ! UTILITY ROUTINES (called by the user routines):
   !
   !     RECURSIVE SUBROUTINE Divide
   !               SUBROUTINE Element_Out
   !               SUBROUTINE Unitize
   !----------------------------------------------------------------------------

CONTAINS

    SUBROUTINE Make_Global_Grid (n_slice, &           ! only input(!)
                               & numnod, node_uvec, & ! output: number of nodes, unit vectors of nodes,
                               & numel, nodes)        !         number of elements, element definitions
       !---------------------------------------------------------------
        !Generate a finite element grid of spherical triangles by subdivision of the icosahedron:
        !Level 0 (generated here) has 12 vertices, 30 edges, and 20 triangles.
        !Then, subdivide each face n_slice times, and output triangular elements,
        !  by using SUBR Divide, which calls itself recursively!
       !---------------------------------------------------------------
        IMPLICIT NONE
        INTEGER,                 INTENT(IN)  :: n_slice   ! subdivision level; 0 or higher
        INTEGER,                 INTENT(OUT) :: numnod    ! number of nodes created
        REAL, DIMENSION(:,:),    INTENT(OUT) :: node_uvec ! (3, numnod)
        INTEGER,                 INTENT(OUT) :: numel     ! number of elements created
        INTEGER, DIMENSION(:,:), INTENT(OUT) :: nodes     ! (3, numel)
       !---------------------------------------------------------------
        INTEGER                            :: m_slice ! copy of n_slice, allowing it to be changed (counted-down)
        INTEGER                            :: facets_done, i, j, k
        INTEGER, DIMENSION(3)              :: node_number
        REAL, DIMENSION(3)                 :: rx, ry, rz, uvec1, uvec2, uvec3
        DOUBLE PRECISION, PARAMETER        :: s = 1.107148719D0
        DOUBLE PRECISION                   :: dot1, dot2, dot3, x1, x2, x3, y1, y2, y3, z1, z2, z3
        DOUBLE PRECISION, DIMENSION(12)    :: lat, lon
        DOUBLE PRECISION, DIMENSION(3)     :: v1, v2, v3
        DOUBLE PRECISION, DIMENSION(3, 12) :: abg  !Cartesian (alpha, beta, gamma) coordinates of these vertices.
        LOGICAL counterclockwise
       !---------------------------------------------------------------
       !generate basic form with a vertex (a 5-fold axis) up; highest symmetry axis
        lat(1) = 1.570796327D0
        lon(1) = 0.0D0
        DO i = 2, 6
            lat(i) = lat(1) - s
            lon(i) = (i - 2.0D0) * 1.256637061D0
        END DO
        DO i = 7, 11
            lat(i) = -lat(1) + s
            lon(i) = (i - 7.0D0) * 1.256637061D0 + .628318531D0
        END DO
        lat(12) = -lat(1)
        lon(12) = 0.0D0
        DO i = 1, 12
            abg(1, i) = COS(lat(i)) * COS(lon(i))
            abg(2, i) = COS(lat(i)) * SIN(lon(i))
            abg(3, i) = SIN(lat(i))
        END DO
       !-------------------------------------------------------
       !create output file for dumping results as they are found:
        OPEN (UNIT = 729, FILE = "Icosahedron.tmp", FORM = "UNFORMATTED") ! unconditional; overwrites any existing file
       !-------------------------------------------------------
       !find all 20 faces and subdivide each into four spherical triangles;
       !WRITE (*, "(' Creating global grid by level-',I1,' subdivision of icosahedron facets:')") n_slice
       !WRITE (*, *) ! advance, because next WRITE will not
        m_slice = n_slice
        facets_done = 0
        DO i = 1, 10
            DO j = (i + 1), 11
                DO k = (j + 1), 12
                    dot1 = abg(1, i) * abg(1, j) + abg(2, i) * abg(2, j) + abg(3, i) * abg(3, j)
                    dot2 = abg(1, j) * abg(1, k) + abg(2, j) * abg(2, k) + abg(3, j) * abg(3, k)
                    dot3 = abg(1, k) * abg(1, i) + abg(2, k) * abg(2, i) + abg(3, k) * abg(3, i)
                    IF ((dot1 > 0.3D0) .AND. (dot2 > 0.3D0) .AND. (dot3 > 0.3D0)) THEN
                        x1 = abg(1, i)
                        x2 = abg(1, j)
                        x3 = abg(1, k)
                        y1 = abg(2, i)
                        y2 = abg(2, j)
                        y3 = abg(2, k)
                        z1 = abg(3, i)
                        z2 = abg(3, j)
                        z3 = abg(3, k)
                        !Note: Divide will call itself, recursively.
                        !Therefore, all inputs are simple numbers (not vectors),
                        !to go on stack as values, not addresses.
                        !Also, note that output is sent to a temporary file, to avoid multiple copies on stack!
                        CALL Divide(x1, y1, z1, x2, y2, z2, x3, y3, z3, m_slice) ! using copy m_slice = n_slice
                        facets_done = facets_done + 1
                       !WRITE (*, "('+',I8,' facets   out of       20 divided into elements')") facets_done
                    END IF
                END DO
            END DO
        END DO
       !-----------------------------------------------------------------------
       !read binary file, extracting groups of 3 uvecs, and assigning to lists:
        CLOSE (UNIT = 729)
        OPEN  (UNIT = 729, FILE = "Icosahedron.tmp", STATUS = "OLD", FORM = "UNFORMATTED")
        numnod = 0  !initialization
        IF (n_slice == 0) THEN
            numel = 20
        ELSE
            numel = 20 * (4**n_slice)
        END IF
        WRITE (*, *) ! advance, because next WRITE will not
        DO i = 1, numel
            READ (729) rx(1), ry(1), rz(1), rx(2), ry(2), rz(2), rx(3), ry(3), rz(3)
            DO j = 1, 3
                node_number(j) = 0 ! initialization
                k_loop: DO k = 1, numnod
                    IF (rx(j) == node_uvec(1, k)) THEN        
                        IF (ry(j) == node_uvec(2, k)) THEN        
                            IF (rz(j) == node_uvec(3, k)) THEN        
                               !this node is already defined
                                node_number(j) = k
                                EXIT k_loop
                            END IF
                        END IF
                    END IF
                END DO k_loop ! k = 1, numnod
                IF (node_number(j) == 0) THEN
                   !no match was found; define a new node
                    numnod = numnod + 1
                    node_number(j) = numnod
                    node_uvec(1, numnod) = rx(j)
                    node_uvec(2, numnod) = ry(j)
                    node_uvec(3, numnod) = rz(j)
                END IF
               !record this element, after checking for counterclockwise ordering!
               !First, determine side from node "#1" to "#2":
                uvec3(1) = rx(2) - rx(1)
                uvec3(2) = ry(2) - ry(1)
                uvec3(3) = rz(2) - rz(1)
               !Next, determine side from node "#2" to "#3":
                uvec1(1) = rx(3) - rx(2)
                uvec1(2) = ry(3) - ry(2)
                uvec1(3) = rz(3) - rz(2)
                CALL Cross(uvec3, uvec1, uvec2)
               !Note: If numbering is counterclockwise, uvec2 should point up.
                counterclockwise = ((uvec2(1)*rx(1) + uvec2(2)*ry(1) + uvec2(3)*rz(1)) > 0.0)
                IF (counterclockwise) THEN
                    nodes(1:3, i) = node_number(1:3)
                ELSE
                    nodes(1, i) = node_number(1)
                    nodes(2, i) = node_number(3)
                    nodes(3, i) = node_number(2)
                END IF
            END DO ! j = 1, 3
            IF (MOD(i, 100) == 0) THEN
               !WRITE (*, "('+',I8,' elements out of ',I8,' scanned for new nodes')") i, numel
            END IF
        END DO ! i = 1, numel
       !WRITE (*, "('+',I8,' elements out of ',I8,' scanned for new nodes')") numel, numel
        CLOSE (UNIT = 729, STATUS = "DELETE")
    END SUBROUTINE Make_Global_Grid

    SUBROUTINE Write_Global_Grid (path, &              ! [Drive:][\Path\] for output (or null)
                                & n_slice, &           
                                & numnod, node_uvec, & 
                                & numel, nodes)        ! all of these are INTENT(IN); get values 
                                                       ! by calling Make_Global_Grid first
        IMPLICIT NONE
        CHARACTER*(*),           INTENT(IN) :: path      ! [Drive:][\Path\] for output (or null)
        INTEGER,                 INTENT(IN) :: n_slice   ! subdivision level; 0 or higher
        INTEGER,                 INTENT(IN) :: numnod    ! number of nodes
        REAL, DIMENSION(:,:),    INTENT(IN) :: node_uvec ! (3, numnod)
        INTEGER,                 INTENT(IN) :: numel     ! number of elements
        INTEGER, DIMENSION(:,:), INTENT(IN) :: nodes     ! (3, numel)
       !---------------------------------------------------------------
        CHARACTER*1  :: c1
        CHARACTER*11 :: grid_file
        CHARACTER*80 :: grid_pathfile
        INTEGER :: i
        REAL, PARAMETER :: Pi = 3.141592654
        REAL, PARAMETER :: degrees_per_radian = 57.2957795
        REAL :: equat, equat2, lat, lon, phi, theta
        REAL, DIMENSION(3) :: uvec
       !---------------------------------------------------------------
        WRITE (c1, "(I1)") n_slice
        grid_file = "global" // c1 // ".feg"
        grid_pathfile = TRIM(path) // grid_file
        WRITE (*, "(' Writing ',A,'...')") TRIM(grid_file)
        OPEN (UNIT = 730, FILE = grid_pathfile) ! unconditional; overwrites any older file
        WRITE (730, "(A)") grid_file
        WRITE (730, "(I6,I6,'     0 1000000 F')") numnod, numnod
        DO i = 1, numnod
            uvec(1:3) = node_uvec(1:3, i)
           !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
           !Note: NOT calling Uvec_2_LonLat, because that is in Sphere or Map_Projections/Adobe_Illustrator,
           !      and in many projects would cause confusing duplication of code.  However, steps exactly the same.
            equat2 = uvec(1)*uvec(1) + uvec(2)*uvec(2)
            IF (equat2 == 0.) THEN
                phi = 0. ! actually undefined; default 0.
                IF (uvec(3) > 0.) THEN
                    theta = 0. ! N pole
                ELSE 
                    theta = Pi ! S pole
                END IF
            ELSE
                equat = SQRT(equat2)
                theta = ATAN2(equat, uvec(3))
                phi = ATAN2(uvec(2), uvec(1))
            END IF
            lat = 90. - degrees_per_radian * ABS(theta)
            lat = MAX (lat, -90.)
            lon = degrees_per_radian * phi
            IF (lon > 180.) lon = lon - 360.
            IF (lon <= -180.) lon = lon + 360.
           !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                WRITE (730, "(I6,F9.3,F8.3,'       0.0       0.0       0.0       0.0')") i, lon, lat
!                write(*,*) i,node_uvec(:,i)
        END DO
        WRITE (730, *) numel
        DO i = 1, numel
            WRITE (730, "(4I6)") i, nodes(1, i), nodes(2, i), nodes(3, i)
        END DO
        WRITE (730, "('     0')")
        CLOSE (UNIT = 730)
        WRITE (*, "('+Writing ',A,'...DONE')") grid_file
    END SUBROUTINE Write_Global_Grid

    SUBROUTINE Cross (a_vec, b_vec, c_vec)
      ! vector cross product: a x b = c
        IMPLICIT NONE
        REAL, DIMENSION(3), INTENT(IN)  :: a_vec, b_vec
        REAL, DIMENSION(3), INTENT(OUT) :: c_vec
        c_vec(1) = a_vec(2)*b_vec(3) - a_vec(3)*b_vec(2)
        c_vec(2) = a_vec(3)*b_vec(1) - a_vec(1)*b_vec(3)
        c_vec(3) = a_vec(1)*b_vec(2) - a_vec(2)*b_vec(1)
    END SUBROUTINE Cross

    RECURSIVE SUBROUTINE Divide (x1, y1, z1, x2, y2, z2, x3, y3, z3, m_slice)
      !CALLed only by SUBR Make_Global_Grid (but called multiple times) and by itself (ditto).
      !Accepts 3 unit vectors of triangle corners as input, subdivides the
      !triangle m_slice times into 4 triangles, and outputs results to a file.
      !NOTE that Divide calls itself recursively, so all the arguments should
      !be simple numbers passed by value, not by address.
      !To avoid complications, output is dumped into an open binary file (UNIT = 729).
      !-------------------------------------------------------
       IMPLICIT NONE
       DOUBLE PRECISION, INTENT(IN)    :: x1, y1, z1, x2, y2, z2, x3, y3, z3
       INTEGER,          INTENT(IN)    :: m_slice
      !-------------------------------------------------------
       INTEGER :: l_slice
       DOUBLE PRECISION :: xe1, xe2, xe3, ye1, ye2, ye3, ze1, ze2, ze3
      !-------------------------------------------------------
       xe1 = 0.5D0 * (x2 + x3)
       xe2 = 0.5D0 * (x3 + x1)
       xe3 = 0.5D0 * (x1 + x2)
       ye1 = 0.5D0 * (y2 + y3)
       ye2 = 0.5D0 * (y3 + y1)
       ye3 = 0.5D0 * (y1 + y2)
       ze1 = 0.5D0 * (z2 + z3)
       ze2 = 0.5D0 * (z3 + z1)
       ze3 = 0.5D0 * (z1 + z2)
       CALL Unitize(xe1, ye1, ze1)
       CALL Unitize(xe2, ye2, ze2)
       CALL Unitize(xe3, ye3, ze3)
       IF (m_slice > 1) THEN
           l_slice = m_slice - 1
           CALL Divide( x1,  y1,  z1, xe2, ye2, ze2, xe3, ye3, ze3, l_slice)
           CALL Divide( x2,  y2,  z2, xe1, ye1, ze1, xe3, ye3, ze3, l_slice)
           CALL Divide( x3,  y3,  z3, xe1, ye1, ze1, xe2, ye2, ze2, l_slice)
           CALL Divide(xe1, ye1, ze1, xe2, ye2, ze2, xe3, ye3, ze3, l_slice)
       ELSE IF (m_slice == 1) THEN
           CALL ElementOut( x1,  y1,  z1, xe2, ye2, ze2, xe3, ye3, ze3)
           CALL ElementOut( x2,  y2,  z2, xe1, ye1, ze1, xe3, ye3, ze3)
           CALL ElementOut( x3,  y3,  z3, xe1, ye1, ze1, xe2, ye2, ze2)
           CALL ElementOut(xe1, ye1, ze1, xe2, ye2, ze2, xe3, ye3, ze3)
       ELSE ! m_slice == 0
           CALL ElementOut(x1, y1, z1, x2, y2, z2, x3, y3, z3)
       END IF
    END SUBROUTINE Divide



    SUBROUTINE ElementOut (x1, y1, z1, x2, y2, z2, x3, y3, z3)
       !CALLed only by SUBR Divide (but called multiple times).
       !Outputs an element defined by 3 unit radius vectors to its vertices.
       !The output (binary) file must already be open, as device 729.
        IMPLICIT NONE
        DOUBLE PRECISION, INTENT(IN) :: x1, y1, z1, x2, y2, z2, x3, y3, z3
        REAL :: rx1, rx2, rx3, ry1, ry2, ry3, rz1, rz2, rz3
       !Convert to single-precision:
        rx1 = x1; rx2 = x2; rx3 = x3
        ry1 = y1; ry2 = y2; ry3 = y3
        rz1 = z1; rz2 = z2; rz3 = z3
       !Write to binary file:
        WRITE (729) rx1, ry1, rz1, rx2, ry2, rz2, rx3, ry3, rz3
    END SUBROUTINE ElementOut



    SUBROUTINE Unitize (x, y, z)
       !Converts any vector in 3-component double-precision form to a unit vector.
        IMPLICIT NONE
        DOUBLE PRECISION, INTENT(INOUT) :: x, y, z
        DOUBLE PRECISION :: r
        r = SQRT(x * x + y * y + z * z)
        IF (r > 0.0D0) THEN
            x = x / r
            y = y / r
            z = z / r
        ELSE
            x = 1.0D0
            y = 0.0D0
            z = 0.0D0
        END IF
    END SUBROUTINE Unitize

END MODULE 
