
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!
!!!
!!! File history:
!!!
!!! Date        Name                Comment
!!! -------------------------------------------------------------------
!!! 2004-02-02  Alvin Brattli       Created
!!!
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine colist(
     &              m, obs_lat, obs_lon, obs_sat, selectsat,
     &              n, em_lat, em_lon, 
     &              obs_ind, irc)
   
      implicit none

      integer, intent(in)    :: m           ! number of observations
      real,    intent(in)    :: obs_lat(m)  ! observation latitudes
      real,    intent(in)    :: obs_lon(m)  ! observation longitudes
      integer, intent(in)    :: obs_sat(m)  ! observation satellite no.
      integer, intent(in)    :: selectsat   ! which satellite to filter for
      integer, intent(in)    :: n           ! number of emissivities
      real,    intent(in)    :: em_lat(n)   ! emissivity latitudes
      real,    intent(in)    :: em_lon(n)   ! emissivity longitudes
      integer, intent(out)   :: obs_ind(m)  ! resultant obs -> em mapping
      integer, intent(inout) :: irc         ! error return code

      real,    parameter    :: maxdeg=0.5 ! max. distance (in degrees)

      ! Define a tree node
      type node2d
         real :: lon, lat           ! Point coordinates
         integer :: ind
         type(node2d), pointer :: l ! Left child node
         type(node2d), pointer :: r ! Right child node
      end type node2d
      
      ! Due to limitations in Fortran 90, we neeed this to make
      ! an array of pointers to tree nodes
      type node2dptr
         type(node2d), pointer :: p
      end type node2dptr
      
      ! Tree main data structure
      type tree2d
         integer N                 ! Number of elements in tree
         type(node2dptr), pointer :: list(:) ! Data points
         type(node2d),    pointer :: root ! Pointer to tree root node
      end type tree2d
      
      character(*), parameter :: myname = 'ASSIGN_OPERATIONS'

      type(tree2d) :: tree
      integer      :: ind
      integer      :: ii, jj, newn
      real         :: dist
      
      !! Put the emissivities in a 2D tree
      call tree_build2(tree, em_lon, em_lat, irc)
      if (irc /= 0) then
         write(*,*) myname, 'Error building tree:', irc
         return
      end if
      write(*,'("Tree depth: ",i7)') tree_depth(tree%root)
      write(*,'("Tree count: ",i7)') tree_count(tree%root)
      
      obs_ind(:) = 0
      do ii = 1, m
         if (obs_sat(ii) == selectsat) then
            dist = 360.0
            call rangesearch2d(dist, ind, tree%root,
     &                         obs_lon(ii), obs_lat(ii), maxdeg, 0)
            if (dist < maxdeg) obs_ind(ii) = ind
         end if
      end do
      
      call tree_destroy(tree)
      
      contains



      real function distance_deg(alon, alat, blon, blat)
      
      implicit none
      
      real :: alon
      real :: alat
      real :: blon
      real :: blat
      real :: d
      
      real, external :: cosdeg, sindeg, acosdeg
      
      d = sindeg(alat)*sindeg(blat) + 
     &     cosdeg(alat)*cosdeg(blat)*cosdeg(alon-blon)
      d = max(-1.0, min(1.0, d))

      distance_deg = abs(acosdeg(d))

      end function distance_deg



      recursive subroutine rangesearch2d(dist, distind, node,
     &                                   lon, lat, maxdist, dir)
      
      implicit none
      
      real,         intent(inout) :: dist
      integer,      intent(inout) :: distind
      type(node2d), pointer    :: node
      real,         intent(in) :: lon
      real,         intent(in) :: lat
      real,         intent(in) :: maxdist
      integer,      intent(in) :: dir
      
      real :: curdist
      logical :: goleft, goright
      integer :: newdir

      if (.not.associated(node)) return

      curdist = distance_deg(node%lon, node%lat, lon, lat)
      if (curdist < dist) then
         dist = curdist
         distind = node%ind
      end if
      
      !!! Determine which subtrees (if any) we must recursively search in
      if (dir == 0) then
         goleft  = lon-maxdist < node%lon
         goright = node%lon <= lon+maxdist
      else
         goleft  = lat-maxdist < node%lat
         goright = node%lat <= lat+maxdist
      end if

      !!! Change sorting direction for next level
      newdir = mod(dir+1, 2)
      if (goleft.and.associated(node%l)) then
         call rangesearch2d(dist, distind, node%l, 
     &                      lon, lat, maxdist, newdir)
      end if
      if (goright.and.associated(node%r)) then
         call rangesearch2d(dist, distind, node%r,
     &                      lon, lat, maxdist, newdir)
      end if
      
      end subroutine rangesearch2d



      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!
      !!! Name: tree_build2
      !!!
      !!! Description:
      !!! 
      !!!   Build a 2D tree.  To optimize speed of both insertion and
      !!!   searching, wee need the tree to be as balanced as possible,
      !!!   with minimum depth.  We attempt to reach this goal by
      !!!   finding the point with median longitude, using this as our
      !!!   root element.  This ensures that left and right sides of
      !!!   the tree contains the same number of elements, plus or minus
      !!!   one or two.  Further, we insert the rest of the points
      !!!   randomly, attempting to balance left and right subtrees
      !!!   properly.
      !!!
      !!! Input:   tree - tree information
      !!!          lon, lat - longitude and latitude information
      !!!
      !!! Output:  tree - tree information, now initialised and
      !!!                 x and y points inserted
      !!!          irc - error return code
      !!!

      subroutine tree_build2(tree, lon, lat, irc)
      
      implicit none
      
      type(tree2d), intent(inout) :: tree ! tree trunk information
      real,         intent(in)    :: lon(:)
      real,         intent(in)    :: lat(size(lon))
      integer,      intent(inout) :: irc
      
      character(*), parameter :: myname = 'TREE_BUILD2'

      integer :: ind(size(lon))
      integer :: i, j, median
      type(node2d), pointer :: median_node
      integer :: n  ! number of entries to sort
      
      n = size(lon)
      
      ! Prepare the tree structure
      tree%N = 0
      nullify(tree%root)
      if (n < 1) return
      allocate(tree%list(n), stat=irc)
      if (irc /= 0) then
         write(*,*) myname, "Couldn't allocate memory for list: ", irc
         irc = 901
         return
      end if
      
      ! First insert element with median longitude.  This is our root node.
      ind(1:n) = (/(i, i = 1, n)/)
      call select_nth(ind, lon, (n+1)/2)
      median = ind((n+1)/2)
      tree%N = tree%N + 1
      call node_insert(median_node, tree,
     &                 lon(median), lat(median), median, irc)
      if (irc /= 0) then
         write(*,*) myname, 'Insertion of first element failed'
         return
      end if
      
      ! In order to get a properly balanced tree, we want to insert the rest
      ! of the elements in random order, using iind() to index ind()
      call randomize_array(ind)
      
      do i = 1, n
         j = ind(i)
         ! Special case for median element
         if (j == median) then
            tree%list(j)%p => median_node
            cycle
         end if
         tree%N = tree%N + 1
         call node_insert(tree%list(j)%p,tree,lon(j),lat(j),j,irc)
!         tree%list(j)%p => node_insert(tree%root,lon(j),lat(j),j,0,irc)
         if (irc /= 0) then
            write(*,*) myname, 'Insertion of element ', i, ' failed'
            return
         end if
      end do
      
      end subroutine tree_build2




      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!
      !!! Name: randomize_array
      !!!
      !!! Description: relocate elements in the array at random
      !!!
      !!! Input: ind - array to randomize
      !!! Output: ind - randomized array
      !!!

      subroutine randomize_array(ind)
      
      implicit none
      
      integer :: ind(:)
      
      real :: rnd
      integer :: n
      integer :: i, j, tmp
      
      n = size(ind)
      do i = 1, n
         call random_number(rnd)
         j = int(1 + real(n)*rnd)
         j = min(max(j, 1), n)  ! Ensure 1 <= j <= n
         tmp = ind(i)
         ind(i) = ind(j)
         ind(j) = tmp
      end do
      
      end subroutine randomize_array

      recursive integer function tree_depth(tree) result(depth)

      implicit none

      type(node2d), pointer :: tree

      if (.not.associated(tree)) then
         depth = 0
         return
      end if

      depth = 1 + max(tree_depth(tree%l), tree_depth(tree%r))

      end function tree_depth


      recursive integer function tree_count(tree) result(res)

      implicit none

      type(node2d), pointer :: tree

      if (.not.associated(tree)) then
         res = 0
         return
      end if

      res = 1 + tree_count(tree%l) + tree_count(tree%r)

      end function tree_count
      
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!
      !!! Name: node_insert
      !!!
      !!! Description: Insert a number pair into tree.
      !!!
      !!! Input:       tree - tree to insert into
      !!!              lon, lat - number pair to insert
      !!!              dir - 0 for lon-branching, 1 for lat-branching
      !!!
      !!! Output:      root - modified tree
      !!!              irc - error return code
      !!!
      !!! Returns:     new - pointer to new node
      !!!

      subroutine node_insert(new, tree, lon, lat, ind, irc)
      
      implicit none
      
      type(node2d), pointer       :: new
      type(tree2d), intent(inout) :: tree
      real,         intent(in)    :: lon
      real,         intent(in)    :: lat
      integer,      intent(in)    :: ind
      integer,      intent(out)   :: irc
      
      character(len=15), parameter :: myname = 'NODE_INSERT'

      type(node2d), pointer :: node, prevnode
      integer               :: dir
      logical               :: goleft
      
      dir = 0
      node => tree%root

      do while (associated(node))
         if (dir == 0) then
            goleft = lon < node%lon
         else
            goleft = lat < node%lat
         end if
         prevnode => node
         if (goleft) then
            node => node%l
         else
            node => node%r
         end if
         dir = mod(dir+1, 2)
      end do

      allocate(new, stat=irc)
      if (irc /= 0) then
         write(*,*) myname, "Couldn't allocate memory: ", irc
         irc = 900
         return
      end if
      new%lon = lon
      new%lat = lat
      new%ind = ind
      nullify(new%l, new%r)

      if (associated(tree%root)) then
         if (goleft) then
            prevnode%l => new
         else
            prevnode%r => new
         end if
      else
         tree%root => new
      end if

      end subroutine node_insert

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!
      !!! Name: tree_destroy
      !!!
      !!! Description: Destroy a 2d tree
      !!!
      !!! Input:       tree - tree to chop down
      !!!
      !!! Output:      tree - chopped down tree
      !!!

      subroutine tree_destroy(tree)
      
      implicit none
      type(tree2d) :: tree
      integer :: irc
      
      tree%N = 0
      if (associated(tree%list)) then
         deallocate(tree%list, stat=irc)
         nullify(tree%list)
      end if
      call branch_chop(tree%root)
      
      end subroutine tree_destroy
      
      
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!
      !!! Name: branch_chop
      !!!
      !!! Description: Recursively traverse and deallocate
      !!!              all branches and nodes in a (sub)tree
      !!!
      !!! Input:       tree - tree to chop down
      !!!
      !!! Output:      tree - chopped down tree
      !!!

      recursive subroutine branch_chop(tree)
      
      implicit none
      type(node2d), pointer :: tree
      integer :: irc
      
      if (associated(tree)) then
         ! Remove left branch
         call branch_chop(tree%l)
         ! Remove right branch
         call branch_chop(tree%r)
         ! Remove root node and disassociate it
         deallocate(tree,stat=irc)
         nullify(tree)
      end if
      
      end subroutine branch_chop


      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!
      !!! Name:  select_nth
      !!!
      !!! Description: Select the Kth smallest element (out of n) in
      !!!              the key array.  Leaves ind() semi-sorted,
      !!!              with the Kth smallest element in key(ind(k)).
      !!!              Leaves ind() in the following state:
      !!!              key(ind(0:k-1)) < key(ind(k)) < key(ind(k+1:n))
      !!!              For more information, see Sedgewick: _Algorithms_
      !!!
      !!! Input:       ind - array of indexes into to the key array
      !!!              key - what to sort by
      !!!              k   - which element to find
      !!!
      !!! Output:      ind - partially sorted array of indexes
      !!!

      subroutine select_nth(ind, key, k)
      
      implicit none
      
      ! Parameters
      integer ind(:)            ! array of indexes
      real key(:)               ! sort key
      integer k
      
      ! Local variables
      integer n, t, i, j, l, r,nv
      real v
      
      n = size(ind)
      l = 1
      r = n
      do while (r > l)
         v = key(ind(r))
         i = l-1
         j = r
         partition: do
            do 
               i = i+1
               if (key(ind(i)) >= v) exit
            end do
            do
               j = j-1
               if (j <= i) exit
               if (key(ind(j)) <= v) exit
            end do
            t = ind(i)
            ind(i) = ind(j)
            ind(j) = t
            if (j <= i) exit partition
         end do partition
         ind(j) = ind(i)
         ind(i) = ind(r)
         ind(r) = t
         if (i >= k) r = i-1
         if (i <= k) l = i+1
      end do

      end subroutine select_nth



      end subroutine colist

#ifdef JNK
      subroutine colist(totobs_sub,
     &     rtm_lat,
     &     rtm_lon,
     &     rtm_satins,
     &     isatins,
     &     maxtemp,
     &     tmp_lat,
     &     tmp_lon,
     &     ind_e,
     &     irc)
C
      implicit none
      integer totobs_sub       ! number of rtm obs
      real rtm_lat(totobs_sub) ! rtm obs latitude
      real rtm_lon(totobs_sub) ! rtm obs longtitude
      integer rtm_satins(totobs_sub) ! rtm obs satellite identifier
      integer isatins          ! emissivity satellite identifier
      integer maxtemp          ! number of emissivity
      real tmp_lat(maxtemp)    ! emissivity latitude
      real tmp_lon(maxtemp)    ! emissivity longtitude
      integer ind_e(totobs_sub) ! OUTPUT index for which emissivity
                               ! corresponds to given a obs, 
                               ! rtm_e(j)=tmp_e(ind_e(j))
      integer irc              ! error return code, 0=ok

C
      real mdist,cdist
      integer iind, ii,jj
C
      do ii=1,totobs_sub
         if (rtm_satins(ii).eq.isatins) then
            mdist=0.0
            iind=0
            do jj=1,maxtemp
               call dddist(rtm_lon(ii),rtm_lat(ii),
     &              tmp_lon(jj),tmp_lat(jj),
     &              cdist,irc)
               if (irc.ne.0) then
                  write(*,*) 'COLIST Error return from DIST.',irc
                  return
               end if
               if (jj.eq.1) then
                  mdist=cdist
                  iind=jj
               else
                  if (mdist.gt.cdist) then
                     mdist=cdist
                     iind=jj
                  end if
               end if
            end do
            ind_e(ii)=iind
         else
            ind_e(ii)=0
         end if
      end do
C
      return
      end
#endif



