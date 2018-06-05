      program nullspace_basis
      implicit none

      type rowvale
      real, pointer :: prow(:)
      end type 
      type(rowvale), dimension(:), allocatable :: lp_vale

      type rowlign
      integer, pointer :: qrow(:)
      end type 
      type(rowlign), dimension(:), allocatable :: lp_lign

      type rowlong
      integer, pointer :: lrow
      end type 
      type(rowlong), dimension(:), allocatable :: lp_long
      real*8, dimension(:), allocatable :: slu_vale,L1_vale
      real*8, dimension(:), allocatable :: lpva, L2_vale, l_Tvalues
      real*8, dimension(:), allocatable ::  values, l_b, tval
      real*8 :: t1, t2, tol, temp
      integer, dimension(:), allocatable :: slu_PermL,slu_PermC
      integer, dimension(:), allocatable :: rowind,colptr, trow,tcolp
      integer, dimension(:), allocatable :: slu_rowind,slu_diagu
      integer, dimension(:), allocatable :: ibid1, ibid2, lpli
      integer, dimension(:), allocatable :: slu_colptr, lp_long1
      integer, dimension(:), allocatable :: L1_colptr, l_Tcolptr
      integer, dimension(:), allocatable :: L1_rowind, L2_colptr
      integer, dimension(:), allocatable :: L2_rowind, l_Trowind
      integer ::  m, n, nnz, nrhs, ldb, info, iopt, nbnz, nnzT
      integer ::  pos, posK, i, nblag, nbphys, l1,j, igap
      integer :: nbz_diag, nbnz_diag, i1, ii,j1,k1
      integer :: nzL1,nzL2,ind_lig, ind_col
      integer :: reste, seuil, inc, nnzK
      integer :: totcrd, ptrcrd, indcrd, valcrd, rhscrd, neltvl
      integer :: totcrdw, ptrcrdw, indcrdw, valcrdw, rhscrdw, neltvlw
      integer :: j2, jpos, x
      logical :: yesa
      character*16 :: indfmt, ptrfmt
      character*16 :: indfmtw, ptrfmtw
      character*8 :: key
      character*8 :: keyw
      character*3 :: mxtype
      character*3 :: mxtypew
      character*20 :: rhsfmt, valfmt
      character*20 :: rhsfmtw, valfmtw
      character*72 :: title
      character*72 :: titlew
      integer*8 :: factors
      integer, pointer :: ZDiag(:) => null()
      integer, pointer :: NZDiag(:) => null()
      integer, pointer :: ZDiag_ind(:) => null()
      integer, pointer :: TempI(:) => null()
      real*8, pointer :: TempR(:) => null()
      integer, pointer :: PermL(:) => null()

!
!      call hbcode1(m, n, nnz, values, rowind, colptr)

!     ================================================================
!     ... READING A SPARSE MATRIX IN STANDARD FORMAT
!     ================================================================

!     ------------------------
!     ... read in header block
!     ------------------------
      
      read ( *, 1000 ) title , key   ,                         &
     &                 totcrd, ptrcrd, indcrd, valcrd, rhscrd, &
     &                 mxtype, m  , n  , nnz, neltvl,          &
     &                 ptrfmt, indfmt, valfmt, rhsfmt          
 1000 format ( a72, a8 / 5i14 / a3, 11x, 4i14 / 2a16, 2a20 )

!    --------------------------
!     ... read matrix structure
!    --------------------------
      allocate(colptr(n+1))
      allocate(rowind(nnz))
      read ( *, ptrfmt ) ( colptr (i), i = 1, n+1 )
      read ( *, indfmt ) ( rowind (i), i = 1, nnz )
      allocate(values(nnz))

      if  ( valcrd .gt. 0 )  then
!    ------------------------
!    ... read matrix values
!    ----------------------
          read ( *, valfmt ) ( values (i), i = 1, nnz )
      endif

      allocate(slu_PermL(m))
      allocate(slu_PermC(n))
      tol=1.D-14
      ldb = n
      allocate(l_b(n)) 



      print *, ''
      print *, '- Rectangular Matrix'
      print *, '- Size :',m,' x ',n
      print *, '- NNZ :',nnz

! =========================================================================
! FACTORISATION/GETTING/CLEANING LU OF THE RECTANGULAR MATRIX             !
! =========================================================================     
      iopt = 1
      nrhs=3

      

      call CPU_TIME( t1 )
      call c_fortran_dgssv(values,rowind,colptr,iopt,m,n,nnz,nrhs,  &
     &            slu_PermC,slu_PermL, l_b, ldb, factors, info )

  


      call CPU_TIME( t2 )
      print *,'- CPU Time recovery and factorization :',t2 - t1

      nbnz=nrhs 
      print *, '- copy of elements :',nbnz

!                         nblag = n / nbphys = m
!      allocate( slu_diagu(n))
      allocate( slu_vale(nbnz))
      allocate( slu_rowind(nbnz))
      allocate( slu_colptr(n+1))

! --------------------------------------------------------------------
!  recovery of the decomposed matrix in compressed column format 
! --------------------------------------------------------------------      
      call CPU_TIME( t1 )
      iopt = 4
      call c_fortran_dgssv(slu_vale,slu_rowind,slu_colptr,iopt,m,n,         &
     &    nnz,nrhs,slu_PermC,slu_PermL,l_b,ldb,factors,info)

! ------------------------------
! SuperLU Object cleaning
! ------------------------------
      
      iopt = 3
      call c_fortran_dgssv(slu_vale,slu_rowind,slu_colptr,iopt,m,n,        &
     &    nnz,nrhs,slu_PermC,slu_PermL,l_b,ldb,factors,info)
      call CPU_TIME( t2 )
      print *, '- Temps CPU Recopie / nettoyage ',t2 - t1
      
! =========================================================================
!       RECOVERY OF L1 AND L2 SUCH THAT PL=[L1 L2]                        !
! ========================================================================= 

! -----------------------------------------------------------
! recovery of the indices of the "null" elements of the diagonal
! -----------------------------------------------------------

      call CPU_TIME( t1 )
      allocate( ZDiag(m))
      allocate( NZDiag(m))
      allocate( ZDiag_ind(m))
      do i1=1,m 
        ZDiag(i1)=0 
        NZDiag(i1)=0 
        ZDiag_ind(i1)=0 
      end do

      nbz_diag=0
      nbnz_diag=0

      do i1=1,n
      if (abs(l_b(i1)) .lt. tol) then
        nbz_diag=nbz_diag+1
        ZDiag(i1)=nbz_diag
        ZDiag_ind(nbz_diag)=i1
      else
        nbnz_diag=nbnz_diag+1
        NZDiag(i1)=nbnz_diag
      endif
      end do

      do i1=n+1,m
        nbz_diag=nbz_diag+1
        ZDiag(i1)=nbz_diag
        ZDiag_ind(nbz_diag)=i1
      end do
      call CPU_TIME( t2 )
      print *, ''
      print *, '*** recovery of non null diagonal elements ***'
      print *, '- CPU Time ',t2 - t1

! --------------------------------------------------------------------
! Construction of the "reverse" permutation to get the right indexes
! --------------------------------------------------------------------

      call CPU_TIME( t1 )
      
      allocate(PermL(m))
      allocate(TempI(m))
      allocate(ibid1(m))
      allocate(ibid2(m))
    
      do i1=1, m
        PermL(i1)=i1
        TempI(i1)=slu_PermL(i1)
      end do  
      


!      call dsorti('LA', .true., m, TempI, PermL)
      

      igap=(m+1)/2
  70  continue
      if (igap .eq. 0) go to 9000
         do 90 i = igap, m
            j = i-igap
   80       continue
!
            if (j.lt.0) go to 90
!           
            if (TempI(j).gt.TempI(j+igap)) then
               temp = TempI(j)
               TempI(j) = TempI(j+igap)
               TempI(j+igap) = temp
               temp = PermL(j)
               PermL(j) = PermL(j+igap)
               PermL(j+igap) = temp

            else
               go to 90
            endif
            
            j = j-igap
            go to 80
   90    continue
         igap = igap / 2
         go to 70

9000    continue


      call CPU_TIME( t2 )
      print *, ''
      print *, '*** Construction of the inverse permutation ***'
      print *, '- CPU time -- construction of the inverse permutation ',t2 - t1
      
! --------------------------------------------------------------------
! L1 and L2 counting
! we put 1's on the diagonal
! --------------------------------------------------------------------
      
      ! do i1=1,n+1
      !   print *, i1, slu_colptr(i1)
      ! end do

      call CPU_TIME( t1 )
      nzL1=0
      nzL2=0
     
      do i1=1,nbnz_diag
         ! print *, '*****'
         ! print *, i1, NZDiag(i1)
        nnz=slu_colptr(NZDiag(i1)+1)-slu_colptr(NZDiag(i1))
        pos=slu_colptr(NZDiag(i1))
         ! print *, nnz, pos
        do j1=1,nnz
           ind_lig=slu_rowind(pos-1 + j1)
           if ( ZDiag( ind_lig ) .ne. 0 ) then
             nzL2=nzL2+1
           elseif ( ind_lig .ge. NZDiag(i1) ) then
             if ( ind_lig .eq. NZDiag(i1) ) then
               slu_vale(pos-1+j1)=1.D0
             endif
             nzL1=nzL1+1
           endif
        end do
      
      end do

      call CPU_TIME( t2 )
      print *, ''
      print *, '*** nnz counting for L1 et L2 ***'
      print *, '- CPU time -- nnz L1/L2 ',t2 - t1
      print *, '-  Size of  L1',nbnz_diag,' x ',nbnz_diag
      print *, '-  NNZ :',nzL1



      call CPU_TIME( t1 )
      allocate(TempR(m))
      if ( nzL2 .eq. 0) then

! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                   !    
!                CASE 1  :  nzL2 == 0                               !
!                                                                   !  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!                majority case => 
!                The constraints to be eliminated are of type DX = 0
!                and constraints are not redundant nzL2 == 0
!                the core base is a truncated identity, 
!                 no further calculation required
! --------------------------------------------------------------------
         print *, '- Size of  L2 :',nbz_diag,' x ',nbnz_diag
         print *, '- NNZ :',0
! --------------------------------------------------------------------
!   Storing the truncated identity 
! --------------------------------------------------------------------


         allocate(lp_vale(nbz_diag))
         allocate(lp_lign(nbz_diag))
         allocate(lp_long(nbz_diag+1))
 



         nnzT=0
         do ind_col=1, nbz_diag 

             nbnz=0
             lp_long(ind_col)%lrow=nbnz+1

             nnzT=nnzT+nbnz+1
! --------------------------------------------------------------------
!             call codent(ind_col, 'D0', k8bid)
! --------------------------------------------------------------------
             allocate(lp_vale(ind_col)%prow (nbnz+1)) 
             allocate(lp_lign(ind_col)%qrow (nbnz+1)) 

             lp_vale(ind_col)%prow(nbnz+1)=1.D0
             lp_lign(ind_col)%qrow(nbnz+1)=PermL( nbnz_diag + ind_col )
                     
         end do 

      elseif ( nzL2 .gt. 0) then

! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                   !    
!                CASE 2  :  nzL2 neq 0                               !
!                                                                   !  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         print *, ''
         print *, '- Taille de L2 :',nbz_diag,' x ',nbnz_diag
         print *, '- Nombre de termes nons nuls :',nzL2

         allocate(L1_vale(nzL1))
         allocate(L1_rowind(nzL1))
         allocate(L1_colptr(n+1))

         allocate(L2_vale(nzL2))
         allocate(L2_rowind(nzL2))
         allocate(L2_colptr(n+1))

! --------------------------------------------------------------------      
! Filling L1 et L2
! --------------------------------------------------------------------

         nzL1=0
         nzL2=0
         L1_colptr(1)=1
         L2_colptr(1)=1
         do i1=1,nbnz_diag
             nnz=slu_colptr(NZDiag(i1)+1)-slu_colptr(NZDiag(i1))
             pos=slu_colptr(NZDiag(i1))
              do j1=1,nnz
                ind_lig=slu_rowind(pos-1 + j1)
                if ( ZDiag( ind_lig ) .ne. 0 ) then
                  L2_vale(nzL2+1)=slu_vale(pos-1 + j1)
                  L2_rowind(nzL2+1)=ind_lig
                  nzL2=nzL2+1
                elseif ( ind_lig .ge. NZDiag(i1) ) then
                  L1_vale(nzL1+1)=slu_vale(pos-1 + j1)
                  L1_rowind(nzL1+1)=ind_lig
                  nzL1=nzL1+1
                endif

             end do
             L1_colptr(i1+1)=nzL1+1
             L2_colptr(i1+1)=nzL2+1
         end do




         print *, ''
         print *, '*** Filling of L1 and L2 ***'

! =========================================================================
! Kernel construction (Factorisation of L1 / Transposition of L2 ...)           !
! ========================================================================= 


! ---------------------------
!   ---  factorization of L1
! ---------------------------
         print *, ''
         print *, '*** Factorization of L1 ***'
         iopt=1
         nrhs=3

         call c_fortran_dgssv(L1_vale,L1_rowind,L1_colptr,           &
     &       iopt,nbnz_diag,nbnz_diag,nzL1,nrhs,                     &
     &   slu_PermC,slu_PermL, l_b, nbnz_diag, factors, info)


! ------------------------------------
!   --- construction of Transpose(L2)
! ------------------------------------
     !     iopt=5
     !     nrhs=1
 
     ! !     call c_fortran_dgssv(L2_vale,L2_rowind,L2_colptr,          &
     ! ! &                   iopt,m,nbnz_diag,nzL2,nrhs,                &
     ! ! &      slu_PermC,slu_PermL, l_b, nbnz_diag, factors, info )




         allocate(tval(nzL2))
         allocate(trow(nzL2))
         allocate(tcolp(m+1))    


!  Initialize TCOLP

         do 100 i=1,m
             tcolp(i) = 0
         100 continue



!         store length of each row (column in the transpose) in tcolp
            do 200 i=1,nzL2
             tcolp(L2_rowind(i)) = tcolp(L2_rowind(i))+1
         200 continue



!        tcolp points to position after end of columns
            tcolp(1) = tcolp(1)+1
            do 300 i=2,m
              tcolp(i) = tcolp(i-1)+tcolp(i)
         300 continue



!        determine trow and tval (if yesa = .true.)
            if (.true.) then
              do 500 i=nbnz_diag,1,-1
                j1 = L2_colptr(i)
                j2 = L2_colptr(i+1) - 1
                do 400 j = j1,j2
                  x = L2_rowind(j)
                  jpos = tcolp(x)-1
                  trow(jpos) = i
                  tval(jpos) = L2_vale(j)
                  tcolp(x) = jpos
         400     continue
         500   continue
            else
              do 700 i=nbnz_diag,1,-1
                j1 = L2_colptr(i)
                j2 = L2_colptr(i+1) - 1
                do 600 j = j1,j2
                  x = L2_rowind(j)
                  jpos = tcolp(x)-1
                  trow(jpos) = i
                  tcolp(x) = jpos
         600     continue
         700   continue
            endif

!  set value of tcolp(m+1)
         tcolp(m+1) = nzL2+1




         call CPU_TIME( t2 )
         print *, ''
         print *, '***Filling /Factorization of L1 /Transposition of L2***'
         print *, '- Temps CPU ',t2 - t1

!        boucle  sur les 2nds membres   

   
         call CPU_TIME(t1)

         allocate(lp_vale(nbz_diag))
         allocate(lp_lign(nbz_diag))
         allocate(lp_long1(nbz_diag+1))


         nnzT=0

         do ind_col=1, nbz_diag
             i1=ZDiag_ind(ind_col)

!         Right hand side recovery
             
             nnz=tcolp(i1+1)-tcolp(i1)
             pos=tcolp(i1)


             do j1=1,nnz
               ind_lig=trow(pos-1 + j1)
               l_b(ind_lig)=-tval(pos-1 + j1)
             end do 
 

!         Resolution

             if ( nnz .gt. 0 ) then
!             !  signe "-" parce qu'on resoud avec la transposee 
               iopt=-2
               nrhs=1
               info=0


               call c_fortran_dgssv(L1_vale,L1_rowind,L1_colptr,         &
     &               iopt,nbnz_diag,nbnz_diag,nzL1,nrhs,                 &
     &          slu_PermC,slu_PermL, l_b, nbnz_diag, factors, info )



             endif



!        Storing

             nbnz=0
             if (nnz .gt. 0) then
!        nnz of the solution
                do k1=1,nbnz_diag
                  
                  if (abs(l_b(k1)) .gt. tol) then 
                    nbnz=nbnz+1
                    TempI(nbnz)=k1
                    TempR(nbnz)=l_b(k1)
                  endif
                  l_b(k1)=0.D0
                end do
             endif

            
             lp_long1(ind_col)=nbnz+1

             nnzT=nnzT+nbnz+1
!             call codent(ind_col, 'D0', k8bid)

             allocate(lp_vale(ind_col)%prow(nbnz+1)) 
             allocate(lp_lign(ind_col)%qrow(nbnz+1)) 
  


!         Filling (-L1**-T * L2**T) and reinitialization of the second right hand side
             if (nnz .gt. 0) then
               do l1=1,nbnz
                 k1=TempI(l1)
                 
                 lp_vale(ind_col)%prow(l1)=TempR(l1)
                 lp_lign(ind_col)%qrow(l1)=PermL (k1)

               end do
             endif
          
!         Filling with the identity matrix
      
             lp_vale(ind_col)%prow(nbnz+1)=1.D0
             lp_lign(ind_col)%qrow(nbnz+1)=PermL( nbnz_diag + ind_col )


         end do 

!        Cleaning SuperLU objects
         iopt=3

         call c_fortran_dgssv(L1_vale,L1_rowind,L1_colptr,             &
     &               iopt,nbnz_diag,nbnz_diag,nzL1,nrhs,               &
     &      slu_PermC,slu_PermL, l_b, nbnz_diag, factors, info )


! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                   !    
!               END OF BOTH CASES                                   !
!                                                                   !  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      endif 




      call CPU_TIME(t2)
      print *, ''
      print *, '*** End of the resolution and storage ***'
      print *, '- Temps CPU ',t2 - t1


! =========================================================================
!   Storing the null basis T                           !
! =========================================================================
      print *, ''
      print *, '*** Storing the null basis T ***'
      print *, '- nnz of T :',nnzT
      print *, '- number of columns of the basis T :',nbz_diag

      call CPU_TIME(t1)

      allocate(l_Tvalues(nnzT))
      allocate(l_Trowind(nnzT))
      allocate(l_Tcolptr(nbz_diag+1)) 

      l_Tcolptr(1)=1
      nnzT=0
      
      do i1=1,nbz_diag
        nnz= lp_long1(i1)  
        do j1=1,nnz
          l_Tvalues(nnzT + j1)=lp_vale(i1)%prow(j1)
          l_Trowind(nnzT + j1)=lp_lign(i1)%qrow(j1)
        end do
        l_Tcolptr(i1+1)=l_Tcolptr(i1)+nnz
        nnzT=nnzT+nnz
!        call codent(i1, 'D0', k8bid)
      end do


! =========================================================================      
!    Get the nullbasis with a Matlab Format I,J,Value              !
! =========================================================================
  open(unit=1,file="matZ",form="formatted",status="replace",action="write")
  write(unit=1,fmt=*),m,nbz_diag,nnzT
   do i1=1,nbz_diag
     nnz=l_Tcolptr(i1+1)-l_Tcolptr(i1)
     posK=l_Tcolptr(i1)
     do j1=1,nnz
       write(unit=1,fmt=*),l_Trowind(posK-1 + j1),i1,l_Tvalues(posK-1 + j1)
     end do
     l_Tcolptr(i1+1)=l_Tcolptr(i1)+nnz
   end do
   
close (unit=1)
! =========================================================================
!   Get the nullbasis with a Harwell-Boeing Format   (uncomment)          !
! =========================================================================

      
!       ptrcrdw = ( nbz_diag - 1 ) / 8 + 1
!       indcrdw = ( nnzT - 1 ) / 8 + 1
!       valcrdw = ( nnzT - 1 ) / 5 + 1
!       rhscrdw = 0
!       totcrdw = ptrcrdw + indcrdw + valcrdw + rhscrdw  
!       mxtypew = 'RUA'
!       neltvlw = 0
!       ptrfmtw = '(8I10)'
!       indfmtw = '(8I10)'

!       valfmtw = '(5E14.7)'
!       rhsfmtw = '(10F7.1)'

!       titlew= 'constraint matrix'
!       keyw='key'

!       ================================================================
! !     ... READING A SPARSE MATRIX IN STANDARD FORMAT
! !     ================================================================

! !     ------------------------
! !     ... read in header block
! !     ------------------------
      
!       write ( *, 1000 ) titlew , keyw   ,                         &
!      &                 totcrdw, ptrcrdw, indcrdw, valcrdw, rhscrdw, &
!      &                 mxtypew, m  , n  , nnz, neltvlw,          &
!      &                 ptrfmtw, indfmtw, valfmtw, rhsfmtw          
!  ! 1000 format ( a72, a8 / 5i14 / a3, 11x, 4i14 / 2a16, 2a20 )

! !    --------------------------
! !     ... read matrix structure
! !    --------------------------

!       write ( *, ptrfmt ) ( colptr (i), i = 1, n+1 )
!       write ( *, indfmt ) ( rowind (i), i = 1, nnz )
      
! !    ------------------------
! !    ... read matrix values
! !    ----------------------
!       write ( *, valfmt ) ( values (i), i = 1, nnz )






!       print '(a72,a8)' ,       titlew, keyw
!       open (13, file = "nullspace", status='replace', action='write')

!       write (13, 5) titlew, keyw
! 5     format (a72,a8)
!       write ( 13, 10 ) totcrdw, ptrcrdw, indcrdw, valcrdw, rhscrdw
! 10     format (5i14)
!       write ( 13, 15 ) mxtypew, nnzT, nbz_diag,  nnzT, neltvlw
! 15     format (a3,11x,4i14)     
!       write ( 13, 20) ptrfmtw, indfmtw, valfmtw, rhsfmtw
! 20     format (2a16,2a20)      
      
!       seuil=(nbz_diag+1)/8
!       reste=(nbz_diag+1)-(seuil*8)
!       inc=0
!       do i=1,seuil
!         write (13, 30)  (l_Tcolptr(i1), i1 =inc+1,inc+8 )
!         inc=inc+8
!       end do
!       write(13,30) (l_Tcolptr(i1), i1 =inc+1,inc+reste)
! 30     format (8I10) 
!       seuil=nnzT/8
!       reste=nnzT-(seuil*8)
!       inc=0
!       do i=1,seuil
!         write ( 13,30)  (l_Trowind(i1), i1 =inc+1,inc+8 )
!         inc=inc+8
!       end do
!       write(13,30) (l_Trowind(i1), i1 =inc+1,inc+reste)

!       seuil=nnzT/5
!       reste=nnzT-(seuil*5)
!       inc=0
!       do i=1,seuil
!         write ( 13, 40)  (l_Tvalues(i1), i1 =inc+1,inc+5 )
!         inc=inc+5
!       end do
!       write(13,40) (l_Tvalues(i1), i1 =inc+1,inc+reste)
! 40     format (5E14.7) 
!       close (13) 
      call CPU_TIME(t2)
      print *, '- Temps CPU ',t2 - t1

      

           stop
      end
      