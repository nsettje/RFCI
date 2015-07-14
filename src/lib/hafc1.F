      Subroutine Hafc1(ca,bout,a,buf,x,scr,itop,ct,Nbf,Nmo,Nocc,Nx)
c        *****  Copyright (c)  1993  Polyatomics Research Inst. *****
c        *****  A. Komornicki: Polyatomics Research Institute   ***** 
c        *****  A Hessian integral transformation.              *****
c        *****  Routine reads sorted two(2) electron integrals  ***** 
c        *****  and does the first half of the transformation.  ***** 
c        *****  Sort in this version is done in core.           *****
c        *****  Calls:  Xrewind                                 ***** 
c        *****          Xread                                   ***** 
c        *****          Gather                                  ***** 
c        *****          Sgemm                                   ***** 
c        *****     Ca    - Molecular orbitals                   *****
c        *****     Bout  - Output half transformed integrals    *****
c        *****     A     - Output integrals                     *****
c        *****     Buf   - Input  integrals                     *****
c        *****     X     - Input integrals in square form       *****
c        *****     Scr   - Scratch array (Nbf,Nocc)             *****
c        *****     Itop  - Scatter pointer                      *****
c        *****     Ct    - Transpose of the Ca matrix.          *****
c        *****     Nbf   - Nbasis                               *****
c        *****     Nmo   - Number of vectors                    *****
c        *****     Nocc  - Number of occupied vectors           *****
      Parameter (Zero=0.0E0, One=1.0e0)
      Common/TranIO/Iuao,iusrt,Iumo,Iw
      Dimension ca(Nbf,Nmo),a(Nmo,Nocc),buf(Nx),x(*)
      Dimension bout(Nx,*),scr(*),itop(*),ct(Nmo,Nbf)
#ifndef  cray
#define  Sgemm  Dgemm
#endif
      Ia(kk) = kk*(kk-1)/2
c        *****  Set up scatter pointers:                        *****
      Nbasis= Nbf
      Nddx  = Nbasis*Nbasis
      ij = 0
      Do 120 i=1,Nbasis
      Do 110 j=1,Nbasis
      ij = ij+1
  110 itop(ij) = Ia(max(i,j))+min(i,j)
  120 continue
      Do 140 i=1,Nbasis
      Do 140 j=1,Nmo
  140 ct(j,i) = ca(i,j)
c
      Call Xrewind(Iuao)
c        *****  Read in all kl records of AO integrals.         *****
      Do 300 kl=1,Nx
      Call Xread(Iuao,buf,Nx,ierr)
c        *****  Do the work.                                    *****
c        *****  Put integrals in square form:transform          *****
      Call Gather(Nddx,x,buf,itop)
      Call Sgemm('n','n',nbasis,nocc,nbasis,One,x,nbasis,ca,nbasis,
     1                   Zero,scr,nbasis)
      Call Sgemm('n','n',Nmo,nocc,nbasis,One,ct,Nmo,scr,nbasis,
     1                   Zero,a,Nmo)
c        *****  Store half transformed integrals                *****
      ij = 0
      Do 220 i=1,Nmo
      jmax = min(i,Nocc)
      Do 210 j=1,jmax
      ij = ij +1
  210 bout(kl,ij) = a(i,j)
  220 continue
  300 continue
      Return
      End
