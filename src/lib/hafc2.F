      Subroutine Hafc2(Vec,xbin,a,x,scr,itop,ct,nbf,Nmo,Nocc,Ncore,nx)
      Character*16 S_nrec,S_Ntot,S_nints
c        *****  Copyright (c)  1993  Polyatomics Research Inst. *****
c        *****  A. Komornicki: Polyatomics Research Institute   ***** 
c        *****  Hessian integral transformation.                *****
c        *****  Routine drives second half of transformation    *****
c        *****  Calls:  Gather                                  *****
c        *****          Dgemm()                                 *****
c        *****          Xwrite                                  *****
c        *****          Intcoma                                 *****
c        *****          Xerror                                  *****
c        *****     Vec   - Molecular orbitals                   *****
c        *****     Xbin  - Input half transformed integrals     *****
c        *****     A     - Output (MO) integrals.               *****
c        *****     X     - Input integrals in square form       *****
c        *****     Scr   - Scratch array (Nbf,Nmo)              *****
c        *****     Itop  - Scatter pointer                      *****
c        *****     Ct    - Transpose of the Ca matrix.          *****
c        *****     Nbf   - Nbasis                               *****
c        *****     Nmo   - Number of vectors                    *****
c        *****     Nocc  - Number of occ vectors in trans.      *****
c        *****     Ncore - Number of frozen orbitals.           *****
      Parameter (Zero=0.0e0, One=1.0e0, Tol=4.0e-9)
      Parameter (Maxb=512*96, Maxa=Maxb/2-1)
      Common/Bufmo/ Ncnt,jxx(Maxa),xxx(Maxa),xxdum
      Common/TranIO/Iuao,iusrt,Iumo,Iw
      Dimension Vec(nbf,Nmo),a(Nmo,Nmo),xbin(nx,*),x(*),itop(*),scr(*)
      Dimension ct(Nmo,Nbf)
 1101 Format(10x,'End of Hessian integral transformation ',/,
     &       10x,'Number of records written        ', A,/,
     &       10x,'Number of non-zero integrals   ', A,/,
     &       10x,'Maximum number of integrals    ', A,/ )
 1102 Format(12x,'Actual number and maximum for each type:   ',/,
     &       10x,'(oo//oo)  (vo//oo)  (vv//oo)  (vo//vo)',7x,'Total',/,
     &        9x,4(i8,2x),4x,A,/,9X,4(i8,2x),4x,A, /)
 1103 Format(10x,3(2H*-),'Error in Hessian transformation ',3(2H-*),/
     &       10x,'There appear to be no virtual orbitals - Exit ',/)
c
      Ia(kk) = kk*(kk-1)/2
c        *****  Set up pointer index.                           *****
      noooo = 0
      nvooo = 0
      nvvoo = 0
      nvovo = 0
      Ncnt  = 0
      nrec  = 0
      nvirt = Nmo-Nocc
      Nocc1 = Nocc+1
c        *****  Set up scatter pointers:                        *****
      Nbasis = nbf
      Nddx   = Nbasis*Nbasis
      ij = 0
      Do 120 i=1,Nbasis
      Do 110 j=1,Nbasis
      ij = ij+1
  110 itop(ij) = Ia(max(i,j))+min(i,j)
  120 continue
      Do 140 i=1,Nbasis
      Do 140 j=1,Nmo
  140 ct(j,i) = Vec(i,j)
c        *****  Start the transformation.                       *****
      ij = 0
      Do 500 i=1,Nmo
      jend = min(Nocc,i)
      Do 500 j=1,jend
      ij = ij+1
c        *****  The following loop replaces the sort step:      *****
c        *****  Put integrals in square form:transform          *****
      Call Gather(nddx,x,xbin(1,ij),itop)
      idx  = i+Ncore
      jdx  = j+Ncore
      labij= Ishft(idx,16) +jdx
c        *****  Check for two possible cases.                   *****
      If(i.le.Nocc) then
c        *****  A(Nmo,Nmo) = Vec(n,mo)*x(n,n)*Vec(n,mo)           *****
         Call Dgemm('n','n',Nbasis,Nmo,Nbasis,One,x,Nbasis,Vec,Nbasis,
     &                      Zero,scr,Nbasis)
         Call Dgemm('n','n',Nmo,Nmo,Nbasis,One,ct,Nmo,scr,Nbasis,
     &                      Zero,a,Nmo)
c        *****  Fill the output buffers.                        *****
c        *****  (oc,oc//oc,oc) block first.                     *****
         Do 300 k=1,i
         lmax = k
         If(k.eq.i) lmax=j
         Do 300 l=1,lmax
         val= a(k,l)
         If(abs(val).lt.tol) go to 300
         Ncnt = Ncnt +1
         noooo= noooo+1
         kdx  = k+Ncore
         ldx  = l+Ncore
         label= Ishft(labij,32)+Ishft(kdx,16)+ldx
         jxx(Ncnt) = label
         xxx(Ncnt) = val
         If(Ncnt.ge.Maxa) then
           Call Xwrite(Iumo,Ncnt,Maxb,ierr)
           Ncnt = 0
           nrec = nrec +1
         endif
  300    continue
c        *****  Now the   (oc,oc//virt,all)  block              *****
         Do 320 k=Nocc1,Nmo
         Do 320 l=1,k
         val= a(k,l)
         If(abs(val).lt.tol) go to 320
         Ncnt = Ncnt +1
         If(l.le.Nocc) then
           nvooo = nvooo+1
         else
           nvvoo = nvvoo+1
         endif
         kdx  = k+Ncore
         ldx  = l+Ncore
         label= Ishft(kdx,48)+Ishft(ldx,32)+labij
         jxx(Ncnt) = label
         xxx(Ncnt) = val
         If(Ncnt.ge.Maxa) then
           Call Xwrite(Iumo,Ncnt,Maxb,ierr)
           Ncnt = 0
           nrec = nrec+1
         endif
  320    continue
      else
c        *****  A(nvirt,j) = Vec(n,nvirt)*x(n,n)*Vec(n,j)       *****
         Call Dgemm('n','n',Nbasis,j,Nbasis,One,x,Nbasis,Vec,Nbasis,
     *                      Zero,scr,Nbasis)
         Call Dgemm('t','n',nvirt,j,Nbasis,One,Vec(1,nocc1),Nbasis,
     *                      scr,Nbasis,Zero,a,Nmo)
c        *****  Now fill the output buffers.                    *****
c        *****  (virt,oc//virt,oc) block                        *****
         Do 420 l=1,j
         kend = Nmo
         If(l.eq.j) kend=i
         Do 420 k=Nocc1,kend
         val= a(k-Nocc,l)
         If(abs(val).lt.tol) go to 420
         Ncnt = Ncnt +1
         nvovo= nvovo+1
         kdx  = k+Ncore
         ldx  = l+Ncore
         labkl= Ishft(kdx,16) +ldx
         If(labkl.gt.labij) then
            lab = Ishft(labkl,32) +labij
         else
            lab = Ishft(labij,32) +labkl
         endif
         jxx(Ncnt) = lab
         xxx(Ncnt) = val
         If(Ncnt.ge.Maxa) then
           Call Xwrite(Iumo,Ncnt,Maxb,ierr)
           Ncnt = 0
           nrec = nrec+1
         endif
  420    continue
      endif
  500 continue
c        *****  Final flush of output buffers.                  *****
      Call Xwrite(Iumo,Ncnt,Maxb,ierr)
      nrec = nrec+1
      Ntot = (nrec-1)*Maxa +Ncnt
      ncp  = Ia(Nocc+1)
      nvp  = Ia(nvirt+1)
      moooo= ncp*(ncp+1)/2
      mvooo= ncp*Nocc*nvirt
      mvovo= Nocc*nvirt*(Nocc*nvirt+1)/2
      mvvoo= ncp*nvp
      nints= moooo+ mvooo+ mvvoo+ mvovo
      Call Intcoma(nrec,S_nrec,ierr)
      Call Intcoma(Ntot,S_Ntot,ierr)
      Call Intcoma(nints,S_nints,ierr)
      Write(Iw,1101) S_nrec,S_Ntot,S_nints
      Ntot = noooo+ nvooo+ nvvoo+ nvovo
      Call Intcoma(Ntot,S_Ntot,ierr)
      Write(Iw,1102) noooo, nvooo, nvvoo, nvovo, S_Ntot,
     &               moooo, mvooo, mvvoo, mvovo, S_nints
c        *****  Check for error in transformation:              *****
      If(Ntot.eq.noooo) then
        Write(Iw,1103)
        Call Xerror('Hafc2')
      endif
      Return
      End
