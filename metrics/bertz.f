                           PROGRAM BERTZ                          
* (NonIsoMorphic Subgraphs, Dist EV) This program finds all connected sub-
* structures and the nonisomorphic connected substructures of a structure
* which may contain multiple bonds and/or heteroatoms. Further, it finds
* and counts all subgraphs NtB and the nonisomorphic subgraphs NsB according
* to Bertz; for nonconnected graphs (ensembles of molecules) also
*                          April 22, 2001
*               by Gerta Ruecker and Christoph Ruecker
      implicit integer (a-z)                                            
      dimension matrix(200,200),used(99),path0(0:200),
     &path(0:200),vert1(200),vert2(200),subgr(0:199),
     &avoid(200),num(99,0:200),colsum(99),
     &mult1(200),mult2(200),v1(200),v2(200)
      character*80 form1, form2, trivia
      real acuma, acumb, usdt
      real*8 bond(99,99),dist(99,99),akj,balab,d(99),aux,
     &adja(99,99),r(99),v(99,99),del1,deln,balab0,del10,deln0,
     &color(200),fv1(99),fv2(99),adja2(99,99)
      library 'c:\win32app\salford\demo\fortran\eispack.obj'
      call clock@(acuma)
      open(10,file='bertz.a')
*      open(12,status='old',file='bertz.b')
      open(13,status='append',file='bertza.res')
      open(14,file='bertz.c')
      open(15,file='bertz.d')
*      open(16,status='old',file='bertz.e')
      open(17,file='bertz.g')
      read(5,*) n
      do 2 i=1,n                                                        
      do 2 j=1,n                                                       
2     bond(i,j)=0                                                       
*                                                                       
* INPUT                                                                 
* Trivial name and formats:                                             
*                                                                       
      read (5, 20) trivia, form1, form2                                
20    format (a80/a80/a80)
*                                                                       
* Input of structure, chains:                                           
*                                                                       
30    read(5,*) i,j                                                    
      if(i.eq.0) goto 50                                                
      do 40 k=0,j-i-1                                                 
      bond(i+k,i+k+1)=1                                                 
40    bond(j-k,j-k-1)=1                                                 
      goto 30                                                           
*                                                                       
* Input of structure, bonded pairs:                                     
*                                                                       
50    read (5,*) i,j                                                   
      if (i.eq.0) goto 55                                           
      bond(i,j) = bond(i,j)+1                                           
      bond(j,i) = bond(j,i)+1                                           
      goto 50
*
* Input of structure, heteroelements:
*
55    read(5,*) i,col
      if(i.eq.0) goto 60
      bond(i,i)=col
      goto 55 
* 
* PROGRAM PART I: Finding and characterizing substructures
*
60    do 502 i=1,200
      do 502 j=1,200
502   matrix(i,j)=0
*
* Numbering and counting edges
*
      dim=0
      do 600 i=1,n-1
      do 600 j=i+1,n
      if (bond(i,j).eq.0) goto 600
      dim=dim+1
      vert1(dim)=i
      vert2(dim)=j
      color(dim)=bond(i,j)
600   continue
      if (dim.le.200) goto 605
      write(6,*) 'Warning! Too many edges!'
      stop
*
* Construct the edge adjacency matrix
*
605   do 610 k=1,dim-1
      do 610 l=k+1,dim
      if (vert1(k).eq.vert1(l)) then
      matrix(k,l)=1
      matrix(l,k)=1
      else if ((vert2(k).eq.vert1(l)).or.(vert2(k).eq.vert2(l))) then
      matrix(k,l)=1
      matrix(l,k)=1
      end if
610   continue
*
* Print edge adjacency information
*
      write(6,611) trivia
      write(13,611) trivia
611   format(/a80)
      do 615 edge=1,dim
      write(6,612) edge,vert1(edge),vert2(edge),color(edge)
612   format(' Edge',i4,'   Vertices',i3,' -',i3,'  Edge type',f4.1)
615   continue
*      write(6,*) ' Edge adjacency matrix'
*      do 616 i=1,dim
*616   write(6,725)(matrix(i,j), j=1,dim)
*
* Program Part Path algorithm
* Search for adjacent edge systems
*
      nsub=0
      ssn=0
      subgr(0)=0
      do 618 i=1,n
      subgr(0)=subgr(0)+1
      ssn=ssn+1
618   write(10,1350) ssn,1,0,0.0,bond(i,i),0.0
      do 620 l=1,dim
      path(l)=0
      subgr(l)=0
620   avoid(l)=0
      do 710 i=1,dim
      step=1
      path(1)=i
      subgr(1)=subgr(1)+1
      ssn=ssn+1
      aux=0.1*(bond(vert1(i),vert1(i))+bond(vert2(i),vert2(i)))
      balab=color(i)*(2**(aux*(-0.5)))
      del1=(2**aux)/color(i)
      write(10,1350) ssn,2,1,balab,del1,del1,i
*
* Going forward from the given substructure to an adjacent bond
*
640   do 690 k=i+1,dim
*
* Is edge k eligible?
*
      do 650 l=1,step
      if(k.eq.path(l)) goto 690
650   continue
      if(avoid(k).gt.0) goto 690
      do 660 l=1,step
      if(matrix(path(l),k).eq.1) goto 670
660   continue
      goto 690
*
* Bond k is eligible, it is added to the given substructure,
* resulting in a new substructure
*
670   step=step+1
      path(step)=k
      subgr(step)=subgr(step)+1
      ssn=ssn+1
*
* Analysis of substructure found
*
      do 1100 ii=1,n
      used(ii)=0
      do 1100 jj=1,n
      adja(ii,jj)=0
      dist(ii,jj)=0
1100  continue
      do 1200 l=1,step
      adja(vert1(path(l)),vert2(path(l)))=color(path(l))
      adja(vert2(path(l)),vert1(path(l)))=color(path(l))
      used(vert1(path(l)))=1
      used(vert2(path(l)))=1 
1200  continue
      nn=0
      do 1202 j=1,n
      adja(j,j)=bond(j,j)
1202  nn=nn+used(j)
*
* Condense adjacency matrix into first nn columns/rows
*
      pointer=0
      do 1208 j=1,n
      if(used(j).eq.0) goto 1208
      pointer=pointer+1
*
* Shift column j left to column "pointer"
* 
      do 1204  ii=1,n
1204  adja(ii,pointer)=adja(ii,j)
*
* Shift row j up to row "pointer"
*
      do 1206 jj=1,n
1206  adja(pointer,jj)=adja(j,jj)
1208  continue
*
* Distance matrix for given substructure (Mueller, Knop, Szymanski, Trinajstic)
*
      do 1220 jj=1,nn
      do 1210 ii=1,nn
      if(adja(ii,jj).eq.0) then
      dist(ii,jj)=nn
      else
      dist(ii,jj)=adja(ii,jj)**(-1)
      end if      
1210  continue
1220  dist(jj,jj)=0
      l=1
1230  do 1250 jj=1,nn
      do 1250 ii=1,nn
      do 1240 kk=1,nn
      akj=dist(kk,ii)+dist(ii,jj)
1240  if(dist(kk,jj).gt.akj) dist(kk,jj)=akj
1250  continue
      l=l+l
      if(l.lt.nn-1) goto 1230
*
* Distance matrix of substructure is now complete
*
*      write(6,*) ' Distance matrix'
*      do 1255 ii=1,nn
*1255  write(6,1256) (dist(ii,jj), jj=1,nn)
*1256  format(10f6.3)
*
* Balaban's J for the given substructure
*
      balab=0
      do 1330 ii=1,nn
      d(ii)=0
      do 1320 jj=1,nn
1320  d(ii)=d(ii)+dist(ii,jj)
      d(ii)=d(ii)*2**(0.1*adja(ii,ii))
1330  continue
      do 1340 ii=1,nn-1
      do 1340 jj=ii+1,nn
      if(adja(ii,jj).eq.0) goto 1340
      balab=balab+(d(ii)*d(jj))**(-0.5)
1340  continue
      c=step-(nn-1)
      balab=step*balab/(c+1)
*
* Enriching the distance matrix with heteroatom information
*
      do 1342 ii=1,nn
      do 1342 jj=1,nn
      dist(ii,jj)=dist(ii,jj)*2**(0.1*adja(ii,ii))
1342  continue
      do 1345 jj=1,nn
      do 1345 ii=1,nn
      dist(ii,jj)=dist(ii,jj)*2**(0.1*adja(jj,jj))
1345  continue
*
* Eigenvalues of the heteroenriched distance matrix of the given substructure
*
      call rs(99,nn,dist,r,0,v,fv1,fv2,info)
      if(info.eq.0) goto 1347
      write(6,*) info,' Error!'
      stop
*
* Writing all substructure information in a file
*
1347  write(10,1350) ssn,nn,step,balab,r(nn),
     &-r(1),(path(l),l=1,step)
1350  format(i7, 2i3, f11.8, f12.7, f12.7, 28i3)
      goto 640
690   continue
*
* End of path reached, going back
*
      avoid(path(step))=step
      do 700 l=1,dim
      if(avoid(l).gt.step) avoid(l)=0
700   continue
      path(step)=0
      step=step-1
      if(step.eq.0) goto 710
      goto 640
*
* Path algorithm for this edge completed
*
710   continue
*
* Output of total number of substructures
*
      write(6,720)
720   format(' Number of connected substructures found, ordered by numbe
     &r of bonds')
      write(6,725) (l,l=0,dim)
      write(6,725) (subgr(l),l=0,dim)
725   format(28i6)
      do 730 l=0,dim
730   nsub=nsub+subgr(l)
      write(6,*) ' Total number of connected substructures is', nsub
      write(13,*) 'Total number of connected substructures is', nsub
*
* PROGRAM PART II: Sorting the substructures found
*
      endfile 10
      call cissue('sort /+9 <bertz.a >bertz.b',ifail)
      if(ifail.ne.0) write(6,*)'sort failed, ifail=',ifail
*
* PROGRAM PART III: Diffline reads the sorted file, removes duplicates
*
      ssn=0
      ssn0=0
      number=0
      do 805 i=1,99
      do 805 j=0,200
805   num(i,j)=0
      n0=0
      m0=0
      balab0=0
      del10=0
      deln0=0
      occursum=0
      occur=0
      open(12,status='old',file='bertz.b')
810   read(12,1350,end=840) ssn,n1,m,balab,del1,deln,(path(l),l=1,m)
      if(deln.ne.deln0) goto 830
      if(del1.ne.del10) goto 830
      if(balab.ne.balab0) goto 830
      if(m.ne.m0) goto 830
      if(n1.ne.n0) goto 830
      occur=occur+1
      goto 810
*
* New nonisomorphic substructure read from the file
*
830   if(ssn0.ne.0) then
      write(6,831)occur,ssn0,n0,m0,balab0,del10,deln0,(path0(l),l=1,m0)
      write(14,831)occur,ssn0,n0,m0,balab0,del10,deln0,(path0(l),l=1,m0)
831   format(2i7, 2i3, f11.8, f12.7, f12.7, 28i3)
      endif
      occursum=occursum+occur
*
* Setting back the parameters
*
      occur=1
      number=number+1
      num(n1,m)=num(n1,m)+1
      ssn0=ssn
      n0=n1
      m0=m
      balab0=balab
      del10=del1
      deln0=deln
      do 835 l=1,m0
835   path0(l)=path(l)
      goto 810
*
* All substructures are read now, the last one is written
*
840   write(6,831)occur,ssn,n1,m,balab,del1,deln,(path(l),l=1,m)
      write(14,831)occur,ssn,n1,m,balab,del1,deln,(path(l),l=1,m)
      occursum=occursum+occur
*
* Output of number of nonisomorphic substructures
*
      write(6,*) 'Number of nonisomorphic substructures is', number
      write(13,*) '  Number of nonisomorphic substructures is', number
      write(6,*) '         Number of vertices'
      write(6,842)(i,i=1,12),'  rsum'
842   format(8x,12i6,a6)
      write(6,*) 'Number'
      write(6,*) 'of edges'
      do 846 j=0,dim
      rowsum=0
      do 843 i=1,n
843   rowsum=rowsum+num(i,j)
      write(6,844) j, (num(i,j),i=1,n), rowsum
844   format(i8,28i6)
846   continue
      do 848 i=1,n
      colsum(i)=0
      do 848 j=0,dim
      colsum(i)=colsum(i)+num(i,j)
848   continue
      write(6,849) ' sum',(colsum(i),i=1,n)
849   format(/4x,a4,28i6)
      if(occursum.ne.nsub) then
      write(6,*)'Discrepancy: occursum=',occursum
      write(13,*)'Discrepancy: occursum=',occursum
      endif
*
* PROGRAM PART IV: Calculation of total number of connected subgraphs
* (Bertz, ntb), and of number of nonisomorphic subgraphs (Bertz, nsb)
*
      rewind 14
      ntb=0
      nsb=0
850   read(14,831,end=1730)occur,ssn,nn,mm,balab,del1,deln,(path(l),l=1,
     &mm)
      write(15,831)occur,ssn,nn,mm,balab,del1,deln,(path(l),l=1,mm)
      n2=0
      n3=0
      do 860 l=1,mm
      if(int(color(path(l))).eq.2)n2=n2+1
      if(int(color(path(l))).eq.3)n3=n3+1
860   continue
      inntb=occur*(3**n2)*(7**n3)
      ntb=ntb+inntb
      if((n2.eq.0).and.(n3.eq.0))then
      nsb=nsb+1
      goto 850
      else
*
* Regeneration of adjacency matrix of current substructure
*
      do 1400 ii=1,n
      used(ii)=0
      do 1400 jj=1,n
      adja(ii,jj)=0
      dist(ii,jj)=0
1400  continue
      do 1450 l=1,mm
      adja(vert1(path(l)),vert2(path(l)))=color(path(l))
      adja(vert2(path(l)),vert1(path(l)))=color(path(l))
      used(vert1(path(l)))=1
      used(vert2(path(l)))=1
1450  continue
      do 1452 j=1,n
1452  adja(j,j)=bond(j,j)
*
* Condense adjacency matrix into first nn columns/rows
*
      pointer=0
      do 1458 j=1,n
      if(used(j).eq.0) goto 1458
      pointer=pointer+1
*
* Shift column j left to column "pointer"
*
      do 1454 ii=1,n
1454  adja(ii,pointer)=adja(ii,j)
*
* Shift row j up to row "pointer"
*
      do 1456 jj=1,n
1456  adja(pointer,jj)=adja(j,jj)
1458  continue
*      write(6,*)' Adjacency Matrix'
*      do 1460 ii=1,nn
*1460  write(6,1256)(adja(ii,jj), jj=1,nn)
*
* Preparations for construction of subgraphs of current substructure
*
      do 1470 i=1,nn
      do 1470 j=1,nn
1470  adja2(i,j)=adja(i,j)
      do 1475 l=1,dim
      mult1(l)=0
      mult2(l)=0
      v1(l)=0
1475  v2(l)=0
      step=0
      do 1480 i=1,nn-1
      do 1480 j=i+1,nn
      if(int(adja(i,j)).le.1) goto 1480
      step=step+1
      mult1(step)=int(adja(i,j))
      mult2(step)=int(adja(i,j))
      v1(step)=i
      v2(step)=j
1480  continue
      if(step.ne.n2+n3) then
      write(6,*) 'Error! step=',step
      endif
*
* Construction of all subgraphs of current substructure
*
      step=step-1
      first=1
*
* Forward
*
1500  step=n2+n3
*
* Analysis of current subgraph
*
      k1=mult1(step)-first
      do 1700 k=k1,1,-1
      mult2(step)=k
      do 1520 st=1,step
      adja2(v1(st),v2(st))=mult2(st)
1520  adja2(v2(st),v1(st))=mult2(st)
*      write(15,*)'Adjacency matrix of less unsat. subgraph'
*      do 1530 ii=1,nn
*1530  write(15,1256)(adja2(ii,jj),jj=1,nn)
*
* Distance matrix for current subgraph (Mueller, Knop, Szymanski, Trinajstic)
*
      do 1536 jj=1,nn
      do 1533 ii=1,nn
      if(adja2(ii,jj).eq.0) then
      dist(ii,jj)=nn
      else
      dist(ii,jj)=adja2(ii,jj)**(-1)
      end if
1533  continue
1536  dist(jj,jj)=0
      l=1
1538  do 1550 jj=1,nn
      do 1550 ii=1,nn
      do 1540 kk=1,nn
      akj=dist(kk,ii)+dist(ii,jj)
1540  if(dist(kk,jj).gt.akj) dist(kk,jj)=akj
1550  continue
      l=l+1
      if(l.lt.nn-1) goto 1538
*
* Distance matrix of subgraph is now complete
*
*      write(6,*) 'Distance matrix of subgraph'
*      do 1555 ii=1,nn
*1555  write(6,1256) (dist(ii,jj), jj=1,nn)
*
* Balaban's J for the current subgraph
*
      balab=0
      do 1630 ii=1,nn
      d(ii)=0
      do 1620 jj=1,nn
1620  d(ii)=d(ii)+dist(ii,jj)
      d(ii)=d(ii)*2**(0.1*adja2(ii,ii))
1630  continue
      do 1640 ii=1,nn-1
      do 1640 jj=ii+1,nn
      if(adja2(ii,jj).eq.0) goto 1640
      balab=balab+(d(ii)*d(jj))**(-0.5)
1640  continue
      c=mm-(nn-1)
      balab=mm*balab/(c+1)
*
* Enriching the distance matrix with heteroatom information
*
      do 1642 ii=1,nn
      do 1642 jj=1,nn
      dist(ii,jj)=dist(ii,jj)*2**(0.1*adja2(ii,ii))
1642  continue
      do 1645 jj=1,nn
      do 1645 ii=1,nn
      dist(ii,jj)=dist(ii,jj)*2**(0.1*adja2(jj,jj))
1645  continue
*
* Eigenvalues of the heteroenriched distance matrix of the current subgraph
*
      call rs(99,nn,dist,r,0,v,fv1,fv2,info)
      if(info.eq.0) goto 1647
      write(6,*) info,' Error!'
      stop
*
* Writing all subgraph information in a file
*
1647  factor=1
      do 1648 st=1,step
      if(mult2(st).lt.mult1(st)) factor=factor*mult1(st)
1648  continue
      goccur=factor*occur
      write(15,831) goccur,ssn,nn,mm,balab,r(nn),-r(1),(path(l),l=1,mm)
1700  continue
      first=0
*
* Backwards
*
1710  step=step-1
      if(step.eq.0) goto 1720
      if(mult2(step).eq.1) then
      mult2(step)=mult1(step)
      goto 1710
      else
      mult2(step)=mult2(step)-1
      goto 1500
      endif
*
* End of substructure analysis
*
1720  continue
      goto 850
      endif
*
* Sorting the subgraphs found
*
1730  endfile 15
      call cissue('sort /+16 <bertz.d >bertz.e',ifail)
      if(ifail.ne.0) write(6,*) 'sort2 failed, ifail=',ifail
*
* Diffline reads the sorted file, removes duplicates
*
      write(6,*)'Total number of connected subgraphs is',ntb
      write(13,*)'Total number of connected subgraphs is',ntb
      ssn=0
      ssn0=0
      number=0
*      do 1735 i=1,99
*      do 1735 j=0,200
*1735  num(i,j)=0
      n0=0
      m0=0
      balab0=0
      del10=0
      deln0=0
      goccursum=0
      goccur=0
      open(16,status='old',file='bertz.e')
1740  read(16,831,end=1780)soccur,ssn,n1,m,balab,del1,deln,(path(l),l=1,
     &m)
      if(deln.ne.deln0) goto 1750
      if(del1.ne.del10) goto 1750
      if(balab.ne.balab0) goto 1750
      if(m.ne.m0) goto 1750
      if(n1.ne.n0) goto 1750
      goccur=goccur+soccur
      goto 1740
*
* New nonisomorphic subgraph read from the file
*
1750  if(ssn0.ne.0) then
      write(6,831)goccur,ssn0,n0,m0,balab0,del10,deln0,(path0(l),l=1,m0)
      write(17,831)goccur,ssn0,n0,m0,balab0,del10,deln0,(path0(l),l=1,m0
     &)
      endif
      goccursum=goccursum+goccur
*
* Setting back the parameters
*
      goccur=soccur
      number=number+1
*      num(n1,m)=num(n1,m)+1
      ssn0=ssn
      n0=n1
      m0=m
      balab0=balab
      del10=del1
      deln0=deln
      do 1770 l=1,m0
1770  path0(l)=path(l)
      goto 1740
*
* All subgraphs are read now, the last one is written
*
1780  write(6,831)goccur,ssn,n1,m,balab,del1,deln,(path(l),l=1,m)
      write(17,831)goccur,ssn,n1,m,balab,del1,deln,(path(l),l=1,m)
      goccursum=goccursum+goccur
      if(goccursum.ne.ntb) then
      write(6,*)'Error: goccursum=',goccursum
      write(13,*)'Error: goccursum=',goccursum
      endif
*
* Final output
*
      nsb=number
      write(6,*)'Number of nonisomorphic subgraphs is',nsb
      write(13,*)'  ','Number of nonisomorphic subgraphs is',nsb
      call clock@(acumb)
      usdt=acumb-acuma
      write(6,1850) usdt
      write(13,1850) usdt
1850  format(f10.3, ' seconds elapsed time')
      stop                                                              
      end                                                               

